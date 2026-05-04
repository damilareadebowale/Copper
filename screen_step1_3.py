#!/usr/bin/env python3
import os, glob, csv
import argparse
import functools
from multiprocessing import Pool
from pathlib import Path
from Bio.PDB import PDBParser, Polypeptide
import numpy as np
from scipy.spatial import cKDTree

# -----------------------
# Benchmarks
# -----------------------
BENCH = {
    "O00244": {"name": "ATOX1", "pairs": [(12, 15)]},
    "O14618": {"name": "CCS",   "pairs": [(22, 25), (141, 227), (144, 227), (244, 246)]},
    "O43819": {"name": "SCO1",  "pairs": [(133, 137)]},
    "O75880": {"name": "SCO2",  "pairs": [(169, 173)]},
    "Q14061": {"name": "COX17", "pairs": [(23, 24), (26, 55), (36, 45)]},
    "Q49B96": {"name": "COX19", "pairs": []},
}

def uniprot_from_af_filename(fn: str) -> str:
    # AF-O00244-F1-model_v6.pdb -> O00244
    parts = Path(fn).name.split("-")
    return parts[1] if len(parts) > 2 and parts[0] == "AF" else Path(fn).stem

def extract_cys_sg(model):
    residues = [r for r in model.get_residues() if Polypeptide.is_aa(r)]
    cys = []
    for r in residues:
        if r.get_resname() == "CYS" and ("SG" in r):
            sg = r["SG"]
            cys.append({
                "chain": r.get_parent().id,
                "resnum": int(r.id[1]),
                "icode": (r.id[2] or "").strip(),
                "coord": sg.get_coord().astype(float),
            })
    return residues, cys

# Top-level for multiprocessing
def process_one(pdb_path, min_d, max_d):
    fn = os.path.basename(pdb_path)
    pid = uniprot_from_af_filename(fn)
    bench_name = BENCH.get(pid, {}).get("name", "")

    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure(pid, pdb_path)
        model = structure[0]
    except Exception:
        # Parsed failed; keep step2_row absent so parent doesn't count it.
        return None

    residues, cys = extract_cys_sg(model)
    L = len(residues)
    n_cys = len(cys)
    step2_row = (pid, bench_name, fn, L, n_cys)

    if n_cys < 2:
        bench_update = {"pid": pid, "ge2cys": False, "pair_hit": set()}
        return (step2_row, None, [], bench_update)

    n_pairs_total = n_cys * (n_cys - 1) // 2
    coords = np.array([x["coord"] for x in cys], dtype=float)
    tree = cKDTree(coords)
    cand = tree.query_pairs(r=max_d, output_type="set")

    kept = []
    pair_hit = set()

    for i, j in cand:
        d = float(np.linalg.norm(coords[i] - coords[j]))
        if min_d <= d <= max_d:
            kept.append((i, j, d))
            if pid in BENCH and BENCH[pid]["pairs"]:
                a = (cys[i]["resnum"], cys[j]["resnum"])
                b = (cys[j]["resnum"], cys[i]["resnum"])
                for (r1, r2) in BENCH[pid]["pairs"]:
                    if a == (r1, r2) or b == (r1, r2):
                        pair_hit.add((r1, r2))

    if not kept:
        bench_update = {"pid": pid, "ge2cys": True, "pair_hit": pair_hit}
        step3_row = (pid, n_pairs_total, 0)
        return (step2_row, step3_row, [], bench_update)

    kept.sort(key=lambda t: (t[2], t[0], t[1]))
    n_pairs_kept = len(kept)
    step3_row = (pid, n_pairs_total, n_pairs_kept)

    pair_rows = []
    for (i, j, d) in kept:
        pair_rows.append({
            "pdb_id": pid, "bench_name": bench_name, "filename": fn,
            "length": L, "n_cys": n_cys,
            "n_pairs_total": n_pairs_total, "n_pairs_kept": n_pairs_kept,
            "chain1": cys[i]["chain"], "res1": cys[i]["resnum"], "icode1": cys[i]["icode"],
            "chain2": cys[j]["chain"], "res2": cys[j]["resnum"], "icode2": cys[j]["icode"],
            "sg_sg_dist": d,
        })

    bench_update = {"pid": pid, "ge2cys": True, "pair_hit": pair_hit}
    return (step2_row, step3_row, pair_rows, bench_update)

def main(pdb_dir, out_dir, min_d, max_d, workers):
    os.makedirs(out_dir, exist_ok=True)

    # (1) Recursive discovery (safe default; works for flat dirs too)
    files = sorted(glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True))
    total_started = len(files)

    # (2) Benchmark presence should reflect "exists in library", independent of parsing
    present_uniprots = set(uniprot_from_af_filename(os.path.basename(p)) for p in files)
    bench_found = {k: {"seen": (k in present_uniprots), "ge2cys": False, "pair_hit": set()} for k in BENCH}

    step2_path = os.path.join(out_dir, "step2_counts.csv")
    step3_path = os.path.join(out_dir, "step3_counts.csv")
    pairs_path = os.path.join(out_dir, "step3_pairs_check.csv")

    step2_ge2 = 0
    step3_ge1pair = 0

    with open(step2_path, "w", newline="") as f2, \
         open(step3_path, "w", newline="") as f3, \
         open(pairs_path, "w", newline="") as fp:

        w2 = csv.writer(f2)
        w2.writerow(["pdb_id", "bench_name", "filename", "length", "n_cys"])

        w3 = csv.writer(f3)
        w3.writerow(["pdb_id", "n_pairs_total", "n_pairs_kept"])

        wp = csv.DictWriter(fp, fieldnames=[
            "pdb_id","bench_name","filename","length",
            "n_cys","n_pairs_total","n_pairs_kept",
            "chain1","res1","icode1","chain2","res2","icode2","sg_sg_dist"
        ])
        wp.writeheader()

        worker_func = functools.partial(process_one, min_d=min_d, max_d=max_d)

        with Pool(processes=workers) as pool:
            result_iter = pool.imap_unordered(worker_func, files, chunksize=25)

            for i, out in enumerate(result_iter):
                if out is None:
                    continue

                step2_row, step3_row, pair_rows, bench_update = out
                w2.writerow(step2_row)

                pid = bench_update["pid"]
                if pid in bench_found:
                    # seen is already determined by existence; don't overwrite it here
                    bench_found[pid]["ge2cys"] = bench_found[pid]["ge2cys"] or bench_update["ge2cys"]
                    bench_found[pid]["pair_hit"].update(bench_update["pair_hit"])

                if step2_row[-1] >= 2:
                    step2_ge2 += 1

                if step3_row is not None:
                    w3.writerow(step3_row)
                    if step3_row[2] >= 1:
                        step3_ge1pair += 1

                for r in pair_rows:
                    wp.writerow(r)

                if (i + 1) % 1000 == 0:
                    print(f"Processed {i+1}/{total_started} | >=2Cys={step2_ge2} | pairs={step3_ge1pair}")

    print(f"\nStep 1: Total proteins started: {total_started}")
    print(f"Step 2: Proteins with >=2 Cys: {step2_ge2}")
    print(f"Step 3: Proteins with valid pairs: {step3_ge1pair}")
    print(f"Saved:\n  {step2_path}\n  {step3_path}\n  {pairs_path}")

    print("\nBenchmark validation:")
    for pid, meta in BENCH.items():
        st = bench_found[pid]
        print(f"  {meta['name']:6s} {pid} | in_library={st['seen']} | >=2Cys={st['ge2cys']} | Pairs={sorted(list(st['pair_hit']))}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--out_dir", default="out_step1_3")
    ap.add_argument("--min_d", type=float, default=2.0)
    ap.add_argument("--max_d", type=float, default=5.5)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    main(args.pdb_dir, args.out_dir, args.min_d, args.max_d, args.workers)


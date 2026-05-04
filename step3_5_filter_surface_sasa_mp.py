#!/usr/bin/env python3
import os, csv, argparse, functools
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool

import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

def build_filename_to_path_map(pdb_dir: str):
    m = {}
    dup = 0
    for p in Path(pdb_dir).rglob("*.pdb"):
        name = p.name
        if name in m and m[name] != str(p):
            dup += 1
            # keep first; warn later via dup count
        else:
            m[name] = str(p)
    return m, dup

def sg_sasa_map_from_pdb(pdb_path: str, probe_radius=1.4, n_points=100):
    """
    Returns dict: (chain, resseq, icode) -> SG atom sasa (Å^2).
    Prefers altLoc blank over others if present.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("X", pdb_path)
    model = structure[0]

    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    sr.compute(model, level="A")  # populate atom.sasa

    rank = {" ": 0, "A": 1}
    best_rank = {}
    sg_map = {}

    for atom in model.get_atoms():
        if atom.get_name() != "SG":
            continue
        res = atom.get_parent()
        if res.get_resname() != "CYS":
            continue
        chain = res.get_parent().id
        resseq = int(res.id[1])
        icode = (res.id[2] or "").strip()
        altloc = atom.get_altloc() if hasattr(atom, "get_altloc") else " "
        r = rank.get(altloc, 2)

        key = (chain, resseq, icode)
        sasa = getattr(atom, "sasa", None)
        if sasa is None:
            continue

        if key not in sg_map or r < best_rank.get(key, 99):
            sg_map[key] = float(sasa)
            best_rank[key] = r

    return sg_map

def process_file(task, probe_radius, n_points, sasa_thr, require_both):
    filename, rows, pdb_path = task
    out_rows = []
    fail_rows = []

    try:
        sg_map = sg_sasa_map_from_pdb(pdb_path, probe_radius=probe_radius, n_points=n_points)
    except Exception as e:
        for r in rows[:1]:
            fail_rows.append({
                "pdb_id": r["pdb_id"], "filename": filename,
                "reason": "sasa_failed", "detail": str(e)[:250]
            })
        return out_rows, fail_rows, 0, 0

    kept = 0
    total = 0

    for r in rows:
        total += 1
        chain1 = r["chain1"]; chain2 = r["chain2"]
        res1 = int(r["res1"]); res2 = int(r["res2"])
        icode1 = (r.get("icode1","") or "").strip()
        icode2 = (r.get("icode2","") or "").strip()

        s1 = sg_map.get((chain1, res1, icode1), 0.0)
        s2 = sg_map.get((chain2, res2, icode2), 0.0)

        ok = (s1 >= sasa_thr and s2 >= sasa_thr) if require_both else (s1 >= sasa_thr or s2 >= sasa_thr)
        if not ok:
            continue

        rr = dict(r)
        rr["res1_sg_sasa"] = f"{s1:.3f}"
        rr["res2_sg_sasa"] = f"{s2:.3f}"
        rr["pair_best_sg_sasa"] = f"{max(s1,s2):.3f}"
        out_rows.append(rr)
        kept += 1

    return out_rows, fail_rows, total, kept

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--pairs_csv", required=True)
    ap.add_argument("--out_dir", default="out_step3_5")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunksize", type=int, default=10)

    ap.add_argument("--probe_radius", type=float, default=1.4)  # water probe
    ap.add_argument("--n_points", type=int, default=100)        # accuracy/speed
    ap.add_argument("--sasa_thr", type=float, default=1.0)      # Å^2; use 0.1–2.0 typical
    ap.add_argument("--require_both", action="store_true")      # require both SG exposed
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # group by filename (treat each model separately)
    by_file = defaultdict(list)
    with open(args.pairs_csv, newline="") as f:
        r = csv.DictReader(f)
        in_fields = r.fieldnames
        for row in r:
            by_file[row["filename"]].append(row)

    fname2path, dup = build_filename_to_path_map(args.pdb_dir)
    if dup:
        print(f"WARNING: detected {dup} duplicate filenames in pdb_dir; using first occurrence per name.")

    tasks = []
    missing = 0
    for filename, rows in by_file.items():
        pdb_path = fname2path.get(filename)
        if pdb_path is None:
            missing += 1
            continue
        tasks.append((filename, rows, pdb_path))

    out_csv = os.path.join(args.out_dir, "step3_5_surface_pairs.csv")
    fail_csv = os.path.join(args.out_dir, "step3_5_failures.csv")

    out_fields = list(in_fields) + ["res1_sg_sasa","res2_sg_sasa","pair_best_sg_sasa"]
    fail_fields = ["pdb_id","filename","reason","detail"]

    total_pairs = 0
    kept_pairs = 0
    models_with_kept = 0

    worker = functools.partial(
        process_file,
        probe_radius=args.probe_radius,
        n_points=args.n_points,
        sasa_thr=args.sasa_thr,
        require_both=args.require_both
    )

    with open(out_csv, "w", newline="") as fo, open(fail_csv, "w", newline="") as ff:
        wo = csv.DictWriter(fo, fieldnames=out_fields)
        wf = csv.DictWriter(ff, fieldnames=fail_fields)
        wo.writeheader(); wf.writeheader()

        with Pool(processes=args.workers) as pool:
            for out_rows, fail_rows, t, k in pool.imap_unordered(worker, tasks, chunksize=args.chunksize):
                total_pairs += t
                kept_pairs += k
                if k > 0:
                    models_with_kept += 1
                for r in out_rows:
                    wo.writerow(r)
                for r in fail_rows:
                    wf.writerow(r)

    print("Step 3.5 (SG SASA surface proxy) complete.")
    print(f"Input models (files): {len(by_file)}  | missing pdb files: {missing}")
    print(f"Processed models:     {len(tasks)}")
    print(f"Pairs evaluated:      {total_pairs}")
    print(f"Pairs kept:           {kept_pairs}")
    print(f"Models with ≥1 kept pair: {models_with_kept}")
    print(f"Saved:\n  {out_csv}\n  {fail_csv}")

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
import os, csv, argparse, subprocess, shutil
from collections import defaultdict
from pathlib import Path
import functools
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree


def which_or_die(exe: str):
    p = shutil.which(exe)
    if p is None:
        raise RuntimeError(f"Required executable not found in PATH: {exe}")
    return p


def build_filename_to_path_map(pdb_dir: str):
    """
    Map exact basename (e.g., AF-O00244-F1-model_v6.pdb) -> full path.
    This avoids UniProt collisions when multiple AF files exist per UniProt.
    """
    m = {}
    for p in Path(pdb_dir).rglob("*.pdb"):
        m.setdefault(p.name, str(p))
    return m


def parse_msms_vert(vert_path: str) -> np.ndarray:
    verts = []
    with open(vert_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            verts.append((x, y, z))
    if not verts:
        raise RuntimeError("No vertices parsed from .vert")
    return np.asarray(verts, dtype=float)


def parse_sg_map_fast(pdb_path: str):
    """
    Map (chain, resseq, icode) -> SG coord.
    For robustness: keep the first SG seen per residue; prefer altLoc=' ' if present.
    """
    sg_map = {}
    alt_rank = {}  # track preference; lower is better
    rank = {" ": 0, "A": 1}  # blank preferred, then A, then others

    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            altloc    = line[16:17]
            resname   = line[17:20].strip()
            chain     = line[21:22].strip() or "?"
            resseq_s  = line[22:26].strip()
            icode     = (line[26:27].strip() or "")
            if resname != "CYS" or atom_name != "SG":
                continue
            try:
                resseq = int(resseq_s)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue

            key = (chain, resseq, icode)
            r = rank.get(altloc, 2)  # accept others but lower priority
            if key not in sg_map or r < alt_rank.get(key, 99):
                sg_map[key] = np.array([x, y, z], dtype=float)
                alt_rank[key] = r

    return sg_map


def run_pdb_to_xyzr(pdb_to_xyzr_exe: str, pdb_path: str, xyzr_path: str):
    with open(xyzr_path, "w") as out:
        r = subprocess.run([pdb_to_xyzr_exe, pdb_path], stdout=out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"pdb_to_xyzr failed: {r.stderr[:300]}")


def run_msms(msms_exe: str, xyzr_path: str, out_prefix: str, density: float):
    r = subprocess.run(
        [msms_exe, "-if", xyzr_path, "-of", out_prefix, "-density", str(density)],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    if r.returncode != 0:
        raise RuntimeError(f"msms failed: {r.stderr[:300]}")


def compute_depths(treeV: cKDTree, p1: np.ndarray, p2: np.ndarray):
    m = 0.5 * (p1 + p2)
    d1 = float(treeV.query(p1, k=1)[0])
    d2 = float(treeV.query(p2, k=1)[0])
    dm = float(treeV.query(m,  k=1)[0])
    return d1, d2, dm


def safe_tag(filename: str) -> str:
    # for scratch dir names
    return filename.replace("/", "_")


def process_file(task, pdb_to_xyzr_exe, msms_exe, density, scratch_root):
    """
    task = (filename, rows, pdb_path)
    Returns: (filename, out_rows, fail_rows)
    """
    filename, rows, pdb_path = task
    out_rows, fail_rows = [], []

    tag = safe_tag(Path(filename).stem)
    sub = os.path.join(scratch_root, tag)
    os.makedirs(sub, exist_ok=True)

    xyzr_path = os.path.join(sub, f"{tag}.xyzr")
    out_prefix = os.path.join(sub, f"{tag}_msms")
    vert_path = out_prefix + ".vert"

    try:
        run_pdb_to_xyzr(pdb_to_xyzr_exe, pdb_path, xyzr_path)
        run_msms(msms_exe, xyzr_path, out_prefix, density)

        V = parse_msms_vert(vert_path)
        treeV = cKDTree(V)

        sg_map = parse_sg_map_fast(pdb_path)

        for row in rows:
            pdb_id = row["pdb_id"]
            chain1 = row["chain1"]
            chain2 = row["chain2"]
            res1 = int(row["res1"])
            res2 = int(row["res2"])
            icode1 = (row.get("icode1","") or "").strip()
            icode2 = (row.get("icode2","") or "").strip()

            p1 = sg_map.get((chain1, res1, icode1))
            p2 = sg_map.get((chain2, res2, icode2))
            if p1 is None or p2 is None:
                fail_rows.append({
                    "pdb_id": pdb_id,
                    "filename": filename,
                    "reason": "sg_not_found",
                    "detail": f"Missing SG for ({chain1},{res1},{icode1}) or ({chain2},{res2},{icode2})"
                })
                continue

            d1, d2, dm = compute_depths(treeV, p1, p2)

            out = dict(row)
            out["res1_depth"] = f"{d1:.6f}"
            out["res2_depth"] = f"{d2:.6f}"
            out["motif_center_depth"] = f"{dm:.6f}"
            out_rows.append(out)

    except Exception as e:
        # use pdb_id from first row for context
        fail_rows.append({
            "pdb_id": rows[0]["pdb_id"],
            "filename": filename,
            "reason": "msms_depth_failed",
            "detail": str(e)[:300]
        })

    return (filename, out_rows, fail_rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--pairs_csv", required=True)
    ap.add_argument("--out_dir", default="out_step4_rerun")
    ap.add_argument("--density", type=float, default=3.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--chunksize", type=int, default=10)
    args = ap.parse_args()

    pdb_to_xyzr_exe = which_or_die("pdb_to_xyzr")
    msms_exe        = which_or_die("msms")

    os.makedirs(args.out_dir, exist_ok=True)
    scratch = os.path.join(args.out_dir, "_msms_tmp")
    os.makedirs(scratch, exist_ok=True)

    depth_csv = os.path.join(args.out_dir, "step4_depths.csv")
    fail_csv  = os.path.join(args.out_dir, "step4_failures.csv")

    # Group by *filename* (critical)
    by_file = defaultdict(list)
    with open(args.pairs_csv, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            by_file[row["filename"]].append(row)

    fname2path = build_filename_to_path_map(args.pdb_dir)

    tasks = []
    missing = []
    for filename, rows in by_file.items():
        pdb_path = fname2path.get(filename)
        if pdb_path is None:
            missing.append((rows[0]["pdb_id"], filename))
        else:
            tasks.append((filename, rows, pdb_path))

    out_fields = [
        "pdb_id","bench_name","filename","length",
        "n_cys","n_pairs_total","n_pairs_kept",
        "chain1","res1","icode1","chain2","res2","icode2","sg_sg_dist",
        "res1_depth","res2_depth","motif_center_depth"
    ]
    fail_fields = ["pdb_id","filename","reason","detail"]

    with open(depth_csv, "w", newline="") as fo, open(fail_csv, "w", newline="") as ff:
        wo = csv.DictWriter(fo, fieldnames=out_fields)
        wf = csv.DictWriter(ff, fieldnames=fail_fields)
        wo.writeheader()
        wf.writeheader()

        for pdb_id, filename in missing:
            wf.writerow({"pdb_id": pdb_id, "filename": filename, "reason": "pdb_not_found", "detail": "Missing exact filename in pdb_dir"})

        worker_func = functools.partial(
            process_file,
            pdb_to_xyzr_exe=pdb_to_xyzr_exe,
            msms_exe=msms_exe,
            density=args.density,
            scratch_root=scratch
        )

        n_total = len(tasks)
        n_done = 0

        with Pool(processes=args.workers) as pool:
            for filename, out_rows, fail_rows in pool.imap_unordered(worker_func, tasks, chunksize=args.chunksize):
                n_done += 1
                for r in out_rows:
                    wo.writerow(r)
                for fr in fail_rows:
                    wf.writerow(fr)

                if n_done % 250 == 0:
                    print(f"Processed {n_done}/{n_total} files")

    print("\nStep 4 (rerun) complete.")
    print(f"Saved:\n  {depth_csv}\n  {fail_csv}")


if __name__ == "__main__":
    main()


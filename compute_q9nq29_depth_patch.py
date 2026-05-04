#!/usr/bin/env python3
import csv
import numpy as np
from Bio.PDB import PDBParser, Polypeptide
from scipy.spatial import cKDTree

PDB_PATH  = "/home/dami/Downloads/copper_screening/human_proteome/AF-Q9NQ29-F1-model_v6.pdb"
PAIRS_CSV = "Q9NQ29_pairs_step3.csv"
VERT_FILE = "Q9NQ29_msms.vert"
OUT_CSV   = "Q9NQ29_step4_patch.csv"

STEP4_FIELDS = [
    "pdb_id","bench_name","filename","length","n_cys","n_pairs_total","n_pairs_kept",
    "chain1","res1","icode1","chain2","res2","icode2","sg_sg_dist",
    "res1_depth","res2_depth","motif_center_depth"
]

def load_msms_vertices(vert_path: str) -> np.ndarray:
    verts = []
    with open(vert_path, "r") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            # MSMS .vert usually starts each vertex line with x y z ...
            try:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            except Exception:
                continue
            verts.append((x, y, z))
    if not verts:
        raise RuntimeError(f"No vertices parsed from {vert_path}")
    return np.asarray(verts, dtype=float)

def build_sg_map(pdb_path: str):
    """
    Returns dict: (chain, resnum, icode) -> SG coord (np.array shape (3,))
    """
    parser = PDBParser(QUIET=True)
    st = parser.get_structure("Q9NQ29", pdb_path)
    model = st[0]

    sg = {}
    for r in model.get_residues():
        if not Polypeptide.is_aa(r):
            continue
        if r.get_resname() != "CYS":
            continue
        if "SG" not in r:
            continue
        chain = r.get_parent().id
        resnum = int(r.id[1])
        icode = (r.id[2] or "").strip()
        sg[(chain, resnum, icode)] = r["SG"].get_coord().astype(float)
    return sg

def main():
    surface = load_msms_vertices(VERT_FILE)
    tree = cKDTree(surface)
    sg_map = build_sg_map(PDB_PATH)

    n_in, n_out, n_missing = 0, 0, 0

    with open(PAIRS_CSV, newline="") as fin, open(OUT_CSV, "w", newline="") as fout:
        rdr = csv.DictReader(fin)
        w = csv.DictWriter(fout, fieldnames=STEP4_FIELDS)
        w.writeheader()

        for row in rdr:
            n_in += 1
            c1 = row["chain1"]; c2 = row["chain2"]
            r1 = int(row["res1"]); r2 = int(row["res2"])
            i1 = (row.get("icode1","") or "").strip()
            i2 = (row.get("icode2","") or "").strip()

            k1 = (c1, r1, i1)
            k2 = (c2, r2, i2)
            if k1 not in sg_map or k2 not in sg_map:
                n_missing += 1
                continue

            p1 = sg_map[k1]
            p2 = sg_map[k2]
            mid = 0.5 * (p1 + p2)

            d1 = float(tree.query(p1, k=1)[0])
            d2 = float(tree.query(p2, k=1)[0])
            dm = float(tree.query(mid, k=1)[0])

            out = {k: row.get(k, "") for k in STEP4_FIELDS}
            out["res1"] = str(r1)
            out["res2"] = str(r2)
            out["icode1"] = i1
            out["icode2"] = i2
            out["res1_depth"] = f"{d1:.6f}"
            out["res2_depth"] = f"{d2:.6f}"
            out["motif_center_depth"] = f"{dm:.6f}"

            w.writerow(out)
            n_out += 1

    print(f"Wrote: {OUT_CSV}")
    print(f"Pairs in: {n_in} | depth rows out: {n_out} | missing SG pairs skipped: {n_missing}")

if __name__ == "__main__":
    main()


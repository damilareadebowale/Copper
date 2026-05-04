# Human Copper Chaperone Screen from AlphaFold Structures

A structure-based workflow for identifying candidate human copper chaperones from the AlphaFold human proteome using cysteine geometry and surface accessibility.

## Rationale

Copper chaperones often contain cysteine pairs that are:

- present in the same structure
- close enough in space to support metal coordination
- sufficiently surface-accessible to interact with solvent or metal ions

This workflow screens the AlphaFold human proteome for such motifs in four steps:

1. input accounting
2. cysteine extraction
3. SG–SG distance filtering
4. MSMS-based depth calculation

The benchmark set is carried through the workflow to ensure known copper-binding motifs survive the filters.

---

## Workflow

### Step 1. Input and accounting

**Action**  
Iterate through every `.pdb` file in the AlphaFold human proteome library.

**Metric**  
Count total proteins/models started.

**Validation**  
Confirm that benchmark proteins exist in the library:

- ATOX1 (`O00244`)
- CCS (`O14618`)
- SCO1 (`O43819`)
- SCO2 (`O75880`)
- COX17 (`Q14061`)
- COX19 (`Q49B96`)

---

### Step 2. Parse each PDB and extract cysteines with SG atoms

**Action**  
For each structure, identify all residues annotated as `CYS` that contain an `SG` atom.

**Filter**  
Discard proteins/models with fewer than 2 cysteines.

**Metric**  
Count proteins/models with `>= 2` Cys(SG).

**Validation**  
Check that the benchmark proteins pass this step.

---

### Step 3. Pairwise SG–SG distances (2.0-5.5 Å)

**Action**  
For every protein/model that passes Step 2, compute all pairwise distances between cysteine sulfur atoms.

**Filter**  
Keep only pairs with:

- `SG–SG distance >= 2.0 Å`
- `SG–SG distance <= 5.5 Å`

**Metric**  
Count proteins/models with at least one valid SG–SG pair.

**Output**  
Emit one row per retained pair to a CSV file for inspection and downstream analysis.

**Validation**  
Verify that known benchmark motifs survive this step, for example:

- ATOX1: Cys12–Cys15
- CCS: Cys22–Cys25, Cys141–Cys227, Cys144–Cys227, Cys244–Cys246
- SCO1: Cys133–Cys137
- SCO2: Cys169–Cys173
- COX17: Cys23–Cys24, Cys26–Cys55, Cys36–Cys45

---

### Step 4. Depth calculation via MSMS

**Action**  
For proteins/models that survive Step 3:

1. generate `.xyzr`
2. run `msms` to build the solvent-excluded surface
3. build a KDTree from surface vertices
4. compute per-pair depth features

**Computed features**

- `res1_depth`: distance from cysteine 1 SG to nearest surface point
- `res2_depth`: distance from cysteine 2 SG to nearest surface point
- `motif_center_depth`: distance from the SG-pair midpoint to nearest surface point

**Result**  
A CSV containing every geometry-valid cysteine pair annotated with burial depth.

---

## Execution order

Run the workflow in two phases:

### Phase 1
Run **Step 1-3** first to:

- count the input library
- identify proteins with cysteine pairs
- validate benchmark motif recovery

### Phase 2
Run **Step 4** only on Step 3 survivors to:

- compute surface depth
- prioritize surface-accessible cysteine pairs

---

## Benchmark set

```python
BENCH = {
    "O00244": {"name": "ATOX1", "pairs": [(12, 15)]},
    "O14618": {"name": "CCS",   "pairs": [(22, 25), (141, 227), (144, 227), (244, 246)]},
    "O43819": {"name": "SCO1",  "pairs": [(133, 137)]},
    "O75880": {"name": "SCO2",  "pairs": [(169, 173)]},
    "Q14061": {"name": "COX17", "pairs": [(23, 24), (26, 55), (36, 45)]},
    "Q49B96": {"name": "COX19", "pairs": []},
}

Installation
Python dependencies
pip install numpy pandas scipy biopython

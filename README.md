# Structure-Based Screening of the AlphaFold Human Proteome for Candidate Copper-Binding Proteins

A Python-based workflow for proteome-scale identification and prioritization of human proteins containing structurally plausible cysteine environments for copper coordination.

The workflow combines cysteine geometry, structural accessibility, molecular-surface analysis, and benchmark tracking to generate candidate proteins for downstream experimental investigation.

## Scope

This workflow was developed to screen AlphaFold-predicted human protein structures for cysteine environments that may be compatible with copper coordination.

The workflow evaluates:

1. the presence of cysteine residues containing sulfur (`SG`) atoms;
2. pairwise sulfur-sulfur distances;
3. surface accessibility of candidate cysteine residues;
4. residue depth relative to the molecular surface;
5. motif-center depth; and
6. recovery of established human copper-handling proteins as benchmarks.

Structural matches identified by this workflow are **testable hypotheses for experimental investigation**.

They do not, by themselves, demonstrate:

- copper binding;
- copper specificity;
- copper-transfer activity; or
- copper-chaperone function.

Direct biochemical and biophysical validation is required to establish these properties.

## Adaptability

Although this implementation was developed for candidate copper-binding proteins, the workflow is modular.

The same general framework can be adapted to other proteome-scale structure-based screening problems by changing:

- residue identities;
- geometric criteria;
- distance thresholds;
- accessibility criteria;
- benchmark proteins; and
- downstream validation rules.

This makes the workflow potentially useful for investigating other metal-binding environments, catalytic residue arrangements, accessible residue clusters, and related structural motifs.

## Workflow

### Step 1. Input accounting

The workflow recursively identifies `.pdb` files in an AlphaFold human proteome directory.

For each structural model, it records the associated UniProt accession and tracks the total number of models processed.

Known human copper-handling proteins are tracked throughout the workflow as benchmarks.

### Step 2. Cysteine extraction

Each structure is parsed and all cysteine residues containing an `SG` atom are identified.

Proteins containing fewer than two cysteine `SG` atoms are excluded from pairwise analysis.

### Step 3. Pairwise sulfur-sulfur distance analysis

All pairwise distances between cysteine sulfur atoms are calculated.

By default, candidate pairs are retained when:

```text
2.0 Å <= SG-SG distance <= 5.5 Å
```

The minimum and maximum distance thresholds can be changed using command-line arguments.

The output contains one row for each retained cysteine pair.

### Step 3.5. Optional SG surface-accessibility filtering

Candidate pairs can be further evaluated using the solvent-accessible surface area of cysteine sulfur atoms.

The script uses Biopython's Shrake-Rupley implementation to estimate SG atom solvent-accessible surface area.

This step can be used as an additional surface-exposure filter before MSMS depth analysis.

By default:

```text
probe radius = 1.4 Å
surface points = 100
SG SASA threshold = 1.0 Å²
```

The thresholds can be modified through command-line options.

### Step 4. Molecular-surface depth calculation

Candidate structures are processed using MSMS.

For each structural model:

1. the PDB structure is converted to `.xyzr` format;
2. MSMS generates the solvent-excluded molecular surface;
3. surface vertices are loaded into a KD-tree;
4. each candidate cysteine sulfur atom is located relative to the surface; and
5. residue and motif-center depth values are calculated.

The following features are reported:

- `res1_depth`: distance from cysteine 1 SG atom to the nearest molecular-surface point;
- `res2_depth`: distance from cysteine 2 SG atom to the nearest molecular-surface point;
- `motif_center_depth`: distance from the midpoint of the SG-SG pair to the nearest molecular-surface point.

These measurements provide a structural estimate of how exposed or buried each candidate motif is.

## Benchmark proteins

The following established human copper-handling proteins are tracked during the analysis:

| Protein | UniProt accession | Benchmark cysteine pair(s) |
|---|---|---|
| ATOX1 | O00244 | Cys12-Cys15 |
| CCS | O14618 | Cys22-Cys25, Cys141-Cys227, Cys144-Cys227, Cys244-Cys246 |
| SCO1 | O43819 | Cys133-Cys137 |
| SCO2 | O75880 | Cys169-Cys173 |
| COX17 | Q14061 | Cys23-Cys24, Cys26-Cys55, Cys36-Cys45 |
| COX19 | Q49B96 | tracked as benchmark |

The benchmark set is used to determine whether established copper-associated cysteine environments survive the structural filters.

## Requirements

### Python

Python 3 is recommended.

Install the required Python packages using:

```bash
pip install -r requirements.txt
```

The current Python dependencies are:

```text
biopython
numpy
pandas
scipy
```

### External software

The MSMS depth-analysis step requires:

- `msms`
- `pdb_to_xyzr`

Both executables must be available on your system `PATH`.

You can check this using:

```bash
which msms
which pdb_to_xyzr
```

## Input data

The primary input is a directory containing AlphaFold-predicted human protein structures in PDB format.

The scripts search recursively, so the PDB files may be stored in a flat directory or within subdirectories.

AlphaFold-style filenames are expected, for example:

```text
AF-O00244-F1-model_v6.pdb
```

The UniProt accession is extracted from the filename.

## Quick Start

### Step 1 to Step 3: cysteine extraction and geometry screening

```bash
python3 screen_step1_3.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --out_dir out_step1_3 \
  --workers 4
```

Default distance thresholds are:

```text
minimum SG-SG distance: 2.0 Å
maximum SG-SG distance: 5.5 Å
```

They can be changed, for example:

```bash
python3 screen_step1_3.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --out_dir out_step1_3 \
  --min_d 2.0 \
  --max_d 5.5 \
  --workers 4
```

### Optional Step 3.5: SG surface-accessibility filtering

```bash
python3 step3_5_filter_surface_sasa_mp.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --pairs_csv out_step1_3/step3_pairs_check.csv \
  --out_dir out_step3_5 \
  --workers 4
```

To require both cysteine SG atoms to meet the SASA threshold:

```bash
python3 step3_5_filter_surface_sasa_mp.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --pairs_csv out_step1_3/step3_pairs_check.csv \
  --out_dir out_step3_5 \
  --sasa_thr 1.0 \
  --require_both \
  --workers 4
```

### Step 4: MSMS depth calculation

If Step 3.5 is not used:

```bash
python3 step4_depth_msms_mp.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --pairs_csv out_step1_3/step3_pairs_check.csv \
  --out_dir out_step4 \
  --workers 4
```

If Step 3.5 is used, the filtered pair table can instead be supplied:

```bash
python3 step4_depth_msms_mp.py \
  --pdb_dir /path/to/alphafold_human_proteome \
  --pairs_csv out_step3_5/step3_5_surface_pairs.csv \
  --out_dir out_step4 \
  --workers 4
```

## Expected outputs

### Step 1 to Step 3

```text
out_step1_3/
├── step2_counts.csv
├── step3_counts.csv
└── step3_pairs_check.csv
```

`step2_counts.csv`

Contains cysteine counts for parsed protein structures.

`step3_counts.csv`

Contains the total number of possible cysteine pairs and the number retained by the SG-SG distance filter.

`step3_pairs_check.csv`

Contains one row for each geometry-valid cysteine pair.

### Optional Step 3.5

```text
out_step3_5/
├── step3_5_surface_pairs.csv
└── step3_5_failures.csv
```

`step3_5_surface_pairs.csv`

Contains cysteine pairs passing the selected SG surface-accessibility criterion.

### Step 4

```text
out_step4/
├── step4_depths.csv
└── step4_failures.csv
```

`step4_depths.csv`

Contains the geometry-valid candidate pairs together with:

```text
res1_depth
res2_depth
motif_center_depth
```

`step4_failures.csv`

Records structural models that could not be processed successfully during MSMS analysis.

## Repository files

```text
.
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── screen_step1_3.py
├── step3_5_filter_surface_sasa_mp.py
├── step4_depth_msms_mp.py
├── compute_q9nq29_depth_patch.py
├── verify_counts.py
├── failed_ids.txt
└── test.xyzr
```

### Main workflow scripts

`screen_step1_3.py`

Performs input accounting, cysteine extraction, SG-SG distance screening, and benchmark tracking.

`step3_5_filter_surface_sasa_mp.py`

Provides an optional SG solvent-accessibility filter.

`step4_depth_msms_mp.py`

Calculates cysteine residue depth and motif-center depth using MSMS-generated solvent-excluded surfaces.

### Utility and supporting files

`compute_q9nq29_depth_patch.py`

A project-specific repair utility used to recompute depth information for an individual structural model. It is not required for the standard workflow.

`verify_counts.py`

A dataset-accounting utility used to inspect retained pairs, unique proteins, and output structure. The input path in this script should be updated to match the final Step 4 table being evaluated.

`failed_ids.txt`

Records identifiers for structural models that could not be processed successfully during parts of the workflow.

`test.xyzr`

A small XYZR-format test file that can be used to inspect or test MSMS-compatible molecular-surface input.

## Important notes

- Each AlphaFold PDB filename is treated as a separate structural model during calculation.
- Multiple AlphaFold fragments may map to the same UniProt accession.
- Results can therefore be summarized downstream by structural model or collapsed by UniProt accession.
- MSMS failures should be reviewed separately and may be rerun when appropriate.
- Candidate ranking criteria should be chosen according to the biological question being investigated.
- Structural geometry alone is not evidence of metal binding or biological function.

## Interpretation

The purpose of this workflow is candidate prioritization.

A protein that passes the structural filters should be interpreted as containing a cysteine environment that is geometrically and structurally compatible with the criteria used in the screen.

It should not be interpreted as proof that the protein:

- binds copper;
- binds copper selectively;
- transfers copper;
- functions as a copper chaperone; or
- has a copper-dependent biological role.

Candidate proteins require direct experimental validation.

## Reuse for other structural-screening projects

The workflow can be adapted for other questions by modifying the screening rules.

Examples include:

- alternative metal-coordination environments;
- different residue combinations;
- catalytic residue arrangements;
- accessible residue clusters;
- redox-active cysteine environments; and
- other proteome-scale structural motifs.

Users should define new residue, geometry, accessibility, benchmark, and validation criteria appropriate for the new biological question.

## Citation

If you use this workflow, please cite the software record together with the relevant AlphaFold and MSMS resources.

The repository includes a `CITATION.cff` file containing citation metadata.

Software citation:

```text
Adebowale, D. D. (2026).
Structure-Based Screening of the AlphaFold Human Proteome for Candidate Copper-Binding Proteins.
Version 1.0.0.
Zenodo.
https://doi.org/10.5281/zenodo.22349882
```

APA citation:

```text
Adebowale, D. D. (2026). Structure-Based Screening of the AlphaFold Human Proteome for Candidate Copper-Binding Proteins (Version v1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22349882
```

## License

This project is distributed under the MIT License. See the `LICENSE` file for details.

## Author

**Damilare Desmond Adebowale**  
Purdue University

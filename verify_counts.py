import pandas as pd

# Load your final Step 4 output
INPUT_CSV = "out_step4/step4_depths_FIXED_dedup.csv"

def verify():
    print(f"Reading {INPUT_CSV}...")
    df = pd.read_csv(INPUT_CSV)
    
    # 1. Total Pairs (Rows)
    total_pairs = len(df)
    
    # 2. Total Unique Proteins (Unique PDB IDs)
    unique_proteins = df['pdb_id'].nunique()
    
    # 3. Multi-pair analysis
    pair_counts = df['pdb_id'].value_counts()
    multi_pair_proteins = pair_counts[pair_counts > 1]
    
    print("\n" + "="*40)
    print("      DATASET ACCOUNTING")
    print("="*40)
    print(f"Total Retained Pairs (Rows):    {total_pairs:,}")
    print(f"Total Unique Proteins:          {unique_proteins:,}")
    print(f"Proteins with Multiple Pairs:   {len(multi_pair_proteins):,}")
    print("="*40)
    
    # Show examples of proteins with the most Cys-pairs
    print("\nTop 5 Proteins by Cys-Pair Count:")
    print(pair_counts.head(5).to_string())

    # Verify output format (Show first 2 rows)
    print("\nData Format Check (First 2 Rows):")
    print(df[['pdb_id', 'res1', 'res2', 'sg_sg_dist', 'motif_center_depth']].head(2).to_string(index=False))

if __name__ == "__main__":
    verify()

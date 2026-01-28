# eDNA Metabarcoding Pipeline for MinION

This pipeline is designed for dual-marker amplicon sequencing (specifically 18S rRNA and COI) from aquatic or terrestrial samples. It takes raw FASTQ files and produces clean and abundance-based OTU matrices ready for statistical analysis.

## Key Features

* **Marker-Aware:** Automatically separates 18S and COI reads based on amplicon length.
* **Noise Reduction:** Uses `filtlong` for quality filtering and VSEARCH for chimera removal.
* **High Accuracy:** Generates Consensus Sequences for OTUs to correct Nanopore sequencing errors.
* **Biologically Valid:** Filters out PCR artifacts (chimeras) and genomic concatemers.
* **Taxonomy Ready:** Outputs clean consensus sequences compatible with SINTAX or BLAST.

## Pipeline Overview

![eDNA Metabarcoding Pipeline Workflow](figures/pipeline_for_edna.png)

### Pipeline Steps

1. **Quality Filtering** (`1_run_preprocessing.py`) - Filter reads using `filtlong` (Q>20).
2. **Marker Classification** (`2_classify_markers.py`) - Separate 18S (1500-2800bp) and COI (300-1000bp) reads.
3. **Clustering + Chimera Detection** (`3_run_clustering_by_marker.py`) - VSEARCH clustering (95% identity) and consensus generation.
4. **Abundance Matrix Generation** (`4_merge_otu_tables_by_marker.py`) - Create OTU count tables.
5. **Taxonomy Assignment** (`5_assign_taxonomy.py`) - SINTAX classification against SILVA (18S) or MIDORI (COI).
6. **Validation** (`6_blast_top_otus.py`) - Manual BLAST checks for top OTUs (useful for forensics).
7. **Reporting** (`7_comprehensive_taxonomy_summary.py`) - Merges abundance, taxonomy, and BLAST results into a final Master CSV.

## Installation and Setup

### 1. Prerequisites

You need a Unix-based system (macOS or Linux). The pipeline runs in a local Conda environment to ensure reproducibility.

### 2. Set up the Environment

You can create the required environment using the provided script or manually.

**Option A: Automated Setup**

```bash
bash create_env.sh
```

**Option B: Manual Setup**

```bash
conda env create -f environment.yml -p ./env
```

Tools installed: filtlong, vsearch (v2.30+), python (3.10), and standard system tools.

## Database Setup

Before running taxonomy assignment (Step 5), you need to download and prepare reference databases.

### 18S rRNA Database (SILVA)

1. **Download the SILVA 18S database:**
   
   Source: [USEARCH SINTAX Downloads](https://www.drive5.com/usearch/manual/sintax_downloads.html)
   
   ```bash
   cd refs/
   wget https://www.drive5.com/sintax/silva_18s_v123.fa.gz
   gunzip silva_18s_v123.fa.gz
   ```

2. **Convert to VSEARCH UDB format:**
   
   ```bash
   ./env/bin/vsearch --makeudb_usearch silva_18s_v123.fa --output silva_18s_v123.udb
   ```
   
   This creates a binary database (~1.1 GB) optimized for VSEARCH SINTAX.

### COI Database (To be determined)

**Status:** No suitable COI reference database has been identified yet.

Once a database is selected, convert it to UDB format using:
```bash
./env/bin/vsearch --makeudb_usearch database.fasta --output database_COI.udb
```

## How to Run

### Quick Start (Background Mode)

The recommended way to run the full pipeline is in the background. This ensures the process does not stop if you close your terminal.

```bash
# 1. Activate the environment
conda activate ./env

# 2. Run the orchestration script
# Replace --root with the path to your raw 'fastq_pass' folder
nohup bash scripts/run_full_pipeline.sh \
    --root data/Water_eDNA_18S_COI_14_01_26/fastq_pass \
    --threads 12 \
    > pipeline_run.log 2>&1 &

# 3. Monitor progress
tail -f pipeline_run.log
```

**Example Output:** See [logs_example/pipeline_20260128_124748.log](logs_example/pipeline_20260128_124748.log) for a complete pipeline execution log showing all steps, timing, and resource usage.

#### Running Step-by-Step (For Debugging)

If you need to test specific steps or run only one part of the analysis:

```bash
# Step 1: Quality Filter (Per Sample)
# Note: Use --input_files with a wildcard to catch all chunks
python3 scripts/1_run_preprocessing.py \
    --input_files data/sample_01/*.fastq.gz \
    --output_dir out/sample_01

# Step 2: Separate 18S and COI
python3 scripts/2_classify_markers.py --input_dir "out/run_name"

# Step 3: Cluster and Remove Chimeras
python3 scripts/3_run_clustering_by_marker.py \
    --input_dir "out/run_name" \
    --output_dir "out/run_name" \
    --threads 12

# Step 4: Create Abundance Matrices
python3 scripts/4_merge_otu_tables_by_marker.py --input_dir "out/run_name"

# Step 5: Assign Taxonomy (Example for 18S)
python3 scripts/5_assign_taxonomy.py \
    --input_dir "out/run_name" \
    --db_18S refs/silva_18s_v123.udb \
    --threads 12

# Step 7: Generate Final Report
python3 scripts/7_comprehensive_taxonomy_summary.py --input_dir "out/run_name"
```

## Output Files and Interpretation

After the pipeline finishes, your results will be in `out/<run_name>/`.

### 1. The Final Matrices (`out/<run_name>/merged/`)

These are your primary results for statistical analysis (R/Python).

- **`otu_abundance_matrix_18S.csv`:** Raw read counts for 18S. Rows = OTUs, Columns = Samples.
- **`otu_relative_abundance_18S.csv`:** Normalized counts (0 to 1). Use this to compare community composition between samples with different sequencing depths.

Note: The same files exist for the COI marker.

### 2. The Biological Sequences (`out/<run_name>/temp_clustering/`)

Use these files to identify the taxonomy of your OTUs.

- **`consensus_18S_clean.fasta`:** The high-quality consensus sequences for every OTU. Chimeras have already been removed from this file.

Header Format: `>OTU_18S_00001;size=500`. The size indicates how many reads belong to this OTU.

### 3. Quality Control Logs (`out/<run_name>/logs/`)

- **`global_clustering.log`:** Check this to see how many chimeras were removed and how many valid OTUs were found.
- **`barcodeXX.log`:** Check these if a specific sample has very low read counts.

### 4. Taxonomy Summary (`out/<run_name>/taxonomy_summary/`)

These are the master files combining abundance, taxonomy assignments and BLAST validation results.

- **`comprehensive_taxonomy_18S.csv`:** Complete taxonomy table for 18S marker. Each row is an OTU with columns for:
  - OTU ID and total abundance
  - Per-sample abundance (all barcode columns)
  - SILVA taxonomy at all levels (Domain, Phylum, Class, Order, Family, Genus, Species)
  - BLAST results (if validation was run)
- **`comprehensive_taxonomy_COI.csv`:** Same format for COI marker.

These files are ready for direct import into R or Python for statistical analysis and visualization.

## Methodology Summary

- **Filtering:** Uses filtlong to remove reads with mean Quality Score < 20. Keeps top 100% of data (no subsampling).
- **Classification:** Uses Python to separate markers by length (18S: 1500-2800bp, COI: 300-1000bp).
- **Clustering:** Uses vsearch to cluster reads at 95% identity. Generates Consensus Sequences to fix Nanopore errors.
- **De-Noising:** Uses uchime_denovo to detect and remove chimeric sequences (PCR artifacts) using the consensus sequences.
- **Quantification:** Maps reads back to valid OTUs to create abundance matrices.

## Visualization

A Jupyter Notebook [`Results_Analysis.ipynb`](Results_Analysis.ipynb) is provided results exploration and visualization, including:
* Community Composition (Stacked Bar Charts)
* Top Genera/Species (Bar Charts)
* Sample Similarity (Heatmaps)
* Forensics of Failed Markers (BLAST Analysis)

## License and Credits

This workflow was developed for the [Genorobotics](https://make.epfl.ch/projects/14/make-genorobotics-14) semester project (EPFL).

### References

- **ONT-AmpSeq:** The clustering and consensus generation strategy is heavily inspired by the [ONT-AmpSeq pipeline](https://github.com/michoug/ONT-AmpSeq).
- **SILVA Database:** 18S rRNA reference database from [SILVA v123](https://www.arb-silva.de/), pre-formatted for SINTAX available at [USEARCH SINTAX Downloads](https://www.drive5.com/usearch/manual/sintax_downloads.html)

# eDNA Metabarcoding Pipeline for MinION

This pipeline is designed for multi-marker amplicon sequencing (18S rRNA, COI, and JEDI) from aquatic or terrestrial eDNA samples. It takes raw FASTQ files and produces clean and abundance-based OTU matrices ready for statistical analysis.

## Supported Markers

| Marker | Target Gene | Amplicon Size | Typical Use | Reference Database |
|--------|------------|---------------|-------------|-------------------|
| **18S** | 18S rRNA | 1,500–2,800 bp | Water eDNA (eukaryotes) | SILVA v123 |
| **COI** | Cytochrome Oxidase I | 500–900 bp | Water/Soil (invertebrates, Folmer primers) | MIDORI2 / eKOI |
| **JEDI** | COI (short) | 250–500 bp | Soil eDNA (arthropods, JEDI primers ~460 bp) | MIDORI2 / eKOI (same COI db) |

## Next Student: What to Test First

Please review the technical limitations and to-do list in [results_analysis.ipynb](results_analysis.ipynb) and focus on these immediate tests:

1. **COI reference database:** Download and validate the **eKOI COI database** (see the COI Database section below) and confirm that SINTAX/BLAST assignments improve.
2. **Test Soil JEDI data:** Run the pipeline on the Soil dataset using `--markers JEDI,COI`. JEDI uses the same MIDORI2/eKOI database since it also targets COI.
3. **COI yield troubleshooting:** The previous >90% COI reads <300 bp issue may have been JEDI amplicons being misclassified. Verify with the new JEDI size range (250-500 bp).
4. **Blocking primers:** Evaluate blocking primers to reduce ciliate amplification in COI.
5. **Local reference library:** Begin compiling local Swiss taxa for a custom COI/18S reference to reduce mis-assignments.
6. **Dereplication:** Implement a dereplication step before clustering to speed up processing and reduce memory usage.
7. **Primer-based classification:** For datasets where JEDI and COI size ranges overlap, implement primer sequence detection as an alternative to length-based classification.

## Key Features

* **Marker-Aware:** Automatically separates 18S, COI, and JEDI reads based on amplicon length. Use `--markers` to select which markers to search for.
* **Noise Reduction:** Uses `filtlong` for quality filtering and VSEARCH for chimera removal.
* **High Accuracy:** Generates Consensus Sequences for OTUs to correct Nanopore sequencing errors.
* **Biologically Valid:** Filters out PCR artifacts (chimeras) and genomic concatemers.
* **Taxonomy Ready:** Outputs clean consensus sequences compatible with SINTAX or BLAST.

## Pipeline Overview

![eDNA Metabarcoding Pipeline Workflow](figures/pipeline_for_edna.png)

### Pipeline Steps

1. **Quality Filtering** (`1_run_preprocessing.py`) - Filter reads using `filtlong` (Q>20).
2. **Marker Classification** (`2_classify_markers.py`) - Separate reads by marker using `--markers` flag:
   - 18S: 1500–2800 bp | COI: 500–900 bp | JEDI: 250–500 bp
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

### COI Database (eKOI / MIDORI2)

The **eKOI database** is a curated COI reference option: https://academic.oup.com/database/article/doi/10.1093/database/baaf057/8263868

Once downloaded, convert it to UDB format using:
```bash
./env/bin/vsearch --makeudb_usearch eKOI.fasta --output eKOI_COI.udb
```

### JEDI Database

**JEDI primers target COI** (shorter amplicon ~460 bp), so the same COI reference database works. No separate database is needed — point `--db_JEDI` to your MIDORI2 or eKOI COI database:
```bash
# Example: use MIDORI2 for both COI and JEDI
python3 scripts/5_assign_taxonomy.py \
    --input_dir out/run_name \
    --db_COI refs/midori2_COI.udb \
    --db_JEDI refs/midori2_COI.udb
```

## How to Run

### Quick Start (Background Mode)

The recommended way to run the full pipeline is in the background. This ensures the process does not stop if you close your terminal.

```bash
# 1. Activate the environment
conda activate ./env

# 2. Run the orchestration script
# Water eDNA (18S + COI markers)
nohup bash scripts/run_full_pipeline.sh \
    --root data/Water_eDNA_18S_COI_14_01_26/fastq_pass \
    --markers 18S,COI \
    --threads 14 \
    > pipeline_water.log 2>&1 &

# Soil eDNA (JEDI + COI markers)
nohup bash scripts/run_full_pipeline.sh \
    --root data/Soil_eDNA_JEDI_COI_14_01_26/fastq_pass \
    --markers JEDI,COI \
    --threads 14 \
    > pipeline_soil.log 2>&1 &

# 3. Monitor progress
tail -f pipeline_water.log
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

# Step 2: Separate markers (choose your marker set)
python3 scripts/2_classify_markers.py --input_dir "out/run_name" --markers 18S,COI
# For soil data:
python3 scripts/2_classify_markers.py --input_dir "out/run_name" --markers JEDI,COI

# Step 3: Cluster and Remove Chimeras
python3 scripts/3_run_clustering_by_marker.py \
    --input_dir "out/run_name" \
    --output_dir "out/run_name" \
    --markers JEDI,COI \
    --threads 12

# Step 4: Create Abundance Matrices
python3 scripts/4_merge_otu_tables_by_marker.py --input_dir "out/run_name"

# Step 5: Assign Taxonomy (databases matched to markers)
# For water (18S + COI):
python3 scripts/5_assign_taxonomy.py \
    --input_dir "out/run_name" \
    --db_18S refs/silva_18s_v123.udb \
    --db_COI refs/midori2_COI.udb \
    --threads 12
# For soil (JEDI + COI) — JEDI uses the same COI database:
python3 scripts/5_assign_taxonomy.py \
    --input_dir "out/run_name" \
    --db_JEDI refs/midori2_COI.udb \
    --db_COI refs/midori2_COI.udb \
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
- **Classification:** Uses Python to separate markers by length (18S: 1500-2800bp, COI: 500-900bp, JEDI: 250-500bp). Marker set is configurable via `--markers`.
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

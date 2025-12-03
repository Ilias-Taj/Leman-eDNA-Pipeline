#!/usr/bin/env bash
set -eo pipefail
# run_full_pipeline.sh
# Run preprocessing -> alignment -> matrix -> validation for all samples in the raw_reads folder.
# Writes per-sample logs to out/logs/<sample>.log and produces outputs under out/<sample>/ and out/validation/<sample>/

ROOT_DIR="raw_reads"
ENV_PREFIX="./env"
THREADS=4
MIN_READS=20
MAPQ=40
KEEP_PERCENT=90
MIN_LENGTH=100
DORADO_MODEL=""
DORADO_DEVICE=""
DORADO_DRY_RUN=0

usage(){
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --root DIR         Root folder with sample subfolders (default: "$ROOT_DIR")
  --env PREFIX       Conda env prefix to use (default: "$ENV_PREFIX")
  --threads N        Threads for minimap2/samtools (default: $THREADS)
  --min_reads N      Minimum reads to attempt consensus (default: $MIN_READS)
  --mapq N           MAPQ threshold for grouping reads (default: $MAPQ)
  --keep_percent N   Filtlong keep percent (default: $KEEP_PERCENT)
  --min_length N     Minimum read length after trimming (default: $MIN_LENGTH)
  -h, --help         Show this help

Example:
  $(basename "$0") --env ./env --threads 4
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT_DIR="$2"; shift 2;;
    --env) ENV_PREFIX="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --min_reads) MIN_READS="$2"; shift 2;;
    --mapq) MAPQ="$2"; shift 2;;
    --keep_percent) KEEP_PERCENT="$2"; shift 2;;
    --min_length) MIN_LENGTH="$2"; shift 2;;
      --dorado_model) DORADO_MODEL="$2"; shift 2;;
      --dorado_device) DORADO_DEVICE="$2"; shift 2;;
      --dorado_dry_run) DORADO_DRY_RUN=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

mkdir -p out/logs

echo "Root: $ROOT_DIR"
echo "Using conda env prefix: $ENV_PREFIX"
echo "Threads: $THREADS, min_reads: $MIN_READS, mapq: $MAPQ"

# Find all sample directories (immediate children)
samples=()
while IFS= read -r -d $'\0' d; do
  samples+=("$(basename "$d")")
done < <(find "$ROOT_DIR" -maxdepth 1 -mindepth 1 -type d -print0 | sort -z)

if [ ${#samples[@]} -eq 0 ]; then
  echo "No sample directories found under: $ROOT_DIR" >&2
  exit 1
fi

echo "Found samples: ${samples[*]}"

for sample in "${samples[@]}"; do
  echo
  echo "========== Processing: $sample =========="
  sample_dir="$ROOT_DIR/$sample"
  outdir="out/$sample"
  raw_reads_path="$sample_dir/assembly/raw.reads.unsorted"
  logf="out/logs/${sample}.log"

  # create per-sample folders early so basecaller can write into outdir
  mkdir -p "$outdir" "out/logs" "out/validation/$sample"
  echo "Sample: $sample" > "$logf"

  # If raw signal (.pod5) files exist, run Dorado basecalling (step 0) and use
  # the produced FASTQ as the reads input. Otherwise, start from the existing
  # raw reads file (assembly/raw.reads.unsorted).
  pod5_dir="$sample_dir/pod5"
  basecalled="$outdir/basecalled.fastq.gz"
  if [ -d "$pod5_dir" ] && [ -n "$(find "$pod5_dir" -type f -name '*.pod5' -print -quit)" ]; then
    echo "Found .pod5 files; running Dorado basecaller..." | tee -a "$logf"
    # Build dorado args
    dorado_args=(--pod5_dir "$pod5_dir" --output "$basecalled")
    if [ -n "$DORADO_MODEL" ]; then
      dorado_args+=(--model "$DORADO_MODEL")
    fi
    if [ -n "$DORADO_DEVICE" ]; then
      dorado_args+=(--device "$DORADO_DEVICE")
    fi
    if [ "$DORADO_DRY_RUN" -eq 1 ]; then
      dorado_args+=(--dry_run)
    fi
    if ! conda run -p "$ENV_PREFIX" python3 scripts/0_dorado_basecall.py "${dorado_args[@]}" >> "$logf" 2>&1; then
      echo "Basecalling failed for $sample (see $logf). Continuing to next sample." | tee -a "$logf"
      continue
    fi
    if [ -f "$basecalled" ]; then
      reads="$basecalled"
    else
      echo "Basecalled file not found for $sample; skipping sample" | tee -a "$logf"
      continue
    fi
  else
    reads="$raw_reads_path"
  fi
  outdir="out/$sample"
  logf="out/logs/${sample}.log"

  mkdir -p "$outdir" "out/logs" "out/validation/$sample"
  if [ ! -f "$reads" ]; then
    echo "Missing reads file: $reads" | tee -a "$logf"
    continue
  fi

  # Choose primer pair based on sample name / gene
  if [[ "$sample" == *_16S* ]]; then
    FWD='CGCCTGTTTATCAAAAACAT'
    REV='CCGGTTTGAACTCAGATCA'
  elif [[ "$sample" == Leptopelis_vermiculatus_CO1 ]]; then
    FWD='CAATACCAAACCCCCTTRTTYGTWTGATC'
    REV='GCTTCTCARATAATAAATATYAT'
  elif [[ "$sample" == Rhynchocyon_udzungwensis_CO1 ]]; then
    FWD='GGTCAACAAATCATAAAGATATTGG'
    REV='TAAACTTCAGGGTGACCAAAAAATCA'
  else
    # fallback: use 16S primers (safe default for vertebrate 16S)
    FWD='CGCCTGTTTATCAAAAACAT'
    REV='CCGGTTTGAACTCAGATCA'
  fi

  echo "Primers: $FWD / $REV" | tee -a "$logf"

  # 1) Preprocess (cutadapt -> filtlong)
  echo "Running preprocessing (trim + filter) ..." | tee -a "$logf"
  if ! conda run -p "$ENV_PREFIX" python3 scripts/2_run_preprocessing.py --input_file "$reads" --output_dir "$outdir" --primer_f "$FWD" --primer_r "$REV" --min_length "$MIN_LENGTH" --keep_percent "$KEEP_PERCENT" --threads "$THREADS" >> "$logf" 2>&1; then
    echo "Preprocessing failed for $sample (see $logf). Continuing to next sample." | tee -a "$logf"
    continue
  fi

  filtered="$outdir/filtered_reads.fastq.gz"
  if [ ! -f "$filtered" ]; then
    echo "Filtered file not found for $sample; skipping alignment" | tee -a "$logf"
    continue
  fi

  # 2) Align filtered reads to DB
  echo "Running alignment..." | tee -a "$logf"
  if ! conda run -p "$ENV_PREFIX" python3 scripts/3_run_alignment.py --reads_file "$filtered" --db_file refs/sanger_db.fa --output_dir "$outdir" --threads "$THREADS" >> "$logf" 2>&1; then
    echo "Alignment failed for $sample (see $logf). Continuing to next sample." | tee -a "$logf"
    continue
  fi

  paf="$outdir/alignments.paf"
  if [ ! -f "$paf" ]; then
    echo "PAF not produced for $sample; skipping matrix and validation" | tee -a "$logf"
    continue
  fi

  # 3) Generate per-sample species counts
  echo "Generating species counts (matrix) ..." | tee -a "$logf"
  if ! conda run -p "$ENV_PREFIX" python3 scripts/4_generate_matrix.py --paf "$paf" --output "$outdir/species_counts.csv" >> "$logf" 2>&1; then
    echo "generate_matrix failed for $sample (see $logf)" | tee -a "$logf"
  fi

  # 4) Run consensus-based validation
  echo "Running consensus validation ..." | tee -a "$logf"
  if ! conda run -p "$ENV_PREFIX" python3 scripts/5_run_validation.py --paf_file "$paf" --reads_file "$reads" --db_file refs/sanger_db.fa --output_dir "out/validation/$sample" --mapq "$MAPQ" --min_reads "$MIN_READS" --threads "$THREADS" >> "$logf" 2>&1; then
    echo "Validation failed for $sample (see $logf)" | tee -a "$logf"
  fi

  echo "Finished sample: $sample" | tee -a "$logf"
done

# After processing all samples, merge per-sample matrices into a single matrix
echo "Running merge and comparison step..." | tee -a out/logs/run_full_pipeline_merge.log
if ! conda run -p "$ENV_PREFIX" python3 scripts/4_merge_and_compare.py >> out/logs/run_full_pipeline_merge.log 2>&1; then
  echo "Merge step failed (see out/logs/run_full_pipeline_merge.log)" >&2
else
  echo "Merge completed. Outputs in out/merged/" | tee -a out/logs/run_full_pipeline_merge.log
fi

echo "All done. Per-sample logs are in out/logs/"

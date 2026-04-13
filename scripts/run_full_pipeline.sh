#!/usr/bin/env bash
set -eo pipefail
# run_full_pipeline.sh
# Full eDNA metabarcoding pipeline: quality filtering → marker classification → clustering + chimera detection → abundance matrices
# Processes already basecalled & demultiplexed FASTQ files.
# Writes per-sample logs to out/logs/<sample>.log and produces outputs under out/<sample>/

# Default to the delivered eDNA run; override with --root if needed
ROOT_DIR="data/Water_eDNA_18S_COI_14_01_26/fastq_pass" #MODIFY THIS TO CHECK OTHER DATASETS
ENV_PREFIX="./env"
THREADS=12
MIN_READS=20
MAPQ=20 # MAPQ threshold for grouping reads
KEEP_PERCENT=100
MIN_LENGTH=0
MIN_MEAN_Q=20
MARKERS="18S,COI" # Comma-separated markers: 18S,COI,JEDI (e.g. "JEDI,COI" for soil data)

# Detect system resources
if [[ "$OSTYPE" == "darwin"* ]]; then
  CPU_CORES=$(sysctl -n hw.ncpu)
  TOTAL_MEM_GB=$(($(sysctl -n hw.memsize) / 1024 / 1024 / 1024))
  SYSTEM_INFO=$(sysctl -n hw.model)
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
  CPU_CORES=$(nproc)
  TOTAL_MEM_GB=$(($(grep MemTotal /proc/meminfo | awk '{print $2}') / 1024 / 1024))
  SYSTEM_INFO=$(uname -s)
else
  CPU_CORES=4
  TOTAL_MEM_GB=8
  SYSTEM_INFO="Unknown"
fi

usage(){
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  --root DIR         Root folder with sample subfolders (default: "$ROOT_DIR")
  --env PREFIX       Conda env prefix to use (default: "$ENV_PREFIX")
  --markers LIST     Comma-separated markers to search for (default: "$MARKERS")
                     Valid markers: 18S, COI, JEDI
                     Examples: "18S,COI" (water), "JEDI,COI" (soil), "18S,COI,JEDI" (all)
  --threads N        Threads for minimap2/samtools (default: $THREADS)
  --min_reads N      Minimum reads to attempt consensus (default: $MIN_READS)
  --mapq N           MAPQ threshold for grouping reads (default: $MAPQ)
  --keep_percent N   Filtlong keep percent (default: $KEEP_PERCENT)
  --min_length N     Minimum read length filter (default: $MIN_LENGTH, 0 disables)
  --min_mean_q N     Minimum mean read Q-score (default: $MIN_MEAN_Q)
  -h, --help         Show this help

Example:
  # Water eDNA (18S + COI)
  $(basename "$0") --root data/Water_eDNA_18S_COI_14_01_26/fastq_pass --markers 18S,COI
  # Soil eDNA (JEDI + COI)
  $(basename "$0") --root data/Soil_eDNA_JEDI_COI_14_01_26/fastq_pass --markers JEDI,COI
EOF
}

# Function to get current memory usage in MB (macOS)
get_memory_usage() {
  if [[ "$OSTYPE" == "darwin"* ]]; then
    # Get memory used (not wired+active+compressed)
    vm_stat | awk '/Pages active/ {active=$3} /Pages wired/ {wired=$4} /Pages compressed/ {comp=$4} END {printf "%.0f", ((active+wired+comp)*4096)/(1024*1024)}'
  else
    # Linux fallback
    free -m | awk 'NR==2{print $3}'
  fi
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT_DIR="$2"; shift 2;;
    --env) ENV_PREFIX="$2"; shift 2;;
    --markers) MARKERS="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --min_reads) MIN_READS="$2"; shift 2;;
    --mapq) MAPQ="$2"; shift 2;;
    --keep_percent) KEEP_PERCENT="$2"; shift 2;;
    --min_length) MIN_LENGTH="$2"; shift 2;;
    --min_mean_q) MIN_MEAN_Q="$2"; shift 2;;
    -h|--help) usage; exit 0;;
    *) echo "Unknown option: $1" >&2; usage; exit 1;;
  esac
done

mkdir -p out/logs

# Derive run name from the data directory name
RUN_NAME=$(basename "$(dirname "$ROOT_DIR")")
OUTPUT_ROOT="out/$RUN_NAME"

# Create main pipeline log and progress files
PIPELINE_LOG="out/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"
PROGRESS_FILE="out/logs/pipeline_progress.txt"
mkdir -p "$(dirname "$PIPELINE_LOG")"

# Function to log and display messages
log_message() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') | $@" | tee -a "$PIPELINE_LOG"
}

# Function to update progress file
update_progress() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') | $@" >> "$PROGRESS_FILE"
}

# Initialize progress file
echo "Pipeline started: $(date)" > "$PROGRESS_FILE"
echo "Log file: $PIPELINE_LOG" >> "$PROGRESS_FILE"
echo "" >> "$PROGRESS_FILE"

log_message "Root: $ROOT_DIR"
log_message "Run name: $RUN_NAME"
log_message "Output root: $OUTPUT_ROOT"
log_message "Using conda env prefix: $ENV_PREFIX"
log_message "Main log: $PIPELINE_LOG"
log_message ""
log_message "SYSTEM RESOURCES"
log_message "======================================="
log_message "CPU Cores available: $CPU_CORES"
log_message "Total memory: ${TOTAL_MEM_GB} GB"
log_message "System: $SYSTEM_INFO"
log_message ""
log_message "PIPELINE CONFIGURATION"
log_message "======================================="
log_message "Threads/jobs per task: $THREADS"
log_message "Active markers: $MARKERS"
log_message "Quality threshold (min_mean_q): $MIN_MEAN_Q"
log_message "Clustering method: VSEARCH (global sequence alignment)"
log_message "Clustering identity threshold: 0.95"
log_message ""

# Timing tracking - arrays to collect per-barcode stats
declare -a barcode_names
declare -a preproc_times
declare -a total_times
start_time=$(date +%s)

# Initialize global step timing variables
marker_time=0
cluster_time=0
matrix_time=0
taxonomy_time=0
blast_time=0
summary_time=0

# Memory tracking arrays
declare -a step_names
declare -a step_times
declare -a step_memory_start
declare -a step_memory_peak

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
echo ""
echo "======================================="
echo "STARTING FULL PIPELINE PROCESSING"
echo "======================================="
log_message ""
log_message "======================================="
log_message "STARTING FULL PIPELINE PROCESSING"
log_message "======================================="
log_message "Processing all ${#samples[@]} samples with full timing and resource tracking..."
update_progress "Processing all ${#samples[@]} samples..."
echo ""

for sample in "${samples[@]}"; do
  echo
  echo "========== Processing: $sample =========="
  log_message ""
  log_message "========== Processing: $sample =========="
  update_progress "[PREPROC] Starting $sample"
  sample_start=$(date +%s)
  sample_dir="$ROOT_DIR/$sample"
  outdir="$OUTPUT_ROOT/$sample"
  logf="$OUTPUT_ROOT/logs/${sample}.log"

  # Create fresh output directories for this sample
  echo "Creating output directory: $outdir"
  log_message "Creating output directory: $outdir"
  mkdir -p "$outdir" "$OUTPUT_ROOT/logs"
  echo "Sample: $sample - Started at $(date)" > "$logf"
  echo "Output directory: $outdir" >> "$logf"

  fq_list=()
  while IFS= read -r -d $'\0' f; do
    fq_list+=("$f")
  done < <(find "$sample_dir" -type f -name '*.fastq.gz' -print0 | sort -z)

  if [ ${#fq_list[@]} -eq 0 ]; then
    echo "No FASTQ files found under $sample_dir" | tee -a "$logf"
    continue
  fi

  raw_reads_path="$outdir/raw_reads.fastq.gz"
  echo "Concatenating $((${#fq_list[@]})) FASTQ chunks into $raw_reads_path" | tee -a "$logf"
  log_message "Concatenating $((${#fq_list[@]})) FASTQ files for $sample"
  cat "${fq_list[@]}" > "$raw_reads_path"
  raw_reads_mb=$(du -m "$raw_reads_path" | cut -f1)

  step_start=$(date +%s)
  echo "Running preprocessing (quality filter only) ..." | tee -a "$logf"
  log_message "Running preprocessing for $sample (${raw_reads_mb}MB)..."
  export PATH="$ENV_PREFIX/bin:/opt/homebrew/bin:$PATH"
  if ! "$ENV_PREFIX/bin/python3" scripts/1_run_preprocessing.py --input_files "$sample_dir"/*.fastq.gz --output_dir "$outdir" --min_mean_q "$MIN_MEAN_Q" >> "$logf" 2>&1; then
    echo "Preprocessing failed for $sample (see $logf). Continuing to next sample." | tee -a "$logf"
    log_message "✗ Preprocessing FAILED for $sample"
    update_progress "[PREPROC] FAILED $sample"
    continue
  fi
  step_end=$(date +%s)
  preproc_time=$((step_end - step_start))
  echo "[TIMING] $sample - Preprocessing: ${preproc_time}s" | tee -a "$logf"
  log_message "✓ Preprocessing complete for $sample (${preproc_time}s)"
  update_progress "[PREPROC] Complete $sample (${preproc_time}s)"

  filtered="$outdir/filtered_reads.fastq.gz"
  if [ ! -f "$filtered" ]; then
    echo "Filtered file not found for $sample; skipping" | tee -a "$logf"
    log_message "✗ Filtered file not found for $sample"
    update_progress "[PREPROC] No output for $sample"
    continue
  fi
  filtered_reads_mb=$(du -m "$filtered" | cut -f1)

  sample_end=$(date +%s)
  sample_total=$((sample_end - sample_start))
  
  # Collect statistics (clustering timing will be added after global clustering step)
  barcode_names+=("$sample")
  preproc_times+=($preproc_time)
  total_times+=($sample_total)

  echo "[STATS] $sample - Preprocessing: ${preproc_time}s (${raw_reads_mb}MB → ${filtered_reads_mb}MB)" | tee -a "$logf"
  log_message "[STATS] $sample - Preprocessing: ${preproc_time}s (${raw_reads_mb}MB → ${filtered_reads_mb}MB)"
  echo "Finished preprocessing: $sample" | tee -a "$logf"
done

echo ""
echo "======================================="
echo "GLOBAL MARKER CLASSIFICATION & CLUSTERING"
echo "======================================="
echo ""

# Step 2: Classify reads by marker (18S vs COI)
echo "[2/7] Classifying reads by marker..."
log_message ""
log_message "======================================="
log_message "GLOBAL MARKER CLASSIFICATION & CLUSTERING"
log_message "======================================="
log_message "[2/7] Classifying reads by marker..."
update_progress "[MARKER] Starting marker classification..."
marker_start=$(date +%s)
marker_mem_start=$(get_memory_usage)
if ! "$ENV_PREFIX/bin/python3" scripts/2_classify_markers.py \
    --input_dir "$OUTPUT_ROOT" \
    --markers "$MARKERS" > "$OUTPUT_ROOT/logs/marker_classification.log" 2>&1; then
  echo "Marker classification failed (see $OUTPUT_ROOT/logs/marker_classification.log)" >&2
  log_message "✗ Marker classification FAILED"
  exit 1
fi
marker_end=$(date +%s)
marker_mem_end=$(get_memory_usage)
marker_time=$((marker_end - marker_start))
marker_mem_delta=$((marker_mem_end - marker_mem_start))
echo "✓ Marker classification complete (${marker_time}s, ${marker_mem_delta}MB memory)"
log_message "✓ Marker classification complete (${marker_time}s, ${marker_mem_delta}MB memory)"
update_progress "[MARKER] Complete (${marker_time}s)"
echo ""

# Step 3: Run clustering separately for each marker
echo "[3/7] Running global OTU clustering..."
log_message "[3/7] Running global OTU clustering..."
update_progress "[CLUSTERING] Starting global clustering..."
cluster_start=$(date +%s)
cluster_mem_start=$(get_memory_usage)
export PATH="$ENV_PREFIX/bin:/opt/homebrew/bin:$PATH"
if ! "$ENV_PREFIX/bin/python3" scripts/3_run_clustering_by_marker.py \
    --input_dir "$OUTPUT_ROOT" \
    --output_dir "$OUTPUT_ROOT" \
    --identity 0.95 \
    --threads "$THREADS" \
    --markers "$MARKERS" > "$OUTPUT_ROOT/logs/global_clustering.log" 2>&1; then
  echo "Global clustering failed (see $OUTPUT_ROOT/logs/global_clustering.log)" >&2
  log_message "✗ Global clustering FAILED"
  exit 1
fi
cluster_end=$(date +%s)
cluster_mem_end=$(get_memory_usage)
cluster_time=$((cluster_end - cluster_start))
cluster_mem_delta=$((cluster_mem_end - cluster_mem_start))
echo "✓ Global clustering complete (${cluster_time}s, ${cluster_mem_delta}MB memory)"
log_message "✓ Global clustering complete (${cluster_time}s, ${cluster_mem_delta}MB memory)"
update_progress "[CLUSTERING] Complete (${cluster_time}s)"
echo ""

# Step 4: Generate marker-separated abundance matrices
echo "[4/7] Generating abundance matrices..."
log_message "[4/7] Generating abundance matrices..."
update_progress "[MATRICES] Starting abundance matrix generation..."
matrix_start=$(date +%s)
matrix_mem_start=$(get_memory_usage)
if ! "$ENV_PREFIX/bin/python3" scripts/4_merge_otu_tables_by_marker.py \
    --input_dir "$OUTPUT_ROOT" \
    --markers "$MARKERS" > "$OUTPUT_ROOT/logs/abundance_matrices.log" 2>&1; then
  echo "Matrix generation failed (see $OUTPUT_ROOT/logs/abundance_matrices.log)" >&2
  log_message "✗ Matrix generation FAILED"
  exit 1
fi
matrix_end=$(date +%s)
matrix_mem_end=$(get_memory_usage)
matrix_time=$((matrix_end - matrix_start))
matrix_mem_delta=$((matrix_mem_end - matrix_mem_start))
echo "✓ Abundance matrices complete (${matrix_time}s, ${matrix_mem_delta}MB memory)"
log_message "✓ Abundance matrices complete (${matrix_time}s, ${matrix_mem_delta}MB memory)"
update_progress "[MATRICES] Complete (${matrix_time}s)"
echo ""

# Step 5: Assign taxonomy
echo "[5/7] Assigning taxonomy with VSEARCH SINTAX..."
log_message "[5/7] Assigning taxonomy with VSEARCH SINTAX..."
update_progress "[TAXONOMY] Starting taxonomy assignment..."
taxonomy_start=$(date +%s)
taxonomy_mem_start=$(get_memory_usage)

# Build database arguments dynamically based on active markers
TAXONOMY_DB_ARGS=""
IFS=',' read -ra MARKER_ARRAY <<< "$MARKERS"
for m in "${MARKER_ARRAY[@]}"; do
  case "$m" in
    18S) [ -f refs/silva_18s_v123.udb ] && TAXONOMY_DB_ARGS="$TAXONOMY_DB_ARGS --db_18S refs/silva_18s_v123.udb" ;;
    COI) [ -f refs/midori2_COI.udb ] && TAXONOMY_DB_ARGS="$TAXONOMY_DB_ARGS --db_COI refs/midori2_COI.udb" ;;
    JEDI) 
      # JEDI targets COI - use the same COI reference database (MIDORI2/eKOI)
      [ -f refs/midori2_COI.udb ] && TAXONOMY_DB_ARGS="$TAXONOMY_DB_ARGS --db_JEDI refs/midori2_COI.udb"
      ;;
  esac
done

if [ -z "$TAXONOMY_DB_ARGS" ]; then
  echo "WARNING: No reference databases found for active markers ($MARKERS). Skipping taxonomy." >&2
  log_message "⚠ No databases found, skipping taxonomy assignment"
else
  if ! "$ENV_PREFIX/bin/python3" scripts/5_assign_taxonomy.py \
      --input_dir "$OUTPUT_ROOT" \
      $TAXONOMY_DB_ARGS \
      --threads "$THREADS" > "$OUTPUT_ROOT/logs/taxonomy_assignment.log" 2>&1; then
    echo "Taxonomy assignment failed (see $OUTPUT_ROOT/logs/taxonomy_assignment.log)" >&2
    log_message "✗ Taxonomy assignment FAILED"
    exit 1
  fi
fi
taxonomy_end=$(date +%s)
taxonomy_mem_end=$(get_memory_usage)
taxonomy_time=$((taxonomy_end - taxonomy_start))
taxonomy_mem_delta=$((taxonomy_mem_end - taxonomy_mem_start))
echo "✓ Taxonomy assignment complete (${taxonomy_time}s, ${taxonomy_mem_delta}MB memory)"
log_message "✓ Taxonomy assignment complete (${taxonomy_time}s, ${taxonomy_mem_delta}MB memory)"
update_progress "[TAXONOMY] Complete (${taxonomy_time}s)"
echo ""

# Step 6: BLAST validation (optional, top 10 OTUs per active marker)
echo "[6/7] Running BLAST validation (top 10 OTUs)..."
log_message "[6/7] Running BLAST validation (top 10 OTUs)..."
update_progress "[BLAST] Starting BLAST validation..."
blast_start=$(date +%s)
blast_mem_start=$(get_memory_usage)

IFS=',' read -ra BLAST_MARKERS <<< "$MARKERS"
for bm in "${BLAST_MARKERS[@]}"; do
  if [ -f "$OUTPUT_ROOT/merged/otu_relative_abundance_${bm}.csv" ]; then
    echo "  BLASTing top 10 ${bm} OTUs..."
    log_message "  BLASTing top 10 ${bm} OTUs..."
    if "$ENV_PREFIX/bin/python3" scripts/6_blast_top_otus.py \
        --matrix "$OUTPUT_ROOT/merged/otu_relative_abundance_${bm}.csv" \
        --fasta "$OUTPUT_ROOT/temp_clustering/consensus_${bm}_clean.fasta" \
        --otu_assignment "$OUTPUT_ROOT/global_otu_assignment_${bm}.txt" \
        --marker "$bm" \
        --top_n 10 > "$OUTPUT_ROOT/logs/blast_${bm}.log" 2>&1; then
      echo "  ✓ ${bm} BLAST complete"
      log_message "  ✓ ${bm} BLAST complete"
    else
      echo "  ⚠ ${bm} BLAST failed (non-critical, continuing...)"
      log_message "  ⚠ ${bm} BLAST failed (non-critical, continuing...)"
    fi
  fi
done
blast_end=$(date +%s)
blast_mem_end=$(get_memory_usage)
blast_time=$((blast_end - blast_start))
blast_mem_delta=$((blast_mem_end - blast_mem_start))
echo "✓ BLAST validation complete (${blast_time}s, ${blast_mem_delta}MB memory)"
log_message "✓ BLAST validation complete (${blast_time}s, ${blast_mem_delta}MB memory)"
update_progress "[BLAST] Complete (${blast_time}s)"
echo ""

# Step 7: Comprehensive taxonomy summary
echo "[7/7] Generating comprehensive taxonomy summary..."
log_message "[7/7] Generating comprehensive taxonomy summary..."
update_progress "[SUMMARY] Starting comprehensive summary..."
summary_start=$(date +%s)
summary_mem_start=$(get_memory_usage)
if ! "$ENV_PREFIX/bin/python3" scripts/7_comprehensive_taxonomy_summary.py \
    --input_dir "$OUTPUT_ROOT" \
    --markers "$MARKERS" \
    --skip_blast > "$OUTPUT_ROOT/logs/taxonomy_summary.log" 2>&1; then
  echo "Taxonomy summary failed (see $OUTPUT_ROOT/logs/taxonomy_summary.log)" >&2
  log_message "✗ Taxonomy summary FAILED"
  exit 1
fi
summary_end=$(date +%s)
summary_mem_end=$(get_memory_usage)
summary_time=$((summary_end - summary_start))
summary_mem_delta=$((summary_mem_end - summary_mem_start))
echo "✓ Comprehensive summary complete (${summary_time}s, ${summary_mem_delta}MB memory)"
log_message "✓ Comprehensive summary complete (${summary_time}s, ${summary_mem_delta}MB memory)"
update_progress "[SUMMARY] Complete (${summary_time}s)"
echo ""

echo ""
echo "======================================="
echo "FINAL PROCESSING SUMMARY"
echo "======================================="
log_message ""
log_message "======================================="
log_message "FINAL PROCESSING SUMMARY"
log_message "======================================="
update_progress "[COMPLETE] Pipeline finished"
end_time=$(date +%s)
total_time=$((end_time - start_time))
hours=$((total_time / 3600))
mins_rem=$((total_time % 3600 / 60))
secs_rem=$((total_time % 60))

# Calculate total preprocessing time
total_preproc=0
for t in "${preproc_times[@]}"; do
  total_preproc=$((total_preproc + t))
done

# Calculate global pipeline steps time
global_steps_time=$((marker_time + cluster_time + matrix_time + taxonomy_time + blast_time + summary_time))

echo "Total pipeline runtime: ${hours}h ${mins_rem}m ${secs_rem}s (${total_time}s)"
log_message "Total pipeline runtime: ${hours}h ${mins_rem}m ${secs_rem}s (${total_time}s)"
echo ""

# Print per-barcode statistics
if [ ${#barcode_names[@]} -gt 0 ]; then
  echo "COMPUTATIONAL RESOURCES USED"
  echo "======================================="
  echo "CPU Cores allocated: $THREADS per task (of $CPU_CORES available)"
  echo "Total memory available: ${TOTAL_MEM_GB} GB"
  echo "System: $SYSTEM_INFO"
  echo ""
  
  echo "PER-BARCODE PREPROCESSING TIMES"
  echo "======================================="
  printf "%-15s %-15s\n" "Barcode" "Preproc(s)"
  echo "---------------------------------------"
  
  for i in "${!barcode_names[@]}"; do
    printf "%-15s %-15d\n" "${barcode_names[$i]}" "${preproc_times[$i]}"
  done
  
  echo "---------------------------------------"
  num_samples=${#barcode_names[@]}
  avg_preproc=$((total_preproc / num_samples))
  printf "%-15s %-15d\n" "AVERAGE" "$avg_preproc"
  printf "%-15s %-15d\n" "TOTAL" "$total_preproc"
  
  echo ""
  echo "GLOBAL PIPELINE STEPS TIMING"
  echo "======================================="
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "Marker classification" "$marker_time" "$((marker_time * 100 / total_time))" "$marker_mem_delta"
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "Global clustering" "$cluster_time" "$((cluster_time * 100 / total_time))" "$cluster_mem_delta"
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "Abundance matrices" "$matrix_time" "$((matrix_time * 100 / total_time))" "$matrix_mem_delta"
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "Taxonomy assignment" "$taxonomy_time" "$((taxonomy_time * 100 / total_time))" "$taxonomy_mem_delta"
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "BLAST validation" "$blast_time" "$((blast_time * 100 / total_time))" "$blast_mem_delta"
  printf "%-30s %10ds (%3d%%) %+8dMB\n" "Comprehensive summary" "$summary_time" "$((summary_time * 100 / total_time))" "$summary_mem_delta"
  echo "---------------------------------------"
  printf "%-30s %10ds (%3d%%)\n" "Total global steps" "$global_steps_time" "$((global_steps_time * 100 / total_time))"
  
  echo ""
  echo "TIME ALLOCATION ANALYSIS"
  echo "======================================="
  echo "Per-barcode preprocessing: ${total_preproc}s (~$((total_preproc * 100 / total_time))%)"
  echo "Global pipeline steps:     ${global_steps_time}s (~$((global_steps_time * 100 / total_time))%)"
  echo "File I/O & overhead:       $((total_time - total_preproc - global_steps_time))s (~$(((total_time - total_preproc - global_steps_time) * 100 / total_time))%)"
  
  echo ""
  echo "PERFORMANCE NOTES"
  echo "======================================="
  if [ $((cluster_time * 100 / total_time)) -gt 50 ]; then
    echo "⚠️  Global clustering is the bottleneck (>50% of runtime)"
    echo "    Consider: reducing identity threshold or increasing threads"
  elif [ $((total_preproc * 100 / total_time)) -gt 50 ]; then
    echo "⚠️  Preprocessing is the bottleneck (>50% of runtime)"
    echo "    Consider: adjusting quality thresholds"
  else
    echo "✓ Well-balanced resource usage across all pipeline steps"
  fi
  
  echo ""
  echo "OUTPUT FILES"
  echo "======================================="
  echo "Abundance matrices:    $OUTPUT_ROOT/merged/"
  echo "Taxonomy assignments:  $OUTPUT_ROOT/taxonomy/"
  echo "BLAST validation:      $OUTPUT_ROOT/blast_results/"
  echo "Summary CSV:           $OUTPUT_ROOT/taxonomy_summary/"
  echo "Logs:                  $OUTPUT_ROOT/logs/"
fi

echo ""
echo "============================================================"
echo "PIPELINE COMPLETE!"
echo "============================================================"
echo "All done. Results are in: $OUTPUT_ROOT"
log_message ""
log_message "============================================================"
log_message "PIPELINE COMPLETE!"
log_message "============================================================"
log_message "All done. Results are in: $OUTPUT_ROOT"
log_message "Main log file: $PIPELINE_LOG"

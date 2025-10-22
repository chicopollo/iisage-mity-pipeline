#!/usr/bin/env bash
# ===============================================================
# run_mito_pipeline.sh - Master Pipeline Orchestrator
# ---------------------------------------------------------------
# Mitochondrial variant calling pipeline for non-model species
# Author: Paul Ten (Bronikowski Lab, IISAGE)
# Version: 1.1 (2025-10)
# Changed detect_species function
# ===============================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${CONFIG_FILE:-$SCRIPT_DIR/config/pipeline.conf}"
LOG_DIR="$SCRIPT_DIR/logs"

# Source config
if [[ -f "$CONFIG_FILE" ]]; then
  source "$CONFIG_FILE"
else
  echo "ERROR: Config file not found: $CONFIG_FILE"
  exit 1
fi

# ---- Usage ----
usage() {
  cat << EOF
Usage: $(basename "$0") [OPTIONS] [SPECIES...]

Mitochondrial variant calling pipeline for non-model species.

STEPS:
  1. BWA alignment (FASTQ → BAM)
  2. Add read groups (BAM → bam_rg/)
  3. GenBank to GFF3 conversion
  4. MITY variant calling & annotation

OPTIONS:
  -s, --step STEP        Run specific step (1-4) or 'all' (default: all)
  -f, --from-step STEP   Run from step N to end
  --species NAME         Run on specific species (default: all)
  --list                 List available species
  --skip-completed       Skip steps with existing outputs (default: enabled)
  --force                Force re-run all steps
  -h, --help             Show this help

EXAMPLES:
  # Run entire pipeline on all species
  $(basename "$0")

  # Run entire pipeline on one species
  $(basename "$0") --species Chrysemys_picta

  # Run only alignment step
  $(basename "$0") --step 1

  # Resume from variant calling
  $(basename "$0") --from-step 4

  # Force re-run everything
  $(basename "$0") --force

EOF
  exit 0
}

# ---- Utility Functions ----
log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_DIR/pipeline.log"
}

error() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" | tee -a "$LOG_DIR/pipeline.log" >&2
  exit 1
}


detect_species() {
  for dir in */; do
    dir=${dir%/}  # Remove trailing slash
    # Skip special directories
    [[ "$dir" == "logs" ]] && continue
    [[ "$dir" == "scripts" ]] && continue
    [[ "$dir" == "config" ]] && continue
    [[ "$dir" == "tools" ]] && continue
    [[ "$dir" == "."* ]] && continue
    
    # Check if it has ref_mito* subdirectory
    if ls -d "$dir"/ref_mito*/ >/dev/null 2>&1; then
      echo "$dir"
    fi
  done | sort
}


list_species() {
  echo "==> Available species:"
  local count=0
  while IFS= read -r species; do
    ((count++))
    ref_dir=$(find "$species" -maxdepth 1 -type d -name "ref_mito_*" 2>/dev/null | head -n1 || echo "N/A")
    fastq_count=$(find "$species" -name "*_R1_*.fastq.gz" 2>/dev/null | wc -l || echo 0)
    bam_count=$(find "$species/bam_rg" -name "*_with_RG.bam" 2>/dev/null | wc -l || echo 0)
    printf "  %2d. %-30s (ref: %-20s, FASTQ: %2d, BAM_RG: %2d)\n" \
      "$count" "$species" "$(basename "$ref_dir")" "$fastq_count" "$bam_count"
  done < <(detect_species)
  echo
  [[ $count -eq 0 ]] && echo "  No species directories found."
  echo "Total: $count species"
}

check_completed() {
  local species=$1
  local step=$2
  
  case $step in
    1) # BWA alignment
      [[ -n "$(find "$species" -maxdepth 1 -name "*.sorted.bam" 2>/dev/null)" ]]
      ;;
    2) # Read groups
      [[ -d "$species/bam_rg" ]] && [[ -n "$(find "$species/bam_rg" -name "*_with_RG.bam" 2>/dev/null)" ]]
      ;;
    3) # GFF3
      [[ -n "$(find "$species/ref_mito"*/ -name "*.gff3.gz" 2>/dev/null)" ]]
      ;;
    4) # MITY
      [[ -f "$species/mity_output/${species}_all.mity.report.xlsx" ]]
      ;;
    *)
      return 1
      ;;
  esac
}

# ---- Step Functions ----
run_step1_alignment() {
  local species=$1
  log "Step 1: BWA Alignment - $species"
  
  if [[ $SKIP_COMPLETED == "true" ]] && check_completed "$species" 1; then
    log "  ✓ Alignment already completed, skipping"
    return 0
  fi
  
  bash "$SCRIPT_DIR/scripts/01_bwa_align.sh" "$species" 2>&1 | tee -a "$LOG_DIR/${species}_step1.log"
}

run_step2_readgroups() {
  local species=$1
  log "Step 2: Add Read Groups - $species"
  
  if [[ $SKIP_COMPLETED == "true" ]] && check_completed "$species" 2; then
    log "  ✓ Read groups already added, skipping"
    return 0
  fi
  
  bash "$SCRIPT_DIR/scripts/02_add_readgroups.sh" "$species" 2>&1 | tee -a "$LOG_DIR/${species}_step2.log"
}

run_step3_gff3() {
  log "Step 3: Convert GenBank to GFF3"
  
  # This runs on all species at once
  if [[ $SKIP_COMPLETED == "true" ]]; then
    local all_done=true
    while IFS= read -r sp; do
      if ! check_completed "$sp" 3; then
        all_done=false
        break
      fi
    done < <(detect_species)
    
    if [[ $all_done == "true" ]]; then
      log "  ✓ All GFF3 files already exist, skipping"
      return 0
    fi
  fi
  
  bash "$SCRIPT_DIR/scripts/03_genbank_to_gff3.sh" 2>&1 | tee -a "$LOG_DIR/step3_gff3.log"
}

run_step4_mity() {
  local species=$1
  log "Step 4: MITY Variant Calling - $species"
  
  if [[ $SKIP_COMPLETED == "true" ]] && check_completed "$species" 4; then
    log "  ✓ MITY analysis already completed, skipping"
    return 0
  fi
  
  cd "$species"
  bash "$SCRIPT_DIR/scripts/04_mity_pipeline.sh" 2>&1 | tee -a "$LOG_DIR/${species}_step4.log"
  cd ..
}

# ---- Main Pipeline ----
run_pipeline() {
  local species=$1
  local start_step=${2:-1}
  local end_step=${3:-4}
  
  log "=========================================="
  log "Running pipeline for: $species"
  log "Steps: $start_step to $end_step"
  log "=========================================="
  
  for step in $(seq $start_step $end_step); do
    case $step in
      1) run_step1_alignment "$species" ;;
      2) run_step2_readgroups "$species" ;;
      3) 
        # Only run once for all species
        if [[ $step -eq $start_step ]] || [[ ! -f "$LOG_DIR/.gff3_done" ]]; then
          run_step3_gff3
          touch "$LOG_DIR/.gff3_done"
        fi
        ;;
      4) run_step4_mity "$species" ;;
    esac
  done
  
  log "✅ Pipeline completed for $species"
}

# ---- Parse Arguments ----
STEP="all"
FROM_STEP=""
SPECIES_LIST=()
SKIP_COMPLETED="true"

while [[ $# -gt 0 ]]; do
  case $1 in
    -s|--step)
      STEP="$2"
      shift 2
      ;;
    -f|--from-step)
      FROM_STEP="$2"
      shift 2
      ;;
    --species)
      SPECIES_LIST+=("$2")
      shift 2
      ;;
    --list)
      list_species
      exit 0
      ;;
    --skip-completed)
      SKIP_COMPLETED="true"
      shift
      ;;
    --force)
      SKIP_COMPLETED="false"
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      SPECIES_LIST+=("$1")
      shift
      ;;
  esac
done

# ---- Setup ----
mkdir -p "$LOG_DIR"
log "==> Mitochondrial Pipeline Starting"
log "Config: $CONFIG_FILE"
log "Logs: $LOG_DIR"

# Determine species to process
if [[ ${#SPECIES_LIST[@]} -eq 0 ]]; then
  mapfile -t SPECIES_LIST < <(detect_species)
fi

if [[ ${#SPECIES_LIST[@]} -eq 0 ]]; then
  error "No species found to process"
fi

log "Species to process: ${SPECIES_LIST[*]}"

# Determine step range
if [[ -n "$FROM_STEP" ]]; then
  START_STEP=$FROM_STEP
  END_STEP=4
elif [[ "$STEP" == "all" ]]; then
  START_STEP=1
  END_STEP=4
else
  START_STEP=$STEP
  END_STEP=$STEP
fi

# ---- Execute ----
FAILED=()
SUCCESSFUL=()

# Clean GFF3 marker for this run
rm -f "$LOG_DIR/.gff3_done"

for species in "${SPECIES_LIST[@]}"; do
  if [[ ! -d "$species" ]]; then
    log "WARNING: Species directory not found: $species"
    FAILED+=("$species")
    continue
  fi
  
  if run_pipeline "$species" "$START_STEP" "$END_STEP"; then
    SUCCESSFUL+=("$species")
  else
    FAILED+=("$species")
  fi
done

# ---- Summary ----
echo
log "=========================================="
log "PIPELINE SUMMARY"
log "=========================================="
log "Total species: ${#SPECIES_LIST[@]}"
log "Successful: ${#SUCCESSFUL[@]}"
log "Failed: ${#FAILED[@]}"

if [[ ${#SUCCESSFUL[@]} -gt 0 ]]; then
  log ""
  log "✅ Successful species:"
  printf '%s\n' "${SUCCESSFUL[@]}" | sed 's/^/   - /' | tee -a "$LOG_DIR/pipeline.log"
fi

if [[ ${#FAILED[@]} -gt 0 ]]; then
  log ""
  log "❌ Failed species:"
  printf '%s\n' "${FAILED[@]}" | sed 's/^/   - /' | tee -a "$LOG_DIR/pipeline.log"
  exit 1
fi

log "=========================================="
log "Pipeline completed successfully!"

#!/usr/bin/env bash
# =============================================================================
# download_data.sh
# Downloads GSE135893 and GSE136831 from NCBI GEO FTP.
#
# Usage:
#   bash data/download_data.sh               # download both datasets (default)
#   bash data/download_data.sh GSE135893     # download one dataset only
#   bash data/download_data.sh --dry-run     # print what would be downloaded
#
# Requirements: wget or curl, and either tar or gzip (standard on Linux/macOS)
# Estimated download size: ~11 GB total
# =============================================================================

set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/raw"
LOG_FILE="${SCRIPT_DIR}/download.log"

GEO_FTP_BASE="https://ftp.ncbi.nlm.nih.gov/geo/series"

# GEO FTP paths follow a pattern: GSE135893 → GSE135nnn/GSE135893
declare -A DATASETS=(
    ["GSE135893"]="GSE135nnn/GSE135893"
    ["GSE136831"]="GSE136nnn/GSE136831"
)

# ── Argument parsing ───────────────────────────────────────────────────────────

DRY_RUN=false
TARGETS=()

for arg in "$@"; do
    case "$arg" in
        --dry-run) DRY_RUN=true ;;
        GSE*) TARGETS+=("$arg") ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# Default: download both
if [ ${#TARGETS[@]} -eq 0 ]; then
    TARGETS=("GSE135893" "GSE136831")
fi

# ── Helper functions ───────────────────────────────────────────────────────────

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"
}

download_file() {
    local url="$1"
    local dest="$2"

    if [ "$DRY_RUN" = true ]; then
        log "DRY RUN: would download $url → $dest"
        return
    fi

    if [ -f "$dest" ]; then
        log "Already exists, skipping: $dest"
        return
    fi

    log "Downloading: $url"
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$dest" "$url"
    elif command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$dest" "$url"
    else
        log "ERROR: Neither wget nor curl found. Please install one."
        exit 1
    fi
}

check_disk_space() {
    # Warn if less than 30 GB free
    local required_gb=30
    local available_gb

    if command -v df &>/dev/null; then
        available_gb=$(df -BG "$SCRIPT_DIR" | awk 'NR==2 {gsub("G",""); print $4}')
        if [ "$available_gb" -lt "$required_gb" ] 2>/dev/null; then
            log "WARNING: Only ${available_gb}GB free. Recommend at least ${required_gb}GB."
            log "         Consider downloading to an external drive or cloud storage."
            read -r -p "Continue anyway? [y/N] " response
            [[ "$response" =~ ^[Yy]$ ]] || exit 0
        fi
    fi
}

# ── Main download logic ────────────────────────────────────────────────────────

log "=== IPF scRNA-seq Data Download ==="
log "Targets: ${TARGETS[*]}"
log "Output directory: $DATA_DIR"

check_disk_space

mkdir -p "$DATA_DIR"

for GSE_ID in "${TARGETS[@]}"; do

    if [[ -z "${DATASETS[$GSE_ID]+_}" ]]; then
        log "ERROR: Unknown dataset $GSE_ID. Valid options: ${!DATASETS[*]}"
        exit 1
    fi

    FTP_PATH="${DATASETS[$GSE_ID]}"
    DEST_DIR="${DATA_DIR}/${GSE_ID}"
    mkdir -p "$DEST_DIR"

    log "──────────────────────────────────────"
    log "Processing $GSE_ID → $DEST_DIR"

    # ── Download the supplementary files tarball ──────────────────────────────
    # GEO stores processed count matrices as supplementary files
    SUPPL_URL="${GEO_FTP_BASE}/${FTP_PATH}/suppl/${GSE_ID}_RAW.tar"
    SUPPL_TAR="${DEST_DIR}/${GSE_ID}_RAW.tar"

    download_file "$SUPPL_URL" "$SUPPL_TAR"

    if [ "$DRY_RUN" = false ] && [ -f "$SUPPL_TAR" ]; then
        log "Extracting $SUPPL_TAR ..."
        tar -xf "$SUPPL_TAR" -C "$DEST_DIR"
        log "Extraction complete."

        # Decompress any .gz files inside
        find "$DEST_DIR" -name "*.gz" | while read -r gz_file; do
            log "Decompressing: $gz_file"
            gunzip -f "$gz_file"
        done
    fi

    # ── Download the series matrix (metadata) ─────────────────────────────────
    MATRIX_URL="${GEO_FTP_BASE}/${FTP_PATH}/matrix/${GSE_ID}_series_matrix.txt.gz"
    MATRIX_FILE="${DEST_DIR}/${GSE_ID}_series_matrix.txt.gz"

    download_file "$MATRIX_URL" "$MATRIX_FILE"

    if [ "$DRY_RUN" = false ] && [ -f "$MATRIX_FILE" ]; then
        gunzip -f "$MATRIX_FILE"
        log "Series matrix downloaded and decompressed."
    fi

    log "$GSE_ID download complete."

done

log "══════════════════════════════════════"
log "All downloads complete."
log ""
log "Next steps:"
log "  1. Activate your conda environment:  conda activate ipf-scrna"
log "  2. Launch Jupyter:                   jupyter lab"
log "  3. Run notebooks in order starting with 01_qc_preprocessing.ipynb"
log ""
log "If running in Colab or a cloud VM, see data/README.md for"
log "alternative setup instructions."
log "══════════════════════════════════════"

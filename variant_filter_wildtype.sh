#!/bin/bash
#SBATCH --job-name=sbw_filter_v2
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=filter_output_v2_%j.out
#SBATCH --error=filter_output_v2_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR

module load biocontainers
module load bcftools

IN_DIR="${IN_DIR:-${PROJECT_DIR}/06.1.variant_calling_wildtype}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/06.1.variant_calling_wildtype}"

mkdir -p "$OUT_DIR"

IN_BCF="${IN_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/wildtype_diapausing_variants.bcf}"
SITE_BCF="${SITE_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/wildtype_site_filtered.bcf}"
GEMMA_BCF="${GEMMA_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/wildtype_GEMMA_filtered.bcf}"

echo "========================================================="
echo "Starting v3 site and GEMMA variant filtering"
echo "Input: $IN_BCF"
echo "Output: $GEMMA_BCF"
echo
echo "Site filters:"
echo "  QUAL >= 30"
echo "  MQ >= 40"
echo
echo "GEMMA filters:"
echo "  PASS only"
echo "  F_MISSING <= 0.10"
echo "  Biallelic sites only"
echo "  MAF >= 0.05"
echo "  MAC >= 2"
echo "========================================================="
date

# Verify input exists
if [[ ! -f "$IN_BCF" ]]; then
    echo "ERROR: Input BCF not found:"
    echo "$IN_BCF"
    exit 1
fi

echo "Step 1: Applying site-quality filters..."

bcftools filter --threads 4 \
    -e 'QUAL < 30 || MQ < 40' \
    -s 'LOW_QUAL' \
    -Ob \
    -o "$SITE_BCF" \
    "$IN_BCF"

bcftools index -f "$SITE_BCF"

echo "Step 1 complete."
echo

echo "Step 2: Applying GEMMA variant filters..."

bcftools view --threads 4 \
    -f PASS \
    -i 'F_MISSING <= 0.10 && MAC >= 2' \
    -m 2 -M 2 \
    -q 0.05:minor \
    -Ob \
    -o "$GEMMA_BCF" \
    "$SITE_BCF"

bcftools index -f "$GEMMA_BCF"

echo "Step 2 complete."
echo

echo "========================================================="
echo "Filtering complete."
echo "Output files:"
ls -lh "$SITE_BCF" "$SITE_BCF.csi"
ls -lh "$GEMMA_BCF" "$GEMMA_BCF.csi"
echo "========================================================="

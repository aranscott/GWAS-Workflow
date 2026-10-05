#!/bin/bash
#SBATCH --job-name=bcf_parallel_v3
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=variant_call_parallel_v3_%j.out
#SBATCH --error=variant_call_parallel_v3_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR

# Load modules
module load biocontainers
module load bcftools
module load parallel

# Paths
REF_FASTA="${REF_FASTA:-${REFERENCE_DIR}/GCA_025370935.1_NRCan_CFum_1_genomic.fna}"

# Corrected 210-sample BAM list
BAM_LIST="${BAM_LIST:-${PROJECT_DIR}/corrected_210_bams.txt}"

# New output directory
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/06.1.variant_calling_merged_v3}"

mkdir -p "$OUT_DIR"

echo "========================================================="
echo "Starting v3 population variant calling"
echo "Reference: $REF_FASTA"
echo "BAM list:  $BAM_LIST"
echo "Output:    $OUT_DIR"
echo "========================================================="
date

# Verify BAM list
echo
echo "BAM count:"
wc -l "$BAM_LIST"

BAM_COUNT=$(wc -l < "$BAM_LIST")

if [[ "$BAM_COUNT" -ne 210 ]]; then
    echo "ERROR: Expected 210 BAMs but found $BAM_COUNT"
    exit 1
fi

echo "210 BAMs confirmed."

# Verify that every BAM in the list exists
echo
echo "Checking BAM files..."

MISSING=0

while read -r BAM
do
    if [[ ! -f "$BAM" ]]; then
        echo "MISSING BAM: $BAM"
        MISSING=$((MISSING + 1))
    fi
done < "$BAM_LIST"

if [[ "$MISSING" -ne 0 ]]; then
    echo "ERROR: $MISSING BAM files are missing."
    exit 1
fi

echo "All 210 BAMs found."

# Generate reference sequence/scaffold list
cut -f1 "${REF_FASTA}.fai" > "$OUT_DIR/scaffolds.txt"

echo
echo "Reference sequences:"
wc -l "$OUT_DIR/scaffolds.txt"

echo
echo "Starting parallel variant calling..."
date

parallel --jobs 24 \
    "module load biocontainers bcftools; \
     bcftools mpileup \
         -r {} \
         -f $REF_FASTA \
         -b $BAM_LIST \
         -a FORMAT/AD,FORMAT/DP \
         -Ou \
     | bcftools call \
         -mv \
         -Ob \
         -o ${OUT_DIR}/split_{}.bcf" \
    :::: "$OUT_DIR/scaffolds.txt"

echo
echo "Parallel variant calling complete."

# Count split BCFs
SPLIT_COUNT=$(find "$OUT_DIR" \
    -maxdepth 1 \
    -type f \
    -name 'split_*.bcf' \
    | wc -l)

echo "Split BCFs generated: $SPLIT_COUNT"

echo
echo "Concatenating split BCFs..."
date

bcftools concat \
    --threads 24 \
    -O b \
    -o "$OUT_DIR/spruce_budworm_variants.bcf" \
    "$OUT_DIR"/split_*.bcf

echo
echo "Indexing master BCF..."

bcftools index -f \
    "$OUT_DIR/spruce_budworm_variants.bcf"

echo
echo "Removing intermediate split BCFs..."

rm "$OUT_DIR"/split_*.bcf
rm "$OUT_DIR/scaffolds.txt"

echo
echo "========================================================="
echo "Variant calling v3 completed successfully"
echo "Output:"
echo "$OUT_DIR/spruce_budworm_variants.bcf"
echo "========================================================="
date

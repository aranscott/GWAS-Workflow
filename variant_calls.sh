#!/bin/bash

#SBATCH --job-name=bcf_parallel_wt
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=variant_call_wt_%j.out
#SBATCH --error=variant_call_wt_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR


# =========================================================
# LOAD MODULES
# =========================================================

module load biocontainers
module load bcftools
module load parallel


# =========================================================
# PATHS
# =========================================================

REF_FASTA="${REF_FASTA:-${REFERENCE_DIR}/GCA_025370935.1_NRCan_CFum_1_genomic.fna}"

BAM_LIST="${BAM_LIST:-${PROJECT_DIR}/wildtype_195_bams.txt}"

OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/06.1.variant_calling_wildtype}"

mkdir -p "$OUT_DIR"


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "Starting wild-type population variant calling"
echo "========================================================="
echo
echo "Reference: $REF_FASTA"
echo "BAM list:  $BAM_LIST"
echo "Output:    $OUT_DIR"
echo
echo "CPUs:      $SLURM_CPUS_PER_TASK"
echo "Memory:    128G"
echo "Parallel jobs: 24"
echo "========================================================="
date


# =========================================================
# VERIFY REFERENCE
# =========================================================

echo
echo "Checking reference FASTA..."

if [[ ! -f "$REF_FASTA" ]]; then
    echo "ERROR: Reference FASTA not found:"
    echo "$REF_FASTA"
    exit 1
fi

if [[ ! -f "${REF_FASTA}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found:"
    echo "${REF_FASTA}.fai"
    exit 1
fi

echo "Reference FASTA and index found."


# =========================================================
# VERIFY BAM LIST
# =========================================================

echo
echo "BAM count:"

wc -l "$BAM_LIST"

BAM_COUNT=$(wc -l < "$BAM_LIST")

if [[ "$BAM_COUNT" -ne 195 ]]; then
    echo "ERROR: Expected 195 BAMs but found $BAM_COUNT"
    exit 1
fi

echo "195 BAMs confirmed."


# =========================================================
# VERIFY BAM FILES
# =========================================================

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

echo "All 195 BAMs found."


# =========================================================
# VERIFY BAM INDEX FILES
# =========================================================

echo
echo "Checking BAM index files..."

MISSING_INDEX=0

while read -r BAM
do
    if [[ ! -f "${BAM}.bai" ]]; then
        echo "MISSING BAI: ${BAM}.bai"
        MISSING_INDEX=$((MISSING_INDEX + 1))
    fi
done < "$BAM_LIST"

if [[ "$MISSING_INDEX" -ne 0 ]]; then
    echo "ERROR: $MISSING_INDEX BAM index files are missing."
    exit 1
fi

echo "All 195 BAM index files found."


# =========================================================
# CREATE COMPLETE REFERENCE SEQUENCE LIST
#
# IMPORTANT:
# Use every sequence present in the FASTA index.
#
# This includes:
#   nuclear chromosomes
#   mitochondrial sequence
#   any other sequences present in the reference
#
# There is intentionally NO chromosome filter here.
# =========================================================

echo
echo "Creating complete reference sequence list..."

cut -f1 "${REF_FASTA}.fai" \
> "$OUT_DIR/reference_sequences.txt"

SEQUENCE_COUNT=$(wc -l < "$OUT_DIR/reference_sequences.txt")

echo
echo "Reference sequences found: $SEQUENCE_COUNT"

cat "$OUT_DIR/reference_sequences.txt"


if [[ "$SEQUENCE_COUNT" -eq 0 ]]; then
    echo "ERROR: No reference sequences found in FASTA index."
    exit 1
fi


# =========================================================
# START PARALLEL VARIANT CALLING
# =========================================================

echo
echo "Starting parallel variant calling..."
echo "Maximum simultaneous jobs: 24"
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
    :::: "$OUT_DIR/reference_sequences.txt"


echo
echo "Parallel variant calling complete."
date


# =========================================================
# COUNT SPLIT BCFs
# =========================================================

SPLIT_COUNT=$(find "$OUT_DIR" \
    -maxdepth 1 \
    -type f \
    -name 'split_*.bcf' \
    | wc -l)

echo
echo "Split BCFs generated: $SPLIT_COUNT"

if [[ "$SPLIT_COUNT" -ne "$SEQUENCE_COUNT" ]]; then
    echo
    echo "ERROR: Expected $SEQUENCE_COUNT split BCFs but found $SPLIT_COUNT"
    exit 1
fi


# =========================================================
# VERIFY SPLIT BCFs
# =========================================================

echo
echo "Checking split BCF files..."

while read -r SEQ
do
    BCF="${OUT_DIR}/split_${SEQ}.bcf"

    if [[ ! -f "$BCF" ]]; then
        echo "MISSING BCF: $BCF"
        exit 1
    fi
done < "$OUT_DIR/reference_sequences.txt"

echo "All split BCFs found."


# =========================================================
# CONCATENATE ALL REFERENCE SEQUENCES
#
# This combines nuclear + mitochondrial + any other
# sequences represented in reference_sequences.txt.
# =========================================================

echo
echo "Concatenating split BCFs..."
date

bcftools concat \
    --threads 24 \
    -O b \
    -o "$OUT_DIR/wildtype_diapausing_variants.bcf" \
    "$OUT_DIR"/split_*.bcf

echo
echo "Concatenation complete."
date


# =========================================================
# INDEX MASTER BCF
# =========================================================

echo
echo "Indexing master BCF..."

bcftools index -f \
    "$OUT_DIR/wildtype_diapausing_variants.bcf"

echo
echo "Master BCF indexed."


# =========================================================
# REMOVE INTERMEDIATE FILES
# =========================================================

echo
echo "Removing intermediate split BCFs..."

rm "$OUT_DIR"/split_*.bcf
rm "$OUT_DIR/reference_sequences.txt"


# =========================================================
# FINAL CHECK
# =========================================================

if [[ ! -f "$OUT_DIR/wildtype_diapausing_variants.bcf" ]]; then
    echo
    echo "ERROR: Final BCF was not created."
    exit 1
fi


if [[ ! -f "$OUT_DIR/wildtype_diapausing_variants.bcf.csi" &&
      ! -f "$OUT_DIR/wildtype_diapausing_variants.bcf.bai" ]]; then

    echo
    echo "ERROR: Final BCF index was not created."
    exit 1
fi


# =========================================================
# COMPLETE
# =========================================================

echo
echo "========================================================="
echo "Wild-type variant calling completed successfully"
echo "========================================================="
echo
echo "Samples:          195"
echo "Reference seqs:   $SEQUENCE_COUNT"
echo "Parallel jobs:    24"
echo
echo "Final BCF:"
echo "$OUT_DIR/wildtype_diapausing_variants.bcf"
echo
echo "Index:"
ls -lh \
    "$OUT_DIR"/wildtype_diapausing_variants.bcf*
echo
echo "========================================================="
date

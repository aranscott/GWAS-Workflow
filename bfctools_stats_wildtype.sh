#!/bin/bash

#SBATCH --job-name=bcf_stats
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=bcf_stats_%j.out
#SBATCH --error=bcf_stats_%j.err
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


# =========================================================
# PATHS
# =========================================================

BASE_DIR="${BASE_DIR:-${PROJECT_DIR}/06.1.variant_calling_wildtype}"

STATS_DIR=${BASE_DIR}/stats

FILE1=${BASE_DIR}/bcf_filt_qualfilt_monobiallele_snponly_max60.bcf

FILE2=${BASE_DIR}/SBW_whole_genome_qualfilt_monobiallele_snponly_max60.recode.vcf.gz


# =========================================================
# OUTPUT FILES
# =========================================================

STATS1=${STATS_DIR}/bcf_filt_qualfilt_monobiallele_snponly_max60.stats.txt

STATS2=${STATS_DIR}/SBW_whole_genome_qualfilt_monobiallele_snponly_max60.recode.stats.txt

INTERSECTION=${STATS_DIR}/bcftools_stats_intersection.txt


# =========================================================
# CREATE STATS DIRECTORY
# =========================================================

mkdir -p "$STATS_DIR"


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "bcftools statistics and comparison"
echo "========================================================="
echo
echo "Dataset 1:"
echo "$FILE1"
echo
echo "Dataset 2:"
echo "$FILE2"
echo
echo "Stats output directory:"
echo "$STATS_DIR"
echo
date
echo "========================================================="


# =========================================================
# VERIFY INPUT FILES
# =========================================================

for f in "$FILE1" "$FILE2"
do

    if [[ ! -f "$f" ]]; then

        echo
        echo "ERROR: Input file not found:"
        echo "$f"

        exit 1

    fi

done


echo
echo "Both input files found."


# =========================================================
# DATASET 1
# =========================================================

echo
echo "========================================================="
echo "Running bcftools stats on dataset 1"
echo "========================================================="
date

bcftools stats \
    --threads 4 \
    "$FILE1" \
    > "$STATS1"

echo
echo "Dataset 1 stats complete."


# =========================================================
# DATASET 2
# =========================================================

echo
echo "========================================================="
echo "Running bcftools stats on dataset 2"
echo "========================================================="
date

bcftools stats \
    --threads 4 \
    "$FILE2" \
    > "$STATS2"

echo
echo "Dataset 2 stats complete."


# =========================================================
# TWO-FILE COMPARISON
# =========================================================

echo
echo "========================================================="
echo "Running two-file bcftools comparison"
echo "========================================================="
date

bcftools stats \
    --threads 4 \
    "$FILE1" \
    "$FILE2" \
    > "$INTERSECTION"

echo
echo "Two-file comparison complete."


# =========================================================
# BASIC FILE CHECKS
# =========================================================

echo
echo "========================================================="
echo "Output files"
echo "========================================================="

ls -lh \
    "$STATS1" \
    "$STATS2" \
    "$INTERSECTION"


# =========================================================
# DISPLAY SUMMARY STATISTICS
# =========================================================

echo
echo "========================================================="
echo "Dataset 1 summary"
echo "========================================================="

grep '^SN' "$STATS1"


echo
echo "========================================================="
echo "Dataset 2 summary"
echo "========================================================="

grep '^SN' "$STATS2"


echo
echo "========================================================="
echo "Two-file comparison summary"
echo "========================================================="

grep '^SN' "$INTERSECTION"


# =========================================================
# COMPLETE
# =========================================================

echo
echo "========================================================="
echo "bcftools statistics completed successfully"
echo "========================================================="
echo
echo "Stats directory:"
echo "$STATS_DIR"
echo
echo "Individual stats:"
echo "$STATS1"
echo "$STATS2"
echo
echo "Comparison:"
echo "$INTERSECTION"
echo
date

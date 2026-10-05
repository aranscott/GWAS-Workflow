#!/bin/bash

#SBATCH --job-name=filter_wildtype
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=filter_wildtype_%j.out
#SBATCH --error=filter_wildtype_%j.err
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

IN_BCF="${IN_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/wildtype_diapausing_variants.bcf}"

SITE_BCF="${SITE_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/wildtype_site_filtered.bcf}"

FINAL_BCF="${FINAL_BCF:-${PROJECT_DIR}/06.1.variant_calling_wildtype/bcf_filt_qualfilt_monobiallele_snponly_max60.bcf}"


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "Starting wild-type variant filtering"
echo "========================================================="
echo
echo "Input:"
echo "$IN_BCF"
echo
echo "Site-filtered output:"
echo "$SITE_BCF"
echo
echo "Final SNP output:"
echo "$FINAL_BCF"
echo
echo "Filters:"
echo "  QUAL >= 30"
echo "  MQ >= 40"
echo "  PASS only"
echo "  Biallelic sites only"
echo "  SNPs only"
echo "  Missingness <= 60%"
echo "  NO MAF FILTER"
echo "  NO MAC FILTER"
echo "========================================================="
date


# =========================================================
# VERIFY INPUT
# =========================================================

if [[ ! -f "$IN_BCF" ]]; then

    echo
    echo "ERROR: Input BCF not found:"
    echo "$IN_BCF"

    exit 1

fi


# =========================================================
# CHECK SAMPLE COUNT
# =========================================================

echo
echo "Checking sample count..."

SAMPLE_COUNT=$(
    bcftools query \
        -l "$IN_BCF" \
    | wc -l
)

echo "Samples in input BCF: $SAMPLE_COUNT"

if [[ "$SAMPLE_COUNT" -ne 195 ]]; then

    echo
    echo "ERROR: Expected 195 samples but found $SAMPLE_COUNT"

    exit 1

fi

echo "195 samples confirmed."


# =========================================================
# STEP 1: SITE QUALITY FILTERING
#
# QUAL >= 30
# MQ   >= 40
#
# Sites failing either criterion are marked LOW_QUAL.
# They are not yet physically removed.
# =========================================================

echo
echo "Step 1: Applying site-quality filters..."
date

bcftools filter \
    --threads 4 \
    -e 'QUAL < 30 || MQ < 40' \
    -s LOW_QUAL \
    -Ob \
    -o "$SITE_BCF" \
    "$IN_BCF"

bcftools index -f "$SITE_BCF"

echo
echo "Step 1 complete."


# =========================================================
# STEP 2: FINAL SNP FILTERING
#
# PASS
# Biallelic
# SNP only
# Missingness <= 60%
#
# NO MAF FILTER
# NO MAC FILTER
# =========================================================

echo
echo "Step 2: Applying final SNP filters..."
date

bcftools view \
    --threads 4 \
    -f PASS \
    -v snps \
    -m 2 \
    -M 2 \
    -i 'F_MISSING <= 0.60' \
    -Ob \
    -o "$FINAL_BCF" \
    "$SITE_BCF"

bcftools index -f "$FINAL_BCF"

echo
echo "Step 2 complete."


# =========================================================
# FINAL VARIANT COUNT
# =========================================================

echo
echo "Final variant count:"

bcftools view \
    -H \
    "$FINAL_BCF" \
| wc -l


# =========================================================
# FINAL VARIANT COUNT BY CHROMOSOME
# =========================================================

echo
echo "Final SNP count by chromosome:"

bcftools query \
    -f '%CHROM\n' \
    "$FINAL_BCF" \
| sort \
| uniq -c


# =========================================================
# FINAL SAMPLE COUNT
# =========================================================

echo
echo "Final sample count:"

bcftools query \
    -l "$FINAL_BCF" \
| wc -l


# =========================================================
# COMPLETE
# =========================================================

echo
echo "========================================================="
echo "Wild-type variant filtering completed successfully"
echo "========================================================="
echo
echo "Input BCF:"
echo "$IN_BCF"
echo
echo "Site-filtered BCF:"
echo "$SITE_BCF"
echo
echo "Final SNP BCF:"
echo "$FINAL_BCF"
echo
echo "Final BCF files:"
ls -lh "$FINAL_BCF"*
echo
echo "========================================================="
date

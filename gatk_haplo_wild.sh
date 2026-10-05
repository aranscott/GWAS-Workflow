#!/bin/bash

#SBATCH --job-name=hc_gvcf
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=72:00:00
#SBATCH --output=haplocaller_%A_%a.out
#SBATCH --error=haplocaller_%A_%a.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --array=1-195%15

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR


# =========================================================
# LOAD MODULES
# =========================================================

module load biocontainers
module load gatk4/4.2.6.0


# =========================================================
# PATHS
# =========================================================

REF="${REF:-${REF_FASTA:-${REFERENCE_DIR}/GCA_025370935.1_NRCan_CFum_1_genomic.fna}}"

INDIR="${INDIR:-${PROJECT_DIR}/05.marked_duplicates}"

OUTDIR="${OUTDIR:-${PROJECT_DIR}/06.1.variant_calling_wildtype/gatk/gvcfs}"

mkdir -p "$OUTDIR"


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "GATK HaplotypeCaller GVCF"
echo "========================================================="
echo
echo "Reference:"
echo "$REF"
echo
echo "Input BAM directory:"
echo "$INDIR"
echo
echo "Output directory:"
echo "$OUTDIR"
echo
echo "Array task:"
echo "$SLURM_ARRAY_TASK_ID"
echo
echo "CPUs:"
echo "$SLURM_CPUS_PER_TASK"
echo
echo "Memory:"
echo "32G"
echo "========================================================="
date


# =========================================================
# VERIFY REFERENCE
# =========================================================

if [[ ! -f "$REF" ]]; then

    echo
    echo "ERROR: Reference FASTA not found:"
    echo "$REF"

    exit 1

fi


if [[ ! -f "${REF}.fai" ]]; then

    echo
    echo "ERROR: Reference FASTA index not found:"
    echo "${REF}.fai"

    exit 1

fi


if [[ ! -f "${REF%.fna}.dict" ]]; then

    echo
    echo "ERROR: Reference sequence dictionary not found:"
    echo "${REF%.fna}.dict"

    exit 1

fi


echo
echo "Reference FASTA, index, and dictionary found."


# =========================================================
# COLLECT MARKED BAM FILES
#
# Expected files:
#
#   *_marked.bam
#   *_marked.bai
# =========================================================

mapfile -t bams < <(
    find "$INDIR" \
        -maxdepth 1 \
        -type f \
        -name '*_marked.bam' \
        | sort
)


echo
echo "Marked BAM files detected:"
echo "${#bams[@]}"


# =========================================================
# VERIFY SAMPLE COUNT
# =========================================================

if [[ ${#bams[@]} -ne 195 ]]; then

    echo
    echo "ERROR: Expected 195 marked BAM files but found ${#bams[@]}."

    exit 1

fi


echo
echo "195 marked BAMs confirmed."


# =========================================================
# VERIFY BAM INDEX FILES
#
# MarkDuplicates created:
#
#   sample_marked.bam
#   sample_marked.bai
#
# NOT:
#
#   sample_marked.bam.bai
# =========================================================

echo
echo "Checking BAM index files..."

MISSING_INDEX=0


for bam in "${bams[@]}"
do

    bai="${bam%.bam}.bai"

    if [[ ! -f "$bai" ]]; then

        echo "MISSING BAI: $bai"

        MISSING_INDEX=$((MISSING_INDEX + 1))

    fi

done


if [[ "$MISSING_INDEX" -ne 0 ]]; then

    echo
    echo "ERROR: $MISSING_INDEX BAM index files are missing."

    exit 1

fi


echo
echo "All 195 BAM index files found."


# =========================================================
# SELECT BAM FOR THIS ARRAY TASK
# =========================================================

bam="${bams[$((SLURM_ARRAY_TASK_ID - 1))]}"


if [[ -z "${bam:-}" ]]; then

    echo
    echo "ERROR: No BAM found for array task:"
    echo "$SLURM_ARRAY_TASK_ID"

    exit 1

fi


# =========================================================
# VERIFY SELECTED BAM INDEX
# =========================================================

bai="${bam%.bam}.bai"


if [[ ! -f "$bai" ]]; then

    echo
    echo "ERROR: BAM index not found:"
    echo "$bai"

    exit 1

fi


# =========================================================
# DERIVE SAMPLE ID
# =========================================================

sampleID=$(basename "$bam")

sampleID="${sampleID%_marked.bam}"


# =========================================================
# DEFINE OUTPUT
# =========================================================

gvcf="${OUTDIR}/${sampleID}.g.vcf.gz"


# =========================================================
# SKIP COMPLETED SAMPLE
# =========================================================

if [[ -f "$gvcf" && -f "${gvcf}.tbi" ]]; then

    echo
    echo "========================================================="
    echo "GVCF already exists and is indexed"
    echo "========================================================="
    echo
    echo "Sample:"
    echo "$sampleID"
    echo
    echo "GVCF:"
    echo "$gvcf"
    echo
    echo "Skipping sample."
    echo
    date

    exit 0

fi


# =========================================================
# SAMPLE INFORMATION
# =========================================================

echo
echo "========================================================="
echo "HaplotypeCaller"
echo "========================================================="
echo
echo "Array task:"
echo "$SLURM_ARRAY_TASK_ID"
echo
echo "Sample:"
echo "$sampleID"
echo
echo "Input BAM:"
echo "$bam"
echo
echo "Input BAI:"
echo "$bai"
echo
echo "Output GVCF:"
echo "$gvcf"
echo
echo "CPUs:"
echo "$SLURM_CPUS_PER_TASK"
echo
echo "Memory:"
echo "32G"
echo "========================================================="
date


# =========================================================
# RUN HAPLOTYPECALLER
#
# GVCF mode is used so that these 195 samples can
# subsequently be jointly genotyped.
# =========================================================

gatk --java-options "-Xmx28g" HaplotypeCaller \
    -R "$REF" \
    -I "$bam" \
    -O "$gvcf" \
    -ERC GVCF \
    --native-pair-hmm-threads 4


# =========================================================
# VERIFY GVCF
# =========================================================

if [[ ! -f "$gvcf" ]]; then

    echo
    echo "ERROR: HaplotypeCaller did not produce:"
    echo "$gvcf"

    exit 1

fi


echo
echo "HaplotypeCaller completed."


# =========================================================
# VERIFY GVCF INDEX
# =========================================================

if [[ ! -f "${gvcf}.tbi" ]]; then

    echo
    echo "Indexing GVCF..."

    gatk IndexFeatureFile \
        -I "$gvcf"

fi


if [[ ! -f "${gvcf}.tbi" ]]; then

    echo
    echo "ERROR: GVCF index was not created:"
    echo "${gvcf}.tbi"

    exit 1

fi


# =========================================================
# FINAL CHECK
# =========================================================

echo
echo "========================================================="
echo "HaplotypeCaller completed successfully"
echo "========================================================="
echo
echo "Sample:"
echo "$sampleID"
echo
echo "GVCF:"
echo "$gvcf"
echo
echo "Index:"
echo "${gvcf}.tbi"
echo
ls -lh \
    "$gvcf" \
    "${gvcf}.tbi"
echo
date

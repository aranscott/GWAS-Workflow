#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32g
#SBATCH --time=08:00:00
#SBATCH --output=markdup_v2_%A_%a.log
#SBATCH --error=markdup_v2_%A_%a.err
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=ARRAY_TASKS
#SBATCH --array=1-195%15

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR

module load biocontainers
module load gatk4/4.2.6.0

indir="${indir:-${PROJECT_DIR}/04.aligned_reads}"
outdir="${outdir:-${PROJECT_DIR}/05.marked_duplicates}"

mkdir -p "$outdir"

# Collect all aligned BAMs.
# This handles both naming conventions:
#   *_sort.bam
#   *_sorted.bam
mapfile -t bams < <(
    find "$indir" -maxdepth 1 -type f \
        \( -name '*_sorted.bam' -o -name '*_sort.bam' \) \
        | sort
)

echo "Total BAM files detected: ${#bams[@]}"

# Safety check:
# This run is expected to contain exactly 195 BAM files.
if [[ ${#bams[@]} -ne 195 ]]; then
    echo "ERROR: Expected 195 BAM files but found ${#bams[@]}."
    exit 1
fi

# Assign one BAM to this array task.
bam="${bams[$((SLURM_ARRAY_TASK_ID - 1))]}"

if [[ -z "${bam:-}" ]]; then
    echo "ERROR: No BAM found for task $SLURM_ARRAY_TASK_ID"
    echo "Total BAMs detected: ${#bams[@]}"
    exit 1
fi

# Derive sample ID from the BAM filename.
sampleID=$(basename "$bam")
sampleID="${sampleID%_sorted.bam}"
sampleID="${sampleID%_sort.bam}"

echo "========================================================="
echo "MarkDuplicates v2"
echo "Array task: $SLURM_ARRAY_TASK_ID"
echo "Sample: $sampleID"
echo "Input:  $bam"
echo "Output: $outdir"
echo "========================================================="
date

gatk MarkDuplicates \
    -I "$bam" \
    -O "${outdir}/${sampleID}_marked.bam" \
    -M "${outdir}/${sampleID}_metrics.txt" \
    --REMOVE_DUPLICATES false \
    --CREATE_INDEX true \
    --java-options "-Xmx28g"

echo
echo "GATK MarkDuplicates finished processing: $sampleID"
date

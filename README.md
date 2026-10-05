# GWAS-Workflow

Shell scripts used for GWAS and variant calling workflows.

## Scripts

- `bfctools_stats_wildtype.sh`
- `gatk_haplo_wild.sh`
- `gwas_top10_namplot_only.sh`
- `mark_duplicates_merged.sh`
- `run_fisher_manhat.sh`
- `run_fisher_scatter.sh`
- `run_gemma10.sh`
- `variant_caller_para_merged_v3.sh`
- `variant_calls.sh`
- `variant_filter_merged_v3.sh`
- `variant_filter_wildtype.sh`
- `wild_filt_max60.sh`

## Paths and inputs

The scripts derive project paths from `PROJECT_DIR`, which defaults to the current working directory. Submit them from the workflow/data root, or set `PROJECT_DIR` to that directory before submission. Inputs and outputs use the pipeline's numbered subdirectories.

The reference FASTA defaults to `reference/GCA_025370935.1_NRCan_CFum_1_genomic.fna` under `PROJECT_DIR`. Set `REFERENCE_DIR` or `REF_FASTA` to use another location. GEMMA association files, GFF annotation, and BAM lists also have defaults under `PROJECT_DIR`; their corresponding variables can be set in the scripts for a different layout.

The scripts require input datasets and cluster software/modules that are not included here. No sequencing data or generated results are included. Slurm notification email addresses are intentionally not embedded; add a local `#SBATCH --mail-user=...` directive if needed.

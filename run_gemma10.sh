#!/bin/bash
#SBATCH --job-name=gemma_top10_fisher
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=gemma_top10_fisher_%j.out
#SBATCH --error=gemma_top10_fisher_%j.err
#SBATCH --mail-type=ALL

set -euo pipefail

# Set PROJECT_DIR to the workflow/data root when submitting from elsewhere.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/reference}"
export PROJECT_DIR

module load r/4.5.0


# =========================================================
# PATHS
# =========================================================

OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/08.fisher}"

GEMMA_FILE="${GEMMA_FILE:-${PROJECT_DIR}/07.2.gemma_lmm_merged_v3/output/gwas_diapause_results.assoc.txt}"

FISHER_FILE=${OUT_DIR}/fisher_allelic_genomewide_fixed.tsv

GFF_FILE="${GFF_FILE:-${PROJECT_DIR}/GCF_025370935.1_NRCan_CFum_1_genomic_GCA_chrom_names.sorted.gff}"

R_SCRIPT=${OUT_DIR}/make_gemma_top10_vs_fisher.R


# ---------------------------------------------------------
# GEMMA intermediate files
# ---------------------------------------------------------

GEMMA_COMPACT=${OUT_DIR}/gemma_all_variants_compact.tsv
GEMMA_SORTED=${OUT_DIR}/gemma_all_variants_sorted_by_p.tsv


# ---------------------------------------------------------
# Top-10 GEMMA files
# ---------------------------------------------------------

TOP10=${OUT_DIR}/top10_gemma_loci.tsv

TOP10_BED=${OUT_DIR}/top10_gemma_loci.bed

TOP10_SORTED_BED=${OUT_DIR}/top10_gemma_loci.sorted.bed

TOP10_WINDOWS=${OUT_DIR}/top10_gemma_windows.bed


# ---------------------------------------------------------
# Fisher files
# ---------------------------------------------------------

FISHER_BED=${OUT_DIR}/fisher_allelic_variants.bed

FISHER_SORTED_BED=${OUT_DIR}/fisher_allelic_variants.sorted.bed

FISHER_INTERSECT=${OUT_DIR}/top10_gemma_fisher_intersections.tsv

TOP10_FISHER=${OUT_DIR}/top10_gemma_fisher_best.tsv


# ---------------------------------------------------------
# Gene annotation files
# ---------------------------------------------------------

ALL_GENES_BED=${OUT_DIR}/all_nuclear_genes.bed

ALL_GENES_SORTED=${OUT_DIR}/all_nuclear_genes.sorted.bed

TOP10_GENE=${OUT_DIR}/top10_gemma_nearest_gene.tsv


# ---------------------------------------------------------
# Chromosome order
# ---------------------------------------------------------

GENOME_ORDER=${OUT_DIR}/nuclear_chromosomes.genome


# =========================================================
# SETTINGS
# =========================================================

CLUSTER_DISTANCE=50000

FISHER_WINDOW=50000


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "GEMMA top-10 independent loci vs Fisher"
echo "========================================================="

echo
echo "GEMMA:"
echo "$GEMMA_FILE"

echo
echo "Fisher:"
echo "$FISHER_FILE"

echo
echo "GFF:"
echo "$GFF_FILE"

echo
echo "Independent-locus distance:"
echo "$CLUSTER_DISTANCE bp"

echo
echo "Fisher comparison window:"
echo "+/- $FISHER_WINDOW bp"

echo
echo "Output:"
echo "$OUT_DIR"

echo
date

echo "========================================================="


# =========================================================
# VERIFY INPUT FILES
# =========================================================

for f in \
    "$GEMMA_FILE" \
    "$FISHER_FILE" \
    "$GFF_FILE"
do

    if [[ ! -f "$f" ]]; then

        echo
        echo "ERROR: Missing required input:"
        echo "$f"

        exit 1

    fi

done


# =========================================================
# CREATE PROPER BEDTOOLS GENOME FILE
#
# This is generated directly from the GFF region records.
#
# Format:
#
# chromosome    chromosome_length
# =========================================================

echo
echo "Creating chromosome-order genome file..."

awk -F'\t' '
BEGIN {
    OFS="\t"
}

$0 !~ /^#/ &&
$3=="region" &&
$1 ~ /^CM046(102|103|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|119|120|121|122|123|124|125|126|127|128|129|130|131)\.1$/ {

    print $1,$5
}
' "$GFF_FILE" \
> "$GENOME_ORDER"


echo "Chromosomes in genome file:"
wc -l "$GENOME_ORDER"

cat "$GENOME_ORDER"


# =========================================================
# PREPARE GEMMA
#
# Current GEMMA structure:
#
# $1  chromosome
# $3  position
# $10 p-value
#
# Keep only the 30 nuclear chromosomes.
# =========================================================

echo
echo "Preparing compact GEMMA table..."

awk 'BEGIN {
    OFS="\t"
}

$1 ~ /^CM046(102|103|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|119|120|121|122|123|124|125|126|127|128|129|130|131)\.1$/ &&
$3 ~ /^[0-9]+$/ &&
$10 > 0 {

    print $1,$3,$10
}
' "$GEMMA_FILE" \
> "$GEMMA_COMPACT"


echo "Valid nuclear GEMMA variants:"
wc -l "$GEMMA_COMPACT"


# =========================================================
# SORT GEMMA BY P-VALUE
# =========================================================

echo
echo "Sorting GEMMA by P-value..."

LC_ALL=C sort \
    -t $'\t' \
    -k3,3g \
    "$GEMMA_COMPACT" \
> "$GEMMA_SORTED"


# =========================================================
# SELECT TOP 10 INDEPENDENT GEMMA LOCI
#
# Greedy selection:
#
# 1. Select lowest GEMMA P.
# 2. Reject subsequent variants <=50 kb from any
#    already-selected lead on the same chromosome.
# 3. Continue until 10 loci have been selected.
#
# IMPORTANT:
#
# This table remains in GEMMA significance order.
# It is NOT the BEDTools-sorted file.
# =========================================================

echo
echo "Selecting top 10 independent GEMMA loci..."

awk -v dist="$CLUSTER_DISTANCE" '

BEGIN {
    FS=OFS="\t"
}

{

    chr=$1
    pos=$2
    p=$3

    keep=1

    for (i=1; i<=n; i++) {

        if (chr == lead_chr[i]) {

            d=pos-lead_pos[i]

            if (d < 0) {
                d=-d
            }

            if (d <= dist) {

                keep=0
                break
            }
        }
    }


    if (keep && n < 10) {

        n++

        lead_chr[n]=chr
        lead_pos[n]=pos
        lead_p[n]=p

        print n,chr,pos,p
    }


    if (n >= 10) {
        exit
    }
}

' "$GEMMA_SORTED" \
> "${TOP10}.body"


{
    echo -e "RANK\tCHR\tGEMMA_POS\tGEMMA_P"
    cat "${TOP10}.body"
} > "$TOP10"


rm -f "${TOP10}.body"


echo
echo "Top 10 GEMMA loci:"
column -t -s $'\t' "$TOP10"


# =========================================================
# CREATE TOP-10 BED
#
# This BED initially remains in GEMMA rank order.
# =========================================================

awk 'BEGIN {
    OFS="\t"
}

NR>1 {

    print $2,
          $3-1,
          $3,
          "gemma_"$1,
          $4
}

' "$TOP10" \
> "$TOP10_BED"


# =========================================================
# SORT TOP-10 BED FOR BEDTOOLS
#
# THIS is the version used by bedtools closest.
#
# The rank-order BED remains untouched.
# =========================================================

echo
echo "Sorting top-10 GEMMA BED in genomic order..."

bedtools sort \
    -g "$GENOME_ORDER" \
    -i "$TOP10_BED" \
> "$TOP10_SORTED_BED"


echo "Coordinate-sorted top-10 GEMMA BED:"
cat "$TOP10_SORTED_BED"


# =========================================================
# CREATE +/-50 KB GEMMA WINDOWS
# =========================================================

awk -v flank="$FISHER_WINDOW" '
BEGIN {
    OFS="\t"
}

NR>1 {

    start=$3-flank-1

    if (start < 0) {
        start=0
    }

    end=$3+flank

    print $2,
          start,
          end,
          "gemma_"$1
}

' "$TOP10" \
> "${TOP10_WINDOWS}.unsorted"


bedtools sort \
    -g "$GENOME_ORDER" \
    -i "${TOP10_WINDOWS}.unsorted" \
> "$TOP10_WINDOWS"


rm -f "${TOP10_WINDOWS}.unsorted"


# =========================================================
# PREPARE FISHER BED
#
# Fisher columns:
#
# 1 CHR
# 2 SNP
# 3 BP
# 4 A1
# 5 C_A
# 6 C_U
# 7 A2
# 8 P
# 9 OR
# =========================================================

echo
echo "Preparing Fisher BED..."

awk 'BEGIN {
    OFS="\t"
}

NR>1 &&
$3 ~ /^[0-9]+$/ &&
$8 > 0 {

    print $1,
          $3-1,
          $3,
          $8
}

' "$FISHER_FILE" \
> "$FISHER_BED"


# =========================================================
# SORT FISHER BED
# =========================================================

echo
echo "Sorting Fisher BED..."

bedtools sort \
    -g "$GENOME_ORDER" \
    -i "$FISHER_BED" \
> "$FISHER_SORTED_BED"


# =========================================================
# FIND FISHER VARIANTS WITHIN +/-50 KB OF EACH GEMMA LEAD
# =========================================================

echo
echo "Finding Fisher variants within GEMMA windows..."

bedtools intersect \
    -a "$TOP10_WINDOWS" \
    -b "$FISHER_SORTED_BED" \
    -wa \
    -wb \
> "$FISHER_INTERSECT"


# =========================================================
# SELECT STRONGEST FISHER SNP PER GEMMA LOCUS
# =========================================================

echo
echo "Selecting strongest Fisher SNP per GEMMA locus..."

awk '
BEGIN {
    OFS="\t"
}

{

    locus=$4

    fisher_chr=$5
    fisher_pos=$7+0
    fisher_p=$8+0


    if (!(locus in best) || fisher_p < best[locus]) {

        best[locus]=fisher_p

        best_chr[locus]=fisher_chr

        best_pos[locus]=fisher_pos
    }
}


END {

    print "LOCUS",
          "FISHER_CHR",
          "FISHER_POS",
          "FISHER_P"


    for (locus in best) {

        print locus,
              best_chr[locus],
              best_pos[locus],
              best[locus]
    }
}

' "$FISHER_INTERSECT" \
> "$TOP10_FISHER"


# =========================================================
# BUILD ALL-NUCLEAR-GENES BED DIRECTLY FROM NEW GFF
#
# IMPORTANT:
#
# Keep ALL gene records.
#
# This includes:
#   named genes
#   LOC genes
#   uncharacterized genes
#
# That allows the summary table to retain the LOC ID even
# when the plot displays "Uncategorized".
#
# BED columns:
#
# 1 CHR
# 2 START
# 3 END
# 4 GENE_ID
# 5 NAME
# 6 DESCRIPTION
# =========================================================

echo
echo "Building nuclear gene BED directly from GFF..."

awk -F'\t' '
BEGIN {
    OFS="\t"
}

$0 !~ /^#/ &&
$3=="gene" &&
$1 ~ /^CM046(102|103|104|105|106|107|108|109|110|111|112|113|114|115|116|117|118|119|120|121|122|123|124|125|126|127|128|129|130|131)\.1$/ {

    gene_id=""
    name=""
    description=""

    n=split($9,a,";")


    for (i=1; i<=n; i++) {

        if (a[i] ~ /^ID=/) {

            sub(/^ID=/,"",a[i])
            sub(/^gene-/,"",a[i])

            gene_id=a[i]
        }


        else if (a[i] ~ /^Name=/) {

            sub(/^Name=/,"",a[i])

            name=a[i]
        }


        else if (a[i] ~ /^description=/) {

            sub(/^description=/,"",a[i])

            description=a[i]
        }
    }


    print $1,
          $4-1,
          $5,
          gene_id,
          name,
          description
}

' "$GFF_FILE" \
> "$ALL_GENES_BED"


# =========================================================
# VALIDATE GENE BED BEFORE SORTING
# =========================================================

echo
echo "Validating gene BED..."

if awk -F'\t' '
$2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {
    bad=1
    print "BAD GENE BED LINE:",NR,$0 > "/dev/stderr"
    exit
}
END {
    exit bad
}
' "$ALL_GENES_BED"
then

    echo "Gene BED validation passed."

else

    echo "ERROR: Gene BED validation failed."
    exit 1

fi


echo
echo "First five gene BED records:"

head -5 "$ALL_GENES_BED"


# =========================================================
# SORT GENE BED
# =========================================================

echo
echo "Sorting nuclear gene BED..."

bedtools sort \
    -g "$GENOME_ORDER" \
    -i "$ALL_GENES_BED" \
> "$ALL_GENES_SORTED"


# =========================================================
# FIND NEAREST GENE TO EACH GEMMA LEAD
#
# Both files are now explicitly coordinate-sorted.
# =========================================================

echo
echo "Finding nearest gene to each GEMMA lead..."

bedtools closest \
    -a "$TOP10_SORTED_BED" \
    -b "$ALL_GENES_SORTED" \
    -d \
> "$TOP10_GENE"


echo
echo "Nearest-gene results:"

column -t -s $'\t' "$TOP10_GENE"


# =========================================================
# WRITE R SCRIPT
# =========================================================

TOP10_FILE="$TOP10"
TOP10_FISHER_FILE="$TOP10_FISHER"
TOP10_GENE_FILE="$TOP10_GENE"
export GEMMA_FILE FISHER_FILE TOP10_FILE TOP10_FISHER_FILE TOP10_GENE_FILE OUT_DIR

cat > "$R_SCRIPT" << 'EOF'

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(grid)
})


# =========================================================
# FILES
# =========================================================

GEMMA_FILE <-
    Sys.getenv("GEMMA_FILE")

FISHER_FILE <-
    Sys.getenv("FISHER_FILE")

TOP10_FILE <-
    Sys.getenv("TOP10_FILE")

TOP10_FISHER_FILE <-
    Sys.getenv("TOP10_FISHER_FILE")

TOP10_GENE_FILE <-
    Sys.getenv("TOP10_GENE_FILE")

OUT_DIR <-
    Sys.getenv("OUT_DIR")


# =========================================================
# SETTINGS
# =========================================================

PLOT_CUTOFF <- 0.005

NUCLEAR_CHR <- paste0(
    "CM046",
    sprintf(
        "%03d",
        102:131
    ),
    ".1"
)


# =========================================================
# READ GEMMA
# =========================================================

cat("Reading GEMMA...\n")

gemma <- read.table(
    GEMMA_FILE,
    header=TRUE,
    sep="",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


gemma$chr <- as.character(
    gemma$chr
)

gemma$pos <- as.numeric(
    gemma$ps
)

gemma$p <- as.numeric(
    gemma$p_score
)


gemma <- gemma[
    gemma$chr %in% NUCLEAR_CHR &
    is.finite(gemma$pos) &
    is.finite(gemma$p) &
    gemma$p > 0,
]


gemma$chr_factor <- factor(
    gemma$chr,
    levels=NUCLEAR_CHR,
    ordered=TRUE
)

gemma$logp <- -log10(
    gemma$p
)


cat(
    "Valid nuclear GEMMA variants:",
    nrow(gemma),
    "\n"
)


# =========================================================
# READ FISHER
# =========================================================

cat("Reading Fisher...\n")

fisher <- read.table(
    FISHER_FILE,
    header=TRUE,
    sep="",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


fisher$chr <- as.character(
    fisher$CHR
)

fisher$pos <- as.numeric(
    fisher$BP
)

fisher$p <- as.numeric(
    fisher$P
)


fisher <- fisher[
    fisher$chr %in% NUCLEAR_CHR &
    is.finite(fisher$pos) &
    is.finite(fisher$p) &
    fisher$p > 0,
]


fisher$chr_factor <- factor(
    fisher$chr,
    levels=NUCLEAR_CHR,
    ordered=TRUE
)

fisher$logp <- -log10(
    fisher$p
)


cat(
    "Valid nuclear Fisher variants:",
    nrow(fisher),
    "\n"
)


# =========================================================
# COMMON GENOME COORDINATES
# =========================================================

coord <- rbind(

    data.frame(
        chr=as.character(
            gemma$chr_factor
        ),
        pos=gemma$pos
    ),

    data.frame(
        chr=as.character(
            fisher$chr_factor
        ),
        pos=fisher$pos
    )
)


coord$chr <- factor(
    coord$chr,
    levels=NUCLEAR_CHR,
    ordered=TRUE
)


chr_lengths <- tapply(
    coord$pos,
    coord$chr,
    max,
    na.rm=TRUE
)


chr_lengths <- chr_lengths[
    NUCLEAR_CHR
]


offsets <- c(
    0,
    cumsum(
        chr_lengths
    )
)[1:length(chr_lengths)]


names(offsets) <- NUCLEAR_CHR


chr_centers <- data.frame(
    chr=NUCLEAR_CHR,
    center=offsets + chr_lengths/2
)


# =========================================================
# CUMULATIVE POSITIONS
# =========================================================

gemma$cum_pos <-
    gemma$pos +
    offsets[
        as.character(
            gemma$chr_factor
        )
    ]


fisher$cum_pos <-
    fisher$pos +
    offsets[
        as.character(
            fisher$chr_factor
        )
    ]


# =========================================================
# ALTERNATING CHROMOSOME COLORS
# =========================================================

gemma$chr_group <- factor(
    as.numeric(
        gemma$chr_factor
    ) %% 2,
    levels=c(0,1)
)


fisher$chr_group <- factor(
    as.numeric(
        fisher$chr_factor
    ) %% 2,
    levels=c(0,1)
)


# =========================================================
# DISPLAY FILTER
# =========================================================

gemma_plot <- gemma[
    gemma$p <= PLOT_CUTOFF,
]


fisher_plot <- fisher[
    fisher$p <= PLOT_CUTOFF,
]


# =========================================================
# READ TOP 10
# =========================================================

top10 <- read.table(
    TOP10_FILE,
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


top10$CHR <- as.character(
    top10$CHR
)


top10$GEMMA_POS <- as.numeric(
    top10$GEMMA_POS
)


top10$GEMMA_P <- as.numeric(
    top10$GEMMA_P
)


top10$gemma_logp <- -log10(
    top10$GEMMA_P
)


top10$gemma_cum <-
    top10$GEMMA_POS +
    offsets[
        top10$CHR
    ]


# =========================================================
# READ FISHER MATCHES
# =========================================================

top10_fisher <- read.table(
    TOP10_FISHER_FILE,
    header=TRUE,
    sep="\t",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


top10_fisher$rank <- as.numeric(
    sub(
        "^gemma_",
        "",
        top10_fisher$LOCUS
    )
)


top10$FISHER_POS <- NA_real_
top10$FISHER_P <- NA_real_
top10$fisher_logp <- NA_real_
top10$fisher_cum <- NA_real_


for (i in seq_len(nrow(top10))) {

    idx <- which(
        top10_fisher$rank ==
        top10$RANK[i]
    )


    if (length(idx) == 0) {
        next
    }


    idx <- idx[1]


    top10$FISHER_POS[i] <-
        as.numeric(
            top10_fisher$FISHER_POS[idx]
        )


    top10$FISHER_P[i] <-
        as.numeric(
            top10_fisher$FISHER_P[idx]
        )


    top10$fisher_logp[i] <-
        -log10(
            top10$FISHER_P[i]
        )


    top10$fisher_cum[i] <-
        top10$FISHER_POS[i] +
        offsets[
            top10$CHR[i]
        ]
}


# =========================================================
# READ NEAREST-GENE RESULTS
#
# Columns:
#
# 1 GEMMA chr
# 2 GEMMA start
# 3 GEMMA end
# 4 GEMMA ID
# 5 GEMMA P
#
# 6 gene chr
# 7 gene start
# 8 gene end
# 9 gene ID
# 10 gene Name
# 11 description
# 12 distance
# =========================================================

gene_ann <- read.table(
    TOP10_GENE_FILE,
    header=FALSE,
    sep="\t",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


names(gene_ann) <- c(
    "gemma_chr",
    "gemma_start",
    "gemma_end",
    "gemma_id",
    "gemma_p",
    "gene_chr",
    "gene_start",
    "gene_end",
    "gene_id",
    "gene_name",
    "gene_description",
    "gene_distance"
)


gene_ann$rank <- as.numeric(
    sub(
        "^gemma_",
        "",
        gene_ann$gemma_id
    )
)


# =========================================================
# LABEL RULES
#
# Figure:
#
# actual biological gene name
#        ↓
# informative description
#        ↓
# Uncategorized
#
# LOC identifier is preserved in the summary but never
# displayed on the figure.
# =========================================================

good_label <- function(x) {

    if (
        is.na(x) ||
        x == ""
    ) {
        return(FALSE)
    }


    if (
        grepl(
            "^LOC[0-9]+$",
            x,
            ignore.case=TRUE
        )
    ) {
        return(FALSE)
    }


    if (
        grepl(
            "^uncharacterized",
            x,
            ignore.case=TRUE
        )
    ) {
        return(FALSE)
    }


    TRUE
}


get_plot_label <- function(
    name,
    description
) {

    if (
        good_label(name)
    ) {

        return(name)
    }


    if (
        good_label(description)
    ) {

        return(description)
    }


    return(
        "Uncategorized"
    )
}


# =========================================================
# ADD GENE INFORMATION TO TOP 10
# =========================================================

top10$gene_id <- ""
top10$gene_name <- ""
top10$gene_description <- ""
top10$gene_distance <- NA_real_
top10$plot_label <- "Uncategorized"


for (i in seq_len(nrow(top10))) {

    idx <- which(
        gene_ann$rank ==
        top10$RANK[i]
    )


    if (length(idx) == 0) {
        next
    }


    idx <- idx[1]


    top10$gene_id[i] <-
        as.character(
            gene_ann$gene_id[idx]
        )


    top10$gene_name[i] <-
        as.character(
            gene_ann$gene_name[idx]
        )


    top10$gene_description[i] <-
        as.character(
            gene_ann$gene_description[idx]
        )


    top10$gene_distance[i] <-
        as.numeric(
            gene_ann$gene_distance[idx]
        )


    top10$plot_label[i] <-
        get_plot_label(
            top10$gene_name[i],
            top10$gene_description[i]
        )
}


# =========================================================
# SAVE FINAL SUMMARY
# =========================================================

summary_file <- file.path(
    OUT_DIR,
    "gemma_top10_vs_fisher_summary.tsv"
)


write.table(
    top10[
        ,
        c(
            "RANK",
            "CHR",
            "GEMMA_POS",
            "GEMMA_P",
            "FISHER_POS",
            "FISHER_P",
            "gene_id",
            "gene_name",
            "gene_description",
            "gene_distance",
            "plot_label"
        )
    ],
    file=summary_file,
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)


cat("\n")
cat("=========================================================\n")
cat("FINAL TOP-10 SUMMARY\n")
cat("=========================================================\n")


print(
    top10[
        ,
        c(
            "RANK",
            "CHR",
            "GEMMA_POS",
            "GEMMA_P",
            "FISHER_P",
            "gene_id",
            "gene_distance",
            "plot_label"
        )
    ],
    row.names=FALSE
)


# =========================================================
# PLOT RANGE
# =========================================================

y_max <- max(
    c(
        gemma_plot$logp,
        fisher_plot$logp,
        top10$gemma_logp,
        top10$fisher_logp
    ),
    na.rm=TRUE
)


y_upper <- y_max * 1.10


x_limits <- c(
    min(offsets),
    max(
        offsets +
        chr_lengths
    )
)


# =========================================================
# GEMMA PANEL
# =========================================================

gemma_panel <- ggplot(
    gemma_plot,
    aes(
        x=cum_pos,
        y=logp
    )
) +

    geom_point(
        aes(
            color=chr_group
        ),
        size=0.55,
        alpha=0.60
    ) +

    scale_color_manual(
        values=c(
            "0"="black",
            "1"="grey55"
        ),
        guide="none"
    ) +

    geom_point(
        data=top10,
        inherit.aes=FALSE,
        aes(
            x=gemma_cum,
            y=gemma_logp
        ),
        shape=21,
        fill="red3",
        color="black",
        stroke=0.8,
        size=4.8
    ) +

    scale_x_continuous(
        limits=x_limits,
        breaks=chr_centers$center,
        labels=NULL,
        expand=c(0,0)
    ) +

    scale_y_continuous(
        limits=c(
            0,
            y_upper
        ),
        expand=c(0,0)
    ) +

    labs(
        x=NULL,
        y=expression(
            -log[10](P)
        ),
        title="GEMMA linear mixed model"
    ) +

    theme_classic() +

    theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),

        plot.title=element_text(
            size=15,
            face="bold"
        ),

        axis.title.y=element_text(
            size=11
        ),

        plot.margin=margin(
            t=10,
            r=20,
            b=5,
            l=20
        )
    )


# =========================================================
# FISHER PANEL
# =========================================================

fisher_panel <- ggplot(
    fisher_plot,
    aes(
        x=cum_pos,
        y=logp
    )
) +

    geom_point(
        aes(
            color=chr_group
        ),
        size=0.55,
        alpha=0.60
    ) +

    scale_color_manual(
        values=c(
            "0"="black",
            "1"="grey55"
        ),
        guide="none"
    ) +

    geom_point(
        data=top10,
        inherit.aes=FALSE,
        aes(
            x=fisher_cum,
            y=fisher_logp
        ),
        shape=21,
        fill="red3",
        color="black",
        stroke=0.8,
        size=4.8
    ) +

    scale_x_continuous(
        limits=x_limits,
        breaks=chr_centers$center,
        labels=NUCLEAR_CHR,
        expand=c(0,0)
    ) +

    scale_y_continuous(
        limits=c(
            0,
            y_upper
        ),
        expand=c(0,0)
    ) +

    labs(
        x="Chromosome",
        y=expression(
            -log[10](P)
        ),
        title="Allelic Fisher exact test"
    ) +

    theme_classic() +

    theme(
        axis.text.x=element_text(
            angle=45,
            hjust=1,
            size=7
        ),

        plot.title=element_text(
            size=15,
            face="bold"
        ),

        axis.title=element_text(
            size=11
        ),

        plot.margin=margin(
            t=10,
            r=20,
            b=15,
            l=20
        )
    )


# =========================================================
# LABEL STRIP
#
# Every GEMMA top-10 locus gets a label.
#
# Unknown annotation = Uncategorized.
# =========================================================

labels <- top10


labels$label_y <- seq(
    1,
    7,
    length.out=nrow(labels)
)


label_panel <- ggplot(
    labels,
    aes(
        x=gemma_cum,
        y=label_y
    )
) +

    geom_segment(
        aes(
            x=gemma_cum,
            xend=gemma_cum,
            y=0.2,
            yend=label_y-0.25
        ),
        color="grey40",
        linewidth=0.45
    ) +

    geom_text_repel(
        aes(
            label=plot_label
        ),
        direction="y",
        force=4,
        force_pull=0.8,
        box.padding=0.8,
        point.padding=0.2,
        min.segment.length=0,
        segment.color="grey40",
        segment.size=0.45,
        size=3.5,
        fontface="bold",
        max.overlaps=Inf,
        seed=12345
    ) +

    scale_x_continuous(
        limits=x_limits,
        expand=c(0,0)
    ) +

    coord_cartesian(
        ylim=c(
            0,
            8
        ),
        clip="off"
    ) +

    theme_void() +

    theme(
        plot.margin=margin(
            t=5,
            r=20,
            b=5,
            l=20
        )
    )


# =========================================================
# ALIGN PANELS
# =========================================================

gemma_grob <- ggplotGrob(
    gemma_panel
)


fisher_grob <- ggplotGrob(
    fisher_panel
)


label_grob <- ggplotGrob(
    label_panel
)


max_width <- grid::unit.pmax(
    gemma_grob$widths,
    fisher_grob$widths,
    label_grob$widths
)


gemma_grob$widths <- max_width

fisher_grob$widths <- max_width

label_grob$widths <- max_width


# =========================================================
# DRAW FIGURE
# =========================================================

draw_combined <- function() {

    grid.newpage()


    pushViewport(
        viewport(
            layout=grid.layout(
                nrow=3,
                ncol=1,

                heights=unit(
                    c(
                        0.26,
                        0.37,
                        0.37
                    ),
                    "npc"
                )
            )
        )
    )


    pushViewport(
        viewport(
            layout.pos.row=1
        )
    )

    grid.draw(
        label_grob
    )

    popViewport()


    pushViewport(
        viewport(
            layout.pos.row=2
        )
    )

    grid.draw(
        gemma_grob
    )

    popViewport()


    pushViewport(
        viewport(
            layout.pos.row=3
        )
    )

    grid.draw(
        fisher_grob
    )

    popViewport()


    popViewport()
}


# =========================================================
# SAVE PNG
# =========================================================

png_file <- file.path(
    OUT_DIR,
    "gemma_top10_vs_fisher_aligned_final.png"
)


png(
    filename=png_file,
    width=5600,
    height=3900,
    res=300
)


draw_combined()


dev.off()


# =========================================================
# SAVE PDF
# =========================================================

pdf_file <- file.path(
    OUT_DIR,
    "gemma_top10_vs_fisher_aligned_final.pdf"
)


pdf(
    file=pdf_file,
    width=18,
    height=12
)


draw_combined()


dev.off()


# =========================================================
# COMPLETE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("GEMMA top-10 vs Fisher completed successfully\n")
cat("=========================================================\n")


cat(
    "Top GEMMA loci:",
    nrow(top10),
    "\n"
)


cat(
    "Named/informatively annotated:",
    sum(
        top10$plot_label !=
        "Uncategorized"
    ),
    "\n"
)


cat(
    "Uncategorized:",
    sum(
        top10$plot_label ==
        "Uncategorized"
    ),
    "\n"
)


cat("\nPNG:\n")
cat(
    png_file,
    "\n"
)


cat("\nPDF:\n")
cat(
    pdf_file,
    "\n"
)


cat("\nSummary table:\n")
cat(
    summary_file,
    "\n"
)


cat("=========================================================\n")

EOF


# =========================================================
# RUN R
# =========================================================

echo
echo "Running R plotting script..."

Rscript "$R_SCRIPT"


echo
echo "========================================================="
echo "GEMMA top-10 vs Fisher job completed."
echo "========================================================="

date

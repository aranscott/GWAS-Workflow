#!/bin/bash

#SBATCH --job-name=gwas_top10_manhattan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=gwas_top10_manhattan_%j.out
#SBATCH --error=gwas_top10_manhattan_%j.err
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

GFF_FILE="${GFF_FILE:-${PROJECT_DIR}/GCF_025370935.1_NRCan_CFum_1_genomic_GCA_chrom_names.sorted.gff}"

R_SCRIPT=${OUT_DIR}/make_gwas_top10_manhattan.R


# =========================================================
# SETTINGS
# =========================================================

CLUSTER_DISTANCE=50000

PLOT_CUTOFF=0.005


# =========================================================
# START
# =========================================================

echo "========================================================="
echo "GWAS diapause vs. non-diapause"
echo "Top 10 independent GEMMA loci"
echo "========================================================="
echo
echo "GEMMA:"
echo "$GEMMA_FILE"
echo
echo "GFF:"
echo "$GFF_FILE"
echo
echo "Independent-locus distance:"
echo "$CLUSTER_DISTANCE bp"
echo
echo "Plot cutoff:"
echo "P <= $PLOT_CUTOFF"
echo
echo "Output:"
echo "$OUT_DIR"
echo
date
echo "========================================================="


# =========================================================
# VERIFY INPUTS
# =========================================================

for f in \
    "$GEMMA_FILE" \
    "$GFF_FILE"
do

    if [[ ! -f "$f" ]]; then
        echo
        echo "ERROR: Missing input:"
        echo "$f"
        exit 1
    fi

done


# =========================================================
# WRITE R SCRIPT
# =========================================================

export GEMMA_FILE GFF_FILE OUT_DIR

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

GFF_FILE <-
    Sys.getenv("GFF_FILE")

OUT_DIR <-
    Sys.getenv("OUT_DIR")


# =========================================================
# SETTINGS
# =========================================================

CLUSTER_DISTANCE <- 50000

# =========================================================
# CHROMOSOME NAMES
#
# CM046102.1 = CHR 2
# ...
# CM046130.1 = CHR 30
# CM046131.1 = Z
# =========================================================

REF_CHR <- paste0(
    "CM046",
    sprintf(
        "%03d",
        102:131
    ),
    ".1"
)


DISPLAY_CHR <- c(
    paste(
        "CHR",
        2:30
    ),
    "Z"
)


chr_display_map <- setNames(
    DISPLAY_CHR,
    REF_CHR
)


# =========================================================
# READ GEMMA
# =========================================================

cat("Reading GEMMA results...\n")

gemma <- read.table(
    GEMMA_FILE,
    header=TRUE,
    sep="",
    stringsAsFactors=FALSE,
    check.names=FALSE,
    quote="",
    comment.char=""
)


required <- c(
    "chr",
    "ps",
    "p_score"
)


missing <- setdiff(
    required,
    names(gemma)
)


if (length(missing) > 0) {

    stop(
        paste(
            "Missing GEMMA columns:",
            paste(
                missing,
                collapse=", "
            )
        )
    )
}


gemma$chr <- as.character(
    gemma$chr
)


gemma$pos <- as.numeric(
    gemma$ps
)


gemma$p <- as.numeric(
    gemma$p_score
)


# ---------------------------------------------------------
# Keep only the 30 chromosomes of interest.
# This includes Z.
# ---------------------------------------------------------

gemma <- gemma[
    gemma$chr %in% REF_CHR &
    is.finite(gemma$pos) &
    is.finite(gemma$p) &
    gemma$p > 0,
]


gemma$chr_factor <- factor(
    gemma$chr,
    levels=REF_CHR,
    ordered=TRUE
)


gemma$logp <- -log10(
    gemma$p
)


cat(
    "Valid GEMMA variants:",
    nrow(gemma),
    "\n"
)


# =========================================================
# CALCULATE CHROMOSOME LENGTHS
# =========================================================

chr_lengths <- tapply(
    gemma$pos,
    gemma$chr_factor,
    max,
    na.rm=TRUE
)


chr_lengths <- chr_lengths[
    REF_CHR
]


if (
    any(
        !is.finite(chr_lengths)
    )
) {

    stop(
        "One or more chromosomes have no valid GEMMA positions."
    )
}


# =========================================================
# CUMULATIVE GENOME POSITION
# =========================================================

offsets <- c(
    0,
    cumsum(
        chr_lengths
    )
)[1:length(chr_lengths)]


names(offsets) <- REF_CHR


gemma$cum_pos <-
    gemma$pos +
    offsets[
        as.character(
            gemma$chr_factor
        )
    ]


# =========================================================
# CHROMOSOME CENTERS
# =========================================================

chr_centers <- data.frame(
    chr=REF_CHR,
    center=offsets +
        chr_lengths / 2,
    display=DISPLAY_CHR,
    stringsAsFactors=FALSE
)


# =========================================================
# CHROMOSOME GROUPS
# =========================================================

gemma$chr_group <- factor(
    as.numeric(
        gemma$chr_factor
    ) %% 2,
    levels=c(0,1)
)


# =========================================================
# SELECT TOP 10 INDEPENDENT GEMMA LOCI
#
# Start with the strongest GEMMA SNP.
#
# Reject additional SNPs within 50 kb on the same
# chromosome of an already-selected lead.
# =========================================================

cat(
    "Selecting top 10 independent GEMMA loci...\n"
)


gemma_sorted <- gemma[
    order(
        gemma$p
    ),
]


top10_rows <- list()


for (
    i in seq_len(
        nrow(gemma_sorted)
    )
) {

    if (
        length(top10_rows) >= 10
    ) {
        break
    }


    candidate <- gemma_sorted[i,]


    keep <- TRUE


    if (
        length(top10_rows) > 0
    ) {

        for (
            j in seq_along(
                top10_rows
            )
        ) {

            selected <- top10_rows[[j]]


            if (
                candidate$chr ==
                selected$chr
            ) {

                distance <- abs(
                    candidate$pos -
                    selected$pos
                )


                if (
                    distance <=
                    CLUSTER_DISTANCE
                ) {

                    keep <- FALSE

                    break
                }
            }
        }
    }


    if (keep) {

        top10_rows[[length(top10_rows) + 1]] <- candidate
    }
}


if (
    length(top10_rows) < 10
) {

    stop(
        paste0(
            "Only ",
            length(top10_rows),
            " independent GEMMA loci could be selected."
        )
    )
}


top10 <- do.call(
    rbind,
    top10_rows
)


rownames(top10) <- NULL


top10$rank <- seq_len(
    nrow(top10)
)


# =========================================================
# TOP-10 CUMULATIVE POSITIONS
# =========================================================

top10$cum_pos <-
    top10$pos +
    offsets[
        top10$chr
    ]


top10$display_chr <-
    chr_display_map[
        top10$chr
    ]


top10$logp <- -log10(
    top10$p
)


# =========================================================
# EXTRACT GENE RECORDS FROM NEW GFF
#
# Keep ALL genes so that we can preserve LOC identifiers
# even when the visible figure label becomes "Uncategorized".
# =========================================================

cat(
    "Reading GFF gene records...\n"
)


gff_lines <- readLines(
    GFF_FILE
)


gff_lines <- gff_lines[
    !grepl(
        "^#",
        gff_lines
    )
]


gene_lines <- gff_lines[
    vapply(
        strsplit(
            gff_lines,
            "\t",
            fixed=TRUE
        ),
        function(x) {

            length(x) == 9 &&
            x[3] == "gene"

        },
        logical(1)
    )
]


# =========================================================
# GFF ATTRIBUTE PARSER
# =========================================================

parse_attributes <- function(
    attr_string
) {

    pieces <- strsplit(
        attr_string,
        ";",
        fixed=TRUE
    )[[1]]


    result <- list()


    for (
        piece in pieces
    ) {

        kv <- strsplit(
            piece,
            "=",
            fixed=TRUE
        )[[1]]


        if (
            length(kv) >= 2
        ) {

            result[[kv[1]]] <-
                paste(
                    kv[-1],
                    collapse="="
                )
        }
    }


    result
}


# =========================================================
# BUILD GENE TABLE
# =========================================================

genes <- data.frame(
    CHR=character(),
    START=numeric(),
    END=numeric(),
    GENE_ID=character(),
    NAME=character(),
    DESCRIPTION=character(),
    stringsAsFactors=FALSE
)


for (
    line in gene_lines
) {

    fields <- strsplit(
        line,
        "\t",
        fixed=TRUE
    )[[1]]


    attrs <- parse_attributes(
        fields[9]
    )


    gene_id <- ""

    if (
        !is.null(
            attrs[["ID"]]
        )
    ) {

        gene_id <- sub(
            "^gene-",
            "",
            attrs[["ID"]]
        )
    }


    name <- ""

    if (
        !is.null(
            attrs[["Name"]]
        )
    ) {

        name <- attrs[["Name"]]
    }


    description <- ""

    if (
        !is.null(
            attrs[["description"]]
        )
    ) {

        description <-
            attrs[["description"]]
    }


    genes <- rbind(
        genes,
        data.frame(
            CHR=fields[1],
            START=as.numeric(fields[4]),
            END=as.numeric(fields[5]),
            GENE_ID=gene_id,
            NAME=name,
            DESCRIPTION=description,
            stringsAsFactors=FALSE
        )
    )
}


cat(
    "Gene records loaded:",
    nrow(genes),
    "\n"
)


# =========================================================
# LABEL FUNCTION
#
# Priority:
#
# 1. Actual Name
# 2. Informative description
# 3. Uncategorized
#
# LOC, XM, XP, and uncharacterized entries are not displayed.
# =========================================================

usable_name <- function(
    x
) {

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
            "^XM_[0-9.]+$",
            x
        )
    ) {

        return(FALSE)
    }


    if (
        grepl(
            "^XP_[0-9.]+$",
            x
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


make_label <- function(
    name,
    description
) {

    if (
        usable_name(name)
    ) {

        return(name)
    }


    if (
        usable_name(description)
    ) {

        return(description)
    }


    "Uncategorized"
}


# =========================================================
# FIND NEAREST GENE FOR EACH TOP-10 GEMMA LEAD
# =========================================================

top10$gene_id <- ""

top10$gene_name <- ""

top10$gene_description <- ""

top10$gene_distance <- NA_real_

top10$plot_label <- "Uncategorized"


for (
    i in seq_len(
        nrow(top10)
    )
) {

    candidates <- genes[
        genes$CHR == top10$chr[i],
    ]


    if (
        nrow(candidates) == 0
    ) {

        next
    }


    # -----------------------------------------------------
    # First look for a gene containing the lead SNP.
    # -----------------------------------------------------

    containing <- candidates[
        candidates$START <= top10$pos[i] &
        candidates$END >= top10$pos[i],
    ]


    if (
        nrow(containing) > 0
    ) {

        chosen <- containing[1,]


        top10$gene_distance[i] <- 0


    } else {

        # -------------------------------------------------
        # Otherwise choose nearest gene.
        # -------------------------------------------------

        distance <- ifelse(
            top10$pos[i] <
            candidates$START,

            candidates$START -
                top10$pos[i],

            ifelse(
                top10$pos[i] >
                candidates$END,

                top10$pos[i] -
                    candidates$END,

                0
            )
        )


        chosen_idx <- which.min(
            distance
        )


        chosen <- candidates[
            chosen_idx,
        ]


        top10$gene_distance[i] <-
            distance[
                chosen_idx
            ]
    }


    top10$gene_id[i] <-
        chosen$GENE_ID


    top10$gene_name[i] <-
        chosen$NAME


    top10$gene_description[i] <-
        chosen$DESCRIPTION


    top10$plot_label[i] <-
        make_label(
            chosen$NAME,
            chosen$DESCRIPTION
        )
}


# =========================================================
# SAVE TOP-10 SUMMARY
# =========================================================

summary_file <- file.path(
    OUT_DIR,
    "gwas_top10_loci_summary.tsv"
)


write.table(
    top10[
        ,
        c(
            "rank",
            "chr",
            "display_chr",
            "pos",
            "p",
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


# =========================================================
# PRINT TOP 10
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("TOP 10 GEMMA LOCI\n")
cat("=========================================================\n")


print(
    top10[
        ,
        c(
            "rank",
            "display_chr",
            "pos",
            "p",
            "gene_id",
            "gene_distance",
            "plot_label"
        )
    ],
    row.names=FALSE
)


# =========================================================
# PLOT DATA
# =========================================================

plot_data <- gemma


# =========================================================
# Y LIMIT
# =========================================================

y_max <- max(
    c(
        plot_data$logp,
        top10$logp
    ),
    na.rm=TRUE
)


y_upper <- y_max * 1.18


# =========================================================
# MANHATTAN PLOT
# =========================================================

p <- ggplot(
    plot_data,
    aes(
        x=cum_pos,
        y=logp
    )
) +

    # -----------------------------------------------------
    # Background points
    # -----------------------------------------------------

    geom_point(
        aes(
            color=chr_group
        ),
        size=0.35,
        alpha=0.45
    ) +

    scale_color_manual(
        values=c(
            "0"="black",
            "1"="grey55"
        ),
        guide="none"
    ) +

    # -----------------------------------------------------
    # Highlight ten strongest independent GEMMA loci
    # -----------------------------------------------------

    geom_point(
        data=top10,
        inherit.aes=FALSE,
        aes(
            x=cum_pos,
            y=logp
        ),
        shape=21,
        fill="red3",
        color="black",
        stroke=0.8,
        size=5
    ) +

    # -----------------------------------------------------
    # Labels
    #
    # ggrepel draws the connector itself. This avoids the
    # manually drawn connector passing through a label.
    # -----------------------------------------------------

    geom_text_repel(
        data=top10,
        inherit.aes=FALSE,
        aes(
            x=cum_pos,
            y=logp,
            label=plot_label
        ),

        direction="y",

        nudge_y=0.7,

        force=5,

        force_pull=0.8,

        box.padding=0.9,

        point.padding=0.35,

        min.segment.length=0,

        segment.color="grey40",

        segment.size=0.45,

        size=3.5,

        fontface="bold",

        max.overlaps=Inf,

        seed=12345
    ) +

    # -----------------------------------------------------
    # Genome-wide axis
    # -----------------------------------------------------

    scale_x_continuous(
        limits=c(
            min(offsets),
            max(
                offsets +
                chr_lengths
            )
        ),

        breaks=chr_centers$center,

        labels=chr_centers$display,

        expand=c(0,0)
    ) +

scale_y_continuous(
    limits=c(
        0,
        y_upper
    ),
    breaks=seq(
        0,
        ceiling(y_upper / 2) * 2,
        by=2
    ),
    expand=c(0,0)
) +

    labs(
        x="Chromosome",

        y=expression(
            -log[10](P)
        ),

        title="GWAS diapause vs. non-diapause"
    ) +

    theme_classic() +

    theme(

        plot.title=element_text(
            size=18,
            face="bold"
        ),

        axis.title.x=element_text(
            size=12
        ),

        axis.title.y=element_text(
            size=12
        ),

        axis.text.x=element_text(
            size=8,
            angle=45,
            hjust=1
        ),

        axis.text.y=element_text(
            size=9
        ),

        plot.margin=margin(
            t=20,
            r=25,
            b=20,
            l=20
        )
    )


# =========================================================
# SAVE PNG
# =========================================================

png_file <- file.path(
    OUT_DIR,
    "gwas_diapause_vs_nondiapause_top10.png"
)


ggsave(
    filename=png_file,
    plot=p,
    width=18,
    height=8.5,
    units="in",
    dpi=300
)


# =========================================================
# SAVE PDF
# =========================================================

pdf_file <- file.path(
    OUT_DIR,
    "gwas_diapause_vs_nondiapause_top10.pdf"
)


ggsave(
    filename=pdf_file,
    plot=p,
    width=18,
    height=8.5,
    units="in"
)


# =========================================================
# COMPLETE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("GWAS Manhattan plot completed successfully\n")
cat("=========================================================\n")

cat(
    "Top 10 independent loci:",
    nrow(top10),
    "\n"
)

cat(
    "PNG:",
    png_file,
    "\n"
)

cat(
    "PDF:",
    pdf_file,
    "\n"
)

cat(
    "Summary:",
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
echo "GWAS top-10 Manhattan plot completed."
echo "========================================================="

date

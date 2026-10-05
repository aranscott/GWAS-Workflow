#!/bin/bash
#SBATCH --job-name=fisher_gemma_compare
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=fisher_gemma_compare_%j.out
#SBATCH --error=fisher_gemma_compare_%j.err
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

FISHER_DIR="${FISHER_DIR:-${PROJECT_DIR}/08.fisher}"
CONV_DIR=${FISHER_DIR}/fisher_gemma_convergence

CONVERGENCE_FILE=${CONV_DIR}/fisher_gemma_convergence_final.tsv

TIER1_FILE=${CONV_DIR}/fisher_three_model_high_priority_fixed.tsv
TIER2_FILE=${CONV_DIR}/fisher_three_model_tier2_intergenic.tsv
TIER3_FILE=${CONV_DIR}/fisher_three_model_tier3_nearby.tsv

GFF_FILE="${GFF_FILE:-${PROJECT_DIR}/GCF_025370935.1_NRCan_CFum_1_genomic_GCA_chrom_names.sorted.gff}"

OUT_DIR=${FISHER_DIR}

R_SCRIPT=${OUT_DIR}/make_fisher_gemma_comparison.R

echo "========================================================="
echo "Fisher vs GEMMA comparison plot"
echo "========================================================="
echo
echo "Convergence:"
echo "$CONVERGENCE_FILE"
echo
echo "Tier 1:"
echo "$TIER1_FILE"
echo
echo "Tier 2:"
echo "$TIER2_FILE"
echo
echo "Tier 3:"
echo "$TIER3_FILE"
echo
echo "GFF:"
echo "$GFF_FILE"
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
    "$CONVERGENCE_FILE" \
    "$TIER1_FILE" \
    "$TIER2_FILE" \
    "$TIER3_FILE" \
    "$GFF_FILE"
do

    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing required file:"
        echo "$f"
        exit 1
    fi

done


# =========================================================
# WRITE R SCRIPT
# =========================================================

export CONVERGENCE_FILE TIER1_FILE TIER2_FILE TIER3_FILE GFF_FILE OUT_DIR

cat > "$R_SCRIPT" << 'EOF'

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
})


# =========================================================
# SETTINGS
# =========================================================

CONVERGENCE_FILE <-
    Sys.getenv("CONVERGENCE_FILE")

TIER1_FILE <-
    Sys.getenv("TIER1_FILE")

TIER2_FILE <-
    Sys.getenv("TIER2_FILE")

TIER3_FILE <-
    Sys.getenv("TIER3_FILE")

GFF_FILE <-
    Sys.getenv("GFF_FILE")

OUT_DIR <-
    Sys.getenv("OUT_DIR")


# =========================================================
# READ CONVERGENCE TABLE
# =========================================================

cat("Reading Fisher/GEMMA convergence table...\n")

dat <- read.table(
    CONVERGENCE_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

cat(
    "Rows read:",
    nrow(dat),
    "\n"
)


# =========================================================
# CHECK REQUIRED COLUMNS
# =========================================================

required_columns <- c(
    "LOCUS",
    "CHR",
    "BEST_FISHER_P",
    "BEST_GEMMA_P",
    "FISHER_LEAD_POS",
    "GEMMA_LEAD_POS"
)

missing_columns <- setdiff(
    required_columns,
    names(dat)
)

if (length(missing_columns) > 0) {

    stop(
        paste(
            "Missing required columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


# =========================================================
# CLEAN STATISTICS
# =========================================================

dat$chr <- as.character(
    dat$CHR
)

dat$fisher_pos <- as.numeric(
    dat$FISHER_LEAD_POS
)

dat$gemma_pos <- as.numeric(
    dat$GEMMA_LEAD_POS
)

dat$fisher_p <- as.numeric(
    dat$BEST_FISHER_P
)

dat$gemma_p <- as.numeric(
    dat$BEST_GEMMA_P
)

dat$fisher_logp <- -log10(
    dat$fisher_p
)

dat$gemma_logp <- -log10(
    dat$gemma_p
)


# =========================================================
# REMOVE INVALID ROWS
# =========================================================

dat <- dat[
    is.finite(dat$fisher_pos) &
    is.finite(dat$gemma_pos) &
    is.finite(dat$fisher_p) &
    is.finite(dat$gemma_p) &
    dat$fisher_p > 0 &
    dat$gemma_p > 0,
]


# =========================================================
# READ TIERS
# =========================================================

read_tier <- function(
    file,
    tier_name
) {

    x <- read.table(
        file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    x$tier <- tier_name

    x
}


tier1 <- read_tier(
    TIER1_FILE,
    "Tier 1"
)

tier2 <- read_tier(
    TIER2_FILE,
    "Tier 2"
)

tier3 <- read_tier(
    TIER3_FILE,
    "Tier 3"
)


# =========================================================
# BUILD TIER LOOKUP
# =========================================================

tier_map <- rbind(
    tier1[, c("LOCUS", "tier")],
    tier2[, c("LOCUS", "tier")],
    tier3[, c("LOCUS", "tier")]
)

tier_map <- tier_map[
    !duplicated(
        tier_map$LOCUS
    ),
]


tier_idx <- match(
    dat$LOCUS,
    tier_map$LOCUS
)

dat$tier <- tier_map$tier[
    tier_idx
]

dat$tier[
    is.na(dat$tier)
] <- "Other"


# =========================================================
# READ THE NEW GFF
#
# We parse it line-by-line rather than using read.table()
# on the entire GFF.
# =========================================================

cat("Reading GFF gene records...\n")

gff_lines <- readLines(
    GFF_FILE
)

gff_lines <- gff_lines[
    !grepl(
        "^#",
        gff_lines
    )
]


# =========================================================
# EXTRACT ONLY GENE LINES
# =========================================================

gene_lines <- gff_lines[
    vapply(
        strsplit(
            gff_lines,
            "\t",
            fixed = TRUE
        ),
        function(x) {

            length(x) == 9 &&
            x[3] == "gene"

        },
        logical(1)
    )
]


cat(
    "Gene records found:",
    length(gene_lines),
    "\n"
)


# =========================================================
# ATTRIBUTE PARSER
# =========================================================

parse_attributes <- function(
    attr_string
) {

    pieces <- strsplit(
        attr_string,
        ";",
        fixed = TRUE
    )[[1]]

    result <- list()

    for (piece in pieces) {

        kv <- strsplit(
            piece,
            "=",
            fixed = TRUE
        )[[1]]

        if (length(kv) >= 2) {

            key <- kv[1]

            value <- paste(
                kv[-1],
                collapse = "="
            )

            result[[key]] <- value
        }
    }

    return(result)
}


# =========================================================
# BUILD GENE ANNOTATION DATA FRAME
# =========================================================

genes <- data.frame(
    CHR = character(),
    START = numeric(),
    END = numeric(),
    GENE_ID = character(),
    NAME = character(),
    DESCRIPTION = character(),
    LABEL = character(),
    stringsAsFactors = FALSE
)


for (line in gene_lines) {

    fields <- strsplit(
        line,
        "\t",
        fixed = TRUE
    )[[1]]

    attrs <- parse_attributes(
        fields[9]
    )


    # -----------------------------------------------------
    # Extract attributes
    # -----------------------------------------------------

    gene_id <- ""

    if (!is.null(attrs[["ID"]])) {

        gene_id <- sub(
            "^gene-",
            "",
            attrs[["ID"]]
        )
    }


    name <- ""

    if (!is.null(attrs[["Name"]])) {

        name <- attrs[["Name"]]
    }


    description <- ""

    if (!is.null(attrs[["description"]])) {

        description <- attrs[["description"]]
    }


    # -----------------------------------------------------
    # Determine useful biological label
    # -----------------------------------------------------

    label <- ""


    # First preference:
    # a genuine gene name/symbol.
    #
    # Examples:
    # per
    # tim
    # Tret1
    #
    # Do not use:
    # LOC123456
    # XM_...
    # XP_...
    # uncharacterized
    #

    usable_name <-
        name != "" &&
        !grepl(
            "^LOC[0-9]+$",
            name,
            ignore.case = TRUE
        ) &&
        !grepl(
            "^XM_[0-9.]+$",
            name
        ) &&
        !grepl(
            "^XP_[0-9.]+$",
            name
        ) &&
        !grepl(
            "^uncharacterized",
            name,
            ignore.case = TRUE
        )


    if (usable_name) {

        label <- name
    }


    # -----------------------------------------------------
    # Second preference:
    # informative description.
    # -----------------------------------------------------

    if (
        label == "" &&
        description != "" &&
        !grepl(
            "^uncharacterized",
            description,
            ignore.case = TRUE
        ) &&
        !grepl(
            "^LOC[0-9]+$",
            description,
            ignore.case = TRUE
        )
    ) {

        label <- description
    }


    # -----------------------------------------------------
    # Add gene record.
    # -----------------------------------------------------

    genes <- rbind(
        genes,
        data.frame(
            CHR = fields[1],
            START = as.numeric(
                fields[4]
            ),
            END = as.numeric(
                fields[5]
            ),
            GENE_ID = gene_id,
            NAME = name,
            DESCRIPTION = description,
            LABEL = label,
            stringsAsFactors = FALSE
        )
    )
}


cat(
    "Gene records loaded:",
    nrow(genes),
    "\n"
)


# =========================================================
# FIND THE GENE CONTAINING EACH FISHER LEAD SNP
# =========================================================

dat$gene_label <- ""


for (i in seq_len(nrow(dat))) {

    candidates <- genes[
        genes$CHR == dat$chr[i] &
        genes$START <= dat$fisher_pos[i] &
        genes$END >= dat$fisher_pos[i] &
        genes$LABEL != "",
    ]


    if (nrow(candidates) > 0) {

        dat$gene_label[i] <-
            candidates$LABEL[1]
    }
}


# =========================================================
# SELECT THE FIVE STRONGEST NAMED LOCI
#
# This is deliberately based on Fisher P regardless of
# Tier 1/2/3.
# =========================================================

named <- dat[
    !is.na(dat$gene_label) &
    dat$gene_label != "",
]


named <- named[
    order(
        named$fisher_p
    ),
]


if (nrow(named) < 5) {

    stop(
        paste0(
            "Only ",
            nrow(named),
            " loci have usable gene names. ",
            "Five cannot be selected."
        )
    )
}


top5 <- named[
    seq_len(5),
]


# =========================================================
# SAVE TOP-FIVE TABLE
# =========================================================

write.table(
    top5[
        ,
        c(
            "LOCUS",
            "CHR",
            "fisher_pos",
            "fisher_p",
            "gemma_pos",
            "gemma_p",
            "tier",
            "gene_label"
        )
    ],
    file = file.path(
        OUT_DIR,
        "fisher_gemma_top5_named_comparison.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# =========================================================
# PRINT TOP FIVE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("TOP FIVE NAMED FISHER/GEMMA LOCI\n")
cat("=========================================================\n")

print(
    top5[
        ,
        c(
            "LOCUS",
            "CHR",
            "fisher_pos",
            "fisher_p",
            "gemma_pos",
            "gemma_p",
            "tier",
            "gene_label"
        )
    ],
    row.names = FALSE
)


# =========================================================
# TIER COLORS
# =========================================================

tier_colors <- c(
    "Tier 1" = "red3",
    "Tier 2" = "darkorange2",
    "Tier 3" = "dodgerblue3",
    "Other" = "grey65"
)


# =========================================================
# MAIN SCATTER PLOT
# =========================================================

max_value <- max(
    c(
        dat$fisher_logp,
        dat$gemma_logp
    ),
    na.rm = TRUE
)


axis_max <- max_value * 1.08


p <- ggplot(
    dat,
    aes(
        x = fisher_logp,
        y = gemma_logp
    )
) +

    # -----------------------------------------------------
    # 1:1 reference
    # -----------------------------------------------------

    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.60,
        color = "grey45"
    ) +

    # -----------------------------------------------------
    # All loci
    # -----------------------------------------------------

    geom_point(
        aes(
            color = tier
        ),
        size = 2.0,
        alpha = 0.65
    ) +

    # -----------------------------------------------------
    # Top five labels
    # -----------------------------------------------------

    geom_text_repel(
        data = top5,
        aes(
            label = gene_label,
            color = tier
        ),
        size = 3.5,
        fontface = "bold",
        box.padding = 0.75,
        point.padding = 0.35,
        force = 2.5,
        force_pull = 0.5,
        min.segment.length = 0,
        segment.color = "grey40",
        segment.size = 0.45,
        max.overlaps = Inf,
        seed = 12345
    ) +

    scale_color_manual(
        values = tier_colors,
        name = "Convergence tier"
    ) +

    scale_x_continuous(
        limits = c(
            0,
            axis_max
        ),
        expand = c(0,0)
    ) +

    scale_y_continuous(
        limits = c(
            0,
            axis_max
        ),
        expand = c(0,0)
    ) +

    labs(
        x = expression(
            -log[10]("Fisher P")
        ),

        y = expression(
            -log[10]("GEMMA P")
        ),

        title =
            "Fisher vs GEMMA association strength",

        subtitle =
            "703 Fisher-defined loci; strongest five named loci labeled"
    ) +

    theme_classic() +

    theme(
        plot.title = element_text(
            size = 17,
            face = "bold"
        ),

        plot.subtitle = element_text(
            size = 11
        ),

        axis.title = element_text(
            size = 12
        ),

        axis.text = element_text(
            size = 9
        ),

        legend.title = element_text(
            size = 10
        ),

        legend.text = element_text(
            size = 9
        ),

        plot.margin = margin(
            t = 15,
            r = 20,
            b = 15,
            l = 20
        )
    )


# =========================================================
# SAVE PNG
# =========================================================

png_file <- file.path(
    OUT_DIR,
    "fisher_vs_gemma_locus_comparison.png"
)

ggsave(
    filename = png_file,
    plot = p,
    width = 10,
    height = 10,
    units = "in",
    dpi = 300
)


# =========================================================
# SAVE PDF
# =========================================================

pdf_file <- file.path(
    OUT_DIR,
    "fisher_vs_gemma_locus_comparison.pdf"
)

ggsave(
    filename = pdf_file,
    plot = p,
    width = 10,
    height = 10,
    units = "in"
)


# =========================================================
# COMPLETE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("Fisher vs GEMMA comparison complete\n")
cat("=========================================================\n")

cat(
    "Total Fisher-defined loci:",
    nrow(dat),
    "\n"
)

cat(
    "Tier 1:",
    sum(dat$tier == "Tier 1"),
    "\n"
)

cat(
    "Tier 2:",
    sum(dat$tier == "Tier 2"),
    "\n"
)

cat(
    "Tier 3:",
    sum(dat$tier == "Tier 3"),
    "\n"
)

cat(
    "Other:",
    sum(dat$tier == "Other"),
    "\n"
)

cat("\nPNG:\n")
cat(png_file, "\n")

cat("\nPDF:\n")
cat(pdf_file, "\n")

cat("\nTop-five table:\n")
cat(
    file.path(
        OUT_DIR,
        "fisher_gemma_top5_named_comparison.tsv"
    ),
    "\n"
)

cat("=========================================================\n")

EOF


# =========================================================
# RUN R
# =========================================================

echo
echo "Running R comparison plot..."

Rscript "$R_SCRIPT"

echo
echo "========================================================="
echo "Fisher vs GEMMA comparison job completed."
echo "========================================================="
date

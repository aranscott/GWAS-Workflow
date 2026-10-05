#!/bin/bash
#SBATCH --job-name=fisher_gemma_top5
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=fisher_gemma_top5_%j.out
#SBATCH --error=fisher_gemma_top5_%j.err
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

FISHER_FILE=${FISHER_DIR}/fisher_allelic_genomewide.tsv

GEMMA_FILE="${GEMMA_FILE:-${PROJECT_DIR}/07.2.gemma_lmm_merged_v3/output/gwas_diapause_results.assoc.txt}"

TIER1_FILE=${CONV_DIR}/fisher_three_model_high_priority_fixed.tsv
TIER2_FILE=${CONV_DIR}/fisher_three_model_tier2_intergenic.tsv
TIER3_FILE=${CONV_DIR}/fisher_three_model_tier3_nearby.tsv

GFF_FILE="${GFF_FILE:-${PROJECT_DIR}/GCF_025370935.1_NRCan_CFum_1_genomic_GCA_chrom_names.sorted.gff}"

OUT_DIR=${FISHER_DIR}

R_SCRIPT=${OUT_DIR}/make_fisher_gemma_top5_manhattan_v3.R

GFF_GENE_TABLE=${OUT_DIR}/top5_gene_annotation_source.tsv

echo "========================================================="
echo "Fisher + GEMMA top-5 aligned Manhattan plot"
echo "========================================================="
echo
echo "Fisher:"
echo "$FISHER_FILE"
echo
echo "GEMMA:"
echo "$GEMMA_FILE"
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
    "$FISHER_FILE" \
    "$GEMMA_FILE" \
    "$TIER1_FILE" \
    "$TIER2_FILE" \
    "$TIER3_FILE" \
    "$GFF_FILE"
do

    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing input file:"
        echo "$f"
        exit 1
    fi

done


# =========================================================
# EXTRACT GENE FEATURES FROM THE GFF
#
# This avoids having R parse the entire GFF.
#
# Output columns:
#
# CHR
# START
# END
# STRAND
# GENE_ID
# NAME
# DESCRIPTION
#
# Only actual "gene" records are retained.
# =========================================================

echo
echo "Extracting gene annotations from GFF..."

awk -F'\t' '
BEGIN {
    OFS="\t"
    print "CHR","START","END","STRAND",
          "GENE_ID","NAME","DESCRIPTION"
}

$0 !~ /^#/ && $3=="gene" {

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

    print $1,$4,$5,$7,gene_id,name,description
}
' "$GFF_FILE" \
> "$GFF_GENE_TABLE"

echo "Gene annotation records extracted:"
wc -l "$GFF_GENE_TABLE"


# =========================================================
# WRITE R SCRIPT
# =========================================================

export FISHER_FILE GEMMA_FILE TIER1_FILE TIER2_FILE TIER3_FILE GFF_GENE_TABLE OUT_DIR

cat > "$R_SCRIPT" << 'EOF'

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(grid)
})


# =========================================================
# SETTINGS
# =========================================================

FISHER_FILE <-
    Sys.getenv("FISHER_FILE")

GEMMA_FILE <-
    Sys.getenv("GEMMA_FILE")

TIER1_FILE <-
    Sys.getenv("TIER1_FILE")

TIER2_FILE <-
    Sys.getenv("TIER2_FILE")

TIER3_FILE <-
    Sys.getenv("TIER3_FILE")

GFF_GENE_TABLE <-
    Sys.getenv("GFF_GENE_TABLE")

OUT_DIR <-
    Sys.getenv("OUT_DIR")

PLOT_CUTOFF <- 0.005

TOP_N <- 5

NUCLEAR_CHR <- paste0(
    "CM046",
    sprintf("%03d", 102:131),
    ".1"
)


# =========================================================
# READ FISHER
# =========================================================

cat("Reading Fisher results...\n")

fisher <- read.table(
    FISHER_FILE,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

fisher$chr <- as.character(fisher$CHR)
fisher$pos <- as.numeric(fisher$BP)
fisher$p <- as.numeric(fisher$P)

fisher <- fisher[
    fisher$chr %in% NUCLEAR_CHR &
    is.finite(fisher$pos) &
    is.finite(fisher$p) &
    fisher$p > 0,
]

fisher$chr <- factor(
    fisher$chr,
    levels = NUCLEAR_CHR,
    ordered = TRUE
)

fisher$logp <- -log10(fisher$p)

cat(
    "Valid nuclear Fisher variants:",
    nrow(fisher),
    "\n"
)


# =========================================================
# READ GEMMA
# =========================================================

cat("Reading GEMMA results...\n")

gemma <- read.table(
    GEMMA_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

gemma$chr <- as.character(gemma$chr)
gemma$pos <- as.numeric(gemma$ps)
gemma$p <- as.numeric(gemma$p_score)

gemma <- gemma[
    gemma$chr %in% NUCLEAR_CHR &
    is.finite(gemma$pos) &
    is.finite(gemma$p) &
    gemma$p > 0,
]

gemma$chr <- factor(
    gemma$chr,
    levels = NUCLEAR_CHR,
    ordered = TRUE
)

gemma$logp <- -log10(gemma$p)

cat(
    "Valid nuclear GEMMA variants:",
    nrow(gemma),
    "\n"
)


# =========================================================
# COMMON CHROMOSOME COORDINATES
# =========================================================

coord_data <- rbind(
    fisher[, c("chr","pos")],
    gemma[, c("chr","pos")]
)

coord_data$chr <- factor(
    coord_data$chr,
    levels = NUCLEAR_CHR,
    ordered = TRUE
)

chr_lengths <- tapply(
    coord_data$pos,
    coord_data$chr,
    max,
    na.rm = TRUE
)

chr_lengths <- chr_lengths[NUCLEAR_CHR]

if (any(!is.finite(chr_lengths))) {
    stop(
        "One or more nuclear chromosomes lack valid coordinates."
    )
}

offsets <- c(
    0,
    cumsum(chr_lengths)
)[1:length(chr_lengths)]

names(offsets) <- NUCLEAR_CHR

chr_centers <- data.frame(
    chr = NUCLEAR_CHR,
    center = offsets + chr_lengths / 2
)


# =========================================================
# CUMULATIVE POSITIONS
# =========================================================

fisher$cum_pos <-
    fisher$pos +
    offsets[as.character(fisher$chr)]

gemma$cum_pos <-
    gemma$pos +
    offsets[as.character(gemma$chr)]


# =========================================================
# CHROMOSOME GROUPS
# =========================================================

fisher$chr_group <- factor(
    as.numeric(fisher$chr) %% 2,
    levels = c(0,1)
)

gemma$chr_group <- factor(
    as.numeric(gemma$chr) %% 2,
    levels = c(0,1)
)


# =========================================================
# PLOT CUTOFF
# =========================================================

fisher_plot <- fisher[
    fisher$p <= PLOT_CUTOFF,
]

gemma_plot <- gemma[
    gemma$p <= PLOT_CUTOFF,
]


# =========================================================
# READ TIER TABLES
# =========================================================

read_tier <- function(file, tier_name) {

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
# COMBINE TIERS
# =========================================================

all_tiers <- rbind(
    tier1,
    tier2,
    tier3
)

# One row per locus.
all_tiers <- all_tiers[
    !duplicated(all_tiers$LOCUS),
]


# =========================================================
# READ GFF GENE TABLE
# =========================================================

cat(
    "Reading extracted GFF gene annotations...\n"
)

genes <- read.table(
    GFF_GENE_TABLE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

genes$START <- as.numeric(
    genes$START
)

genes$END <- as.numeric(
    genes$END
)

cat(
    "GFF genes loaded:",
    nrow(genes),
    "\n"
)


# =========================================================
# DETERMINE A BIOLOGICAL LABEL
# =========================================================
#
# Priority:
#
# 1. Name= when it is a real gene name
# 2. description= when it is biologically informative
# 3. Otherwise NA
#
# LOCxxxxx, XM_, XP_, and "uncharacterized" are NOT
# accepted as final labels.
# =========================================================

is_bad_name <- function(x) {

    if (is.na(x) || x == "") {
        return(TRUE)
    }

    grepl(
        "^LOC[0-9]+$",
        x,
        ignore.case = TRUE
    ) ||
    grepl(
        "^XM_[0-9.]+$",
        x
    ) ||
    grepl(
        "^XP_[0-9.]+$",
        x
    ) ||
    grepl(
        "^uncharacterized",
        x,
        ignore.case = TRUE
    )
}


genes$LABEL <- NA_character_


for (i in seq_len(nrow(genes))) {

    nm <- genes$NAME[i]
    desc <- genes$DESCRIPTION[i]

    if (!is_bad_name(nm)) {

        genes$LABEL[i] <- nm

    } else if (!is_bad_name(desc)) {

        genes$LABEL[i] <- desc
    }
}


# =========================================================
# FIND A NAMED GENE CONTAINING EACH FISHER LEAD
#
# We work directly from chromosome + position.
# No intermediate list of data frames is needed.
# =========================================================

all_tiers$gene_label <- NA_character_
all_tiers$gene_start <- NA_real_
all_tiers$gene_end <- NA_real_
all_tiers$gene_id <- NA_character_


for (i in seq_len(nrow(all_tiers))) {

    current_chr <- as.character(
        all_tiers$CHR[i]
    )

    current_pos <- as.numeric(
        all_tiers$FISHER_LEAD_POS[i]
    )

    candidates <- genes[
        genes$CHR == current_chr &
        genes$START <= current_pos &
        genes$END >= current_pos &
        !is.na(genes$LABEL) &
        genes$LABEL != "",
    ]

    if (nrow(candidates) > 0) {

        # Use the first matching annotated gene.
        all_tiers$gene_label[i] <-
            candidates$LABEL[1]

        all_tiers$gene_start[i] <-
            candidates$START[1]

        all_tiers$gene_end[i] <-
            candidates$END[1]

        all_tiers$gene_id[i] <-
            candidates$GENE_ID[1]
    }
}


# =========================================================
# RETAIN ONLY LOCI WITH A USEFUL BIOLOGICAL NAME
# =========================================================

eligible <- all_tiers[
    !is.na(all_tiers$gene_label) &
    all_tiers$gene_label != "",
]


cat(
    "Loci with usable biological annotation:",
    nrow(eligible),
    "\n"
)


# =========================================================
# SELECT FIVE STRONGEST ELIGIBLE LOCI
#
# IMPORTANT:
# We exclude loci with no usable biological annotation.
# Then select the five strongest by Fisher P regardless
# of their Tier 1/2/3 status.
# =========================================================

eligible <- all_tiers[
    !is.na(all_tiers$gene_label) &
    all_tiers$gene_label != "",
]


eligible <- eligible[
    order(
        as.numeric(
            eligible$BEST_FISHER_P
        )
    ),
]


if (nrow(eligible) < TOP_N) {

    stop(
        paste0(
            "Only ",
            nrow(eligible),
            " loci have usable gene names; ",
            "cannot select five."
        )
    )
}


top5 <- eligible[
    seq_len(TOP_N),
]


# =========================================================
# WRITE TOP-FIVE TABLE
# =========================================================

write.table(
    top5[
        ,
        c(
            "LOCUS",
            "CHR",
            "FISHER_START",
            "FISHER_END",
            "FISHER_LEAD_POS",
            "BEST_FISHER_P",
            "GEMMA_LEAD_POS",
            "BEST_GEMMA_P",
            "LEAD_DISTANCE_BP",
            "tier",
            "gene_label",
            "gene_id"
        )
    ],
    file = file.path(
        OUT_DIR,
        "top5_fisher_gemma_loci.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# =========================================================
# REPORT TOP FIVE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("TOP FIVE SELECTED LOCI\n")
cat("=========================================================\n")

print(
    top5[
        ,
        c(
            "LOCUS",
            "CHR",
            "FISHER_LEAD_POS",
            "BEST_FISHER_P",
            "BEST_GEMMA_P",
            "LEAD_DISTANCE_BP",
            "tier",
            "gene_label"
        )
    ],
    row.names = FALSE
)


# =========================================================
# PREPARE TOP-FIVE COORDINATES
# =========================================================

top5$fisher_cum <-
    as.numeric(
        top5$FISHER_LEAD_POS
    ) +
    offsets[
        as.character(top5$CHR)
    ]

top5$gemma_cum <-
    as.numeric(
        top5$GEMMA_LEAD_POS
    ) +
    offsets[
        as.character(top5$CHR)
    ]

top5$fisher_logp <-
    -log10(
        as.numeric(
            top5$BEST_FISHER_P
        )
    )

top5$gemma_logp <-
    -log10(
        as.numeric(
            top5$BEST_GEMMA_P
        )
    )


# =========================================================
# TIER COLORS
# =========================================================

tier_colors <- c(
    "Tier 1" = "red3",
    "Tier 2" = "darkorange2",
    "Tier 3" = "dodgerblue3"
)


# =========================================================
# COMMON Y LIMIT
# =========================================================

y_max <- max(
    c(
        fisher_plot$logp,
        gemma_plot$logp,
        top5$fisher_logp,
        top5$gemma_logp
    ),
    na.rm = TRUE
)

y_upper <- y_max * 1.12


# =========================================================
# COMMON X LIMIT
# =========================================================

x_limits <- c(
    min(offsets),
    max(
        offsets + chr_lengths
    )
)


# =========================================================
# FISHER PANEL
# =========================================================

fisher_panel <- ggplot(
    fisher_plot,
    aes(
        x = cum_pos,
        y = logp
    )
) +

    geom_point(
        aes(
            color = chr_group
        ),
        size = 0.55,
        alpha = 0.60
    ) +

    scale_color_manual(
        values = c(
            "0" = "black",
            "1" = "grey55"
        ),
        guide = "none"
    ) +

    geom_point(
        data = top5,
        inherit.aes = FALSE,
        aes(
            x = fisher_cum,
            y = fisher_logp,
            fill = tier
        ),
        shape = 21,
        color = "black",
        stroke = 0.80,
        size = 4.8
    ) +

    scale_fill_manual(
        values = tier_colors,
        name = "Convergence tier"
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = chr_centers$center,
        labels = NULL,
        expand = c(0,0)
    ) +

    scale_y_continuous(
        limits = c(
            0,
            y_upper
        ),
        expand = c(0,0)
    ) +

    labs(
        x = NULL,
        y = expression(-log[10](P)),
        title = "Allelic Fisher exact test"
    ) +

    theme_classic() +

    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),

        plot.title = element_text(
            size = 15,
            face = "bold"
        ),

        plot.margin = margin(
            t = 10,
            r = 20,
            b = 5,
            l = 20
        )
    )


# =========================================================
# GEMMA PANEL
# =========================================================

gemma_panel <- ggplot(
    gemma_plot,
    aes(
        x = cum_pos,
        y = logp
    )
) +

    geom_point(
        aes(
            color = chr_group
        ),
        size = 0.55,
        alpha = 0.60
    ) +

    scale_color_manual(
        values = c(
            "0" = "black",
            "1" = "grey55"
        ),
        guide = "none"
    ) +

    geom_point(
        data = top5,
        inherit.aes = FALSE,
        aes(
            x = gemma_cum,
            y = gemma_logp,
            fill = tier
        ),
        shape = 21,
        color = "black",
        stroke = 0.80,
        size = 4.8
    ) +

    scale_fill_manual(
        values = tier_colors,
        name = "Convergence tier"
    ) +

    scale_x_continuous(
        limits = x_limits,
        breaks = chr_centers$center,
        labels = NUCLEAR_CHR,
        expand = c(0,0)
    ) +

    scale_y_continuous(
        limits = c(
            0,
            y_upper
        ),
        expand = c(0,0)
    ) +

    labs(
        x = "Chromosome",
        y = expression(-log[10](P)),
        title = "GEMMA linear mixed model"
    ) +

    theme_classic() +

    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 7
        ),

        plot.title = element_text(
            size = 15,
            face = "bold"
        ),

        plot.margin = margin(
            t = 10,
            r = 20,
            b = 15,
            l = 20
        ),

        legend.position = "none"
    )


# =========================================================
# TOP-FIVE LABEL STRIP
# =========================================================

# ---------------------------------------------------------
# Start labels at deliberately different vertical levels.
# ggrepel can further separate them.
# ---------------------------------------------------------

labels <- top5

labels$initial_y <- seq(
    1,
    5,
    length.out = nrow(labels)
)


label_panel <- ggplot(
    labels,
    aes(
        x = fisher_cum,
        y = initial_y
    )
) +

    geom_segment(
        aes(
            x = fisher_cum,
            xend = fisher_cum,
            y = 0.15,
            yend = initial_y - 0.20
        ),
        color = "grey40",
        linewidth = 0.45
    ) +

    geom_text_repel(
        aes(
            label = gene_label,
            color = tier
        ),
        direction = "y",

        force = 3,

        force_pull = 0.7,

        box.padding = 0.75,

        point.padding = 0.20,

        min.segment.length = 0,

        segment.color = "grey40",

        segment.size = 0.45,

        size = 3.5,

        fontface = "bold",

        max.overlaps = Inf,

        seed = 12345
    ) +

    scale_color_manual(
        values = tier_colors,
        guide = "none"
    ) +

    scale_x_continuous(
        limits = x_limits,
        expand = c(0,0)
    ) +

    coord_cartesian(
        ylim = c(
            0,
            6
        ),
        clip = "off"
    ) +

    theme_void() +

    theme(
        plot.margin = margin(
            t = 5,
            r = 20,
            b = 5,
            l = 20
        )
    )


# =========================================================
# ALIGN PANELS
# =========================================================

fisher_grob <- ggplotGrob(
    fisher_panel
)

gemma_grob <- ggplotGrob(
    gemma_panel
)

label_grob <- ggplotGrob(
    label_panel
)

max_width <- grid::unit.pmax(
    fisher_grob$widths,
    gemma_grob$widths,
    label_grob$widths
)

fisher_grob$widths <- max_width
gemma_grob$widths <- max_width
label_grob$widths <- max_width


# =========================================================
# DRAW COMBINED FIGURE
# =========================================================

draw_combined <- function() {

    grid.newpage()

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 3,
                ncol = 1,
                heights = unit(
                    c(
                        0.25,
                        0.375,
                        0.375
                    ),
                    "npc"
                )
            )
        )
    )

    pushViewport(
        viewport(
            layout.pos.row = 1
        )
    )

    grid.draw(
        label_grob
    )

    popViewport()


    pushViewport(
        viewport(
            layout.pos.row = 2
        )
    )

    grid.draw(
        fisher_grob
    )

    popViewport()


    pushViewport(
        viewport(
            layout.pos.row = 3
        )
    )

    grid.draw(
        gemma_grob
    )

    popViewport()

    popViewport()
}


# =========================================================
# SAVE PNG
# =========================================================

png_file <- file.path(
    OUT_DIR,
    "aligned_fisher_gemma_top5_named_v3.png"
)

png(
    filename = png_file,
    width = 5600,
    height = 3900,
    res = 300
)

draw_combined()

dev.off()


# =========================================================
# SAVE PDF
# =========================================================

pdf_file <- file.path(
    OUT_DIR,
    "aligned_fisher_gemma_top5_named_v3.pdf"
)

pdf(
    file = pdf_file,
    width = 18,
    height = 12
)

draw_combined()

dev.off()


# =========================================================
# COMPLETE
# =========================================================

cat("\n")
cat("=========================================================\n")
cat("Fisher + GEMMA top-five figure completed successfully\n")
cat("=========================================================\n")

cat(
    "Top-five loci:",
    nrow(top5),
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
    "Top-five table:",
    file.path(
        OUT_DIR,
        "top5_fisher_gemma_loci.tsv"
    ),
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
echo "Job completed successfully."
echo "========================================================="
date

#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2020-01-07 20:53:17
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse")

    if (!require("optparse", quietly = TRUE)) {
        stop("optparse not installed, please install it first.")
    }
}

if (!require("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")

    if (!require("ggplot2", quietly = TRUE)) {
        stop("ggplot2 not installed, please install it first.")
    }
}

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
library(reshape2)

opt_list <- list(
    make_option(
        c("-i", "--indir"),
        type = "character",
        help = "Directory of input(RNAmetrcs output).",
        metavar = "character"
        ),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        help = "Directory of output.",
        metavar = "character"
        ),
    make_option(
        c("--project_name"),
        type = "character",
        help = "Prefix of output files.",
        metavar = "character"
        ),
    make_option(
        c("--suffix"),
        type = "character",
        default = "geneBodyCoverage.txt",
        help = "Suffix of geneBody_coverage.py output, [default \"%default\"].",
        metavar = "character"
        )
)

opts <- parse_args(OptionParser(option_list = opt_list, 
                                usage = "usage: %prog [options]", 
                                add_help_option = TRUE, 
                                prog = "", 
                                description = "Summary CollectRnaMetrics output files and plot it.")
                   )

if (is.null(opts$indir)) {
    stop("-i --indir not set")
} else if (is.null(opts$outdir)) {
    stop("-o --outdir not set")
} else if (is.null(opts$project_name)) {
    stop("-o --project_name not set")
}

merge_files <- function(files) {
    filelist <- lapply(files, function(f) {
        filename <- sub(paste0("[^0-9A-Za-z]*", opts$suffix, "$"), "", basename(f))
        tmp <- read.table(f, header = TRUE, sep = "\t", row.names = 1)
        rownames(tmp) <- filename
        return(tmp)
    })
    read_cover <- do.call(rbind, filelist)
    return(read_cover)
}

files <- list.files(opts$indir, pattern = opts$suffix, full.names = TRUE)

read_cover <- merge_files(files)
read_cover <- read_cover/rowMeans(read_cover)

# ---------------------------------------------------------------------------- #
read_cover$Run <- rownames(read_cover)
read_cover <- melt(read_cover, id.vars = c('Run'))
colnames(read_cover)[2] <- "Pos"
read_cover$Pos <- as.numeric(as.vector(sub("^X", "", read_cover$Pos)))

width <- floor(length(files)/20) + 8

pdf(paste0(opts$project_name, "_RSeQC_geneBodyCoverage.pdf"), width = width, height = 6)
ggplot(read_cover, aes(x = Pos, y = value, color = Run)) + 
    geom_line() + 
    theme_bw() + 
    theme(axis.title = element_text(size = 8,face = "bold"), 
          axis.text = element_text(size = 6,face = "bold"), 
          plot.title = element_text(hjust = 0.5, size = 12,face = "bold"), 
          axis.text.x=element_text(angle=45,hjust=1, vjust=1),
          legend.text = element_text(size=6, face="bold")) + 
    labs(title = "Gene coverage(RSeQC)", x = "Gene body percentile (5'->3')", y = "")
dev.off()

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
        default = "RNAmetrics.txt",
        help = "Suffix of CollectRnaMetrics output, [default \"%default\"].",
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
        tmp <- readLines(f)
        tmp[8] <- paste0(tmp[8], "\t")
        read_dist <- data.frame(do.call(cbind, strsplit(tmp[7:8], "\t")))
        read_dist[, 2] <- as.numeric(as.vector(read_dist[, 2]))
        colnames(read_dist) <- c("Feature", filename)
        read_cover <- data.frame(do.call(rbind, strsplit(tmp[12:112], "\t")))
        read_cover <- data.frame(lapply(read_cover, function(x) {as.numeric(as.vector(x))}))
        colnames(read_cover) <- c("Pos", filename)
        return(list(read_dist = read_dist, read_cover = read_cover))
    })
    read_dist <- lapply(1:length(filelist), function(i) {filelist[[i]][["read_dist"]]})
    read_dist <- Reduce(function(x, y){merge(x, y, by='Feature', all=T)}, read_dist)
    read_cover <- lapply(1:length(filelist), function(i) {filelist[[i]][["read_cover"]]})
    read_cover <- Reduce(function(x, y){merge(x, y, by='Pos', all=T)}, read_cover)
    return(list(read_dist = read_dist, read_cover = read_cover))
}

files <- list.files(opts$indir, pattern = paste0(opts$suffix, "$"), full.names = TRUE)

reads <- merge_files(files)


# ---------------------------------------------------------------------------- #
features <- c("CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES", "NON_ALIGNED_BASES")
read_dist <- reads$read_dist
rownames(read_dist) <- read_dist$Feature
read_dist <- read_dist[, -1]
read_dist <- data.frame(t(read_dist))
read_dist$Run <- rownames(read_dist)
read_dist$NON_ALIGNED_BASES <- read_dist$PF_BASES - read_dist$PF_ALIGNED_BASES
read_dist <- read_dist[, c('Run', features)]
read_dist <- melt(read_dist, id.vars = c('Run'))
colnames(read_dist)[2] <- "Feature"
read_dist$Feature <- factor(read_dist$Feature, 
                            levels = rev(c("CODING_BASES", "UTR_BASES", "INTRONIC_BASES", "INTERGENIC_BASES", "NON_ALIGNED_BASES")))


height <- ceiling(length(files)/5)
height <- ifelse(height <= 4, 4, height)

pdf(paste0(opts$project_name, "_Picard_RNAmetrics.pdf"), width = 8, height = height)
ggplot(read_dist, aes(x = Run, y = value, fill = Feature)) + 
    geom_bar(stat = "identity", position = "fill") + 
    geom_abline(slope = 0, intercept = c(0.25,0.50,0.75), color = "red") + 
    coord_flip() + 
    theme_bw() + 
    theme(axis.title = element_text(size = 8,face = "bold"), 
          axis.text = element_text(size = 6,face = "bold"), 
          plot.title = element_text(hjust = 0.5, size = 12,face = "bold"), 
          axis.text.x=element_text(angle=45,hjust=1, vjust=1),
          legend.text = element_text(size=6, face="bold")) + 
    labs(title = "Aligned reads distribution(Picard)", x = "", y = "Ratio")
dev.off()


# ---------------------------------------------------------------------------- #
read_cover <- melt(reads$read_cover, id.vars = c('Pos'))
colnames(read_cover)[2] <- "Run"

height <- floor(length(files)/20) + 8

pdf(paste0(opts$project_name, "_Picard_geneCover.pdf"), width = 8, height = height)
ggplot(read_cover, aes(x = Pos, y = value, color = Run)) + 
    geom_line() + 
    theme_bw() + 
    theme(axis.title = element_text(size = 8,face = "bold"), 
          axis.text = element_text(size = 6,face = "bold"), 
          plot.title = element_text(hjust = 0.5, size = 12,face = "bold"), 
          axis.text.x=element_text(angle=45,hjust=1, vjust=1),
          legend.text = element_text(size=6, face="bold"),
          legend.position = "bottom") + 
    labs(title = "Gene coverage(Picard)", x = "Gene body percentile (5'->3')", y = "")
dev.off()

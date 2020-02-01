#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2019-11-04 20:10:05
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

if (!require("optparse", quietly = TRUE)) {
    install.packages("optparse")

    if (!require("optparse", quietly = TRUE)) {
        stop("optparse not installed, please install it first.")
    }
}

suppressPackageStartupMessages(library("optparse"))

opt_list <- list(
    make_option(
        c("-i", "--indir"),
        type = "character",
        help = "Path to infer_experiment.py log files, one sample one file.",
        metavar = "character"
        ),
    make_option(
        c("-o", "--outdir"),
        type = "character",
        help = "Directory of merged file names.",
        metavar = "character"
        ),
    make_option(
        c("--suffix"),
        type = "character",
        help = "Suffix of log files.",
        metavar = "character",
        default = "inferExpr.log"
        ),
    make_option(
        c("--project_name"),
        type = "character",
        help = "Prefix of output files.",
        metavar = "character"
        )
)

opts <- parse_args(OptionParser(option_list = opt_list, 
                                usage = "usage: %prog [options]", 
                                add_help_option = TRUE, 
                                prog = "", 
                                description = "Merge infer_experiment out files. infer_experiment.py -i in.Aligned.sortedByCoord.out.bam -r ref.sorted.bed12 -s 2000000 2>&1 | tee in_infer.log"))

if (is.null(opts$indir)) {
    stop("-i --indir not set")
} else if (is.null(opts$outdir)) {
    stop("-o --outdir not set")
}


multiMerge <- function(path,
                       pattern) {
    files <- file.path(path, list.files(path, pattern = paste0(pattern, "$"), recursive = TRUE))
    fileList <- lapply(files, function(x) {tail(readLines(x), 4)})
    names(fileList) <- sub(paste0("[^0-9a-zA-Z]*", pattern, "$"), "", basename(files))
    fileList <- do.call(rbind, fileList)
    fileList <- as.data.frame(cbind(rownames(fileList),
                     ifelse(fileList[, 1] == "This is SingleEnd Data", "SINGLE", "PAIRED"),
                     stringr::str_match(fileList[, 2], "(.*): (.*)")[, 3],
                     stringr::str_match(fileList[, 3], "(.*)\"(.*)\": (.*)")[, 3:4],
                     stringr::str_match(fileList[, 4], "(.*)\"(.*)\": (.*)")[, 3:4]
                     ), stringsAsFactors = FALSE)
    colnames(fileList) <- c("Run", "Layout", "Failed_ratio", "Type1", "Ratio_of_type1", "Type2", "Ratio_of_type2")
    fileList$Type_diff <- abs(as.numeric(fileList$Ratio_of_type1) - as.numeric(fileList$Ratio_of_type2))
    fileList$Stranded_infer <- apply(fileList, 1, function(x) {
                                    diff <- as.numeric(x[5]) - as.numeric(x[7])
                                    absDiff <- abs(diff)
                                    if (diff >= 0 & absDiff > 0.5) {
                                        return("fr-secondstrand")
                                    } else if (diff < 0 & absDiff > 0.5) {
                                        return("fr-firststrand")
                                    } else if (absDiff < 0.1) {
                                        return("fr-unstranded")
                                    } else {
                                        return("not-sure")
                                    }
                                })
    return(fileList)
}


res <- multiMerge(opts$indir, opts$suffix)

write.table(res,
            file = file.path(opts$outdir, paste0(opts$project_name, "_inferExpr.tsv")),
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

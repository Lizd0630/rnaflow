#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2019-10-11 16:24:10
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
suppressPackageStartupMessages(library("data.table"))

opt_list <- list(
    make_option(
        c("-i", "--indir"),
        type = "character",
        help = "RSEM output directory, only contain the results form the same specie.",
        metavar = "character"
        ),
    make_option(
        c("-o", "--out_dir"),
        type = "character",
        help = "Output directory of merged files.",
        metavar = "character"
        ),
    make_option(
        c("--project_name"),
        type = "character",
        help = "Prefix of merged results.",
        metavar = "character"
        )
)

opts <- parse_args(OptionParser(
    option_list = opt_list, 
    usage = "usage: %prog [options]", 
    add_help_option = TRUE, 
    prog = "", 
    description = "Merge RSEM output files from spesific directory."))

if (is.null(opts$indir)) {
    stop("-i --indir not set")
} else if (is.null(opts$out_dir)) {
    stop("-o --out_dir not set")
} else if (is.null(opts$project_name)) {
    stop("-o --project_name not set")
}


geneMerge <- function(path, 
                      type = "expected_count"){
  filenames <- sort(list.files(path=path, pattern="genes.results", full.names=TRUE))
  datalist <- lapply(filenames, function(x){
    tmp <- fread(x, header=TRUE, sep = "\t")
    tmp <- tmp[, c("gene_id", type), with=F]
    samplename <- sub("[\\._]*genes.results$", "", tail(strsplit(x,"/")[[1]], n=1))
    setnames(tmp, type, samplename)
    setkey(tmp, gene_id)
    return(tmp)})
  all <- Reduce(function(x,y) {merge(x, y, all=T, by="gene_id")}, datalist)
  outname <- file.path(opts$out_dir, paste0(opts$project_name, "_Gene_", type, ".tsv"))
  fwrite(all, file = outname, sep = "\t")
}


isoMerge <- function(path, 
                     type = "expected_count"){
  filenames <- sort(list.files(path=path, pattern="isoforms.results", full.names=TRUE))
  datalist <- lapply(filenames, function(x){
    tmp <- fread(x, header=TRUE, sep = "\t")
    tmp <- tmp[, c("transcript_id", "gene_id", type), with=F]
    samplename <- sub("[\\._]*isoforms.results$", "", tail(strsplit(x,"/")[[1]], n=1))
    setnames(tmp, type, samplename)
    setkey(tmp, transcript_id)
    return(tmp)})
  all <- Reduce(function(x,y) {merge(x, y, all=T, by=c("transcript_id", "gene_id"))}, datalist)
  outname <- file.path(opts$out_dir, paste0(opts$project_name, "_Iso_", type, ".tsv"))
  fwrite(all, file = outname, sep = "\t")
}


geneMerge(opts$indir, type = "expected_count")
geneMerge(opts$indir, type = "TPM")
geneMerge(opts$indir, type = "FPKM")

isoMerge(opts$indir, type = "expected_count")
isoMerge(opts$indir, type = "TPM")
isoMerge(opts$indir, type = "FPKM")


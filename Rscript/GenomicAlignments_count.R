#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Date    : 2019-09-06 09:35:20
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
        c("-i", "--in_dir"),
        type = "character",
        help = "STAR aligned output directory.",
        metavar = "character"
        ),
    make_option(
        c("-o", "--out_dir"),
        type = "character",
        help = "Quantification RData output directory.",
        metavar = "character"
        ),
    make_option(
        c("-c", "--counts_type"),
        type = "character",
        default = "gene",
        help = "Counts type: gene, exon, exonbin.",
        metavar = "character"
        ),
    make_option(
        c("-n", "--project_name"),
        type = "character",
        help = "Quantification RData/matrix output name prefix",
        metavar = "character"
        ),
    make_option(
        c("-m", "--meta_file"),
        type = "character",
        help = "Tab seperated sample information (tsv file), contains at least 5 columns: Run, R1, R2, Layout, strand_specificity.",
        metavar = "character"
        ),
    make_option(
        c("-g", "--ref"),
        type = "character",
        help = "The reference annotation of species, gtf or gff.",
        metavar = "character"
        ),
    make_option(
        c("-t", "--reftype"),
        type = "character",
        default = "gtf",
        help = "The reference type of annotation, gff or gtf, [default \"%default\"].",
        metavar = "gff/gtf"
        ),
    make_option(
        c("-p", "--num_threads"),
        type = "integer",
        default = 8,
        help = "Process threads number, [default %default].",
        metavar = "integer"
        ),
    make_option(
        c("-s", "--bam_suffix"),
        type = "character",
        default = "Aligned.sortedByCoord.out.bam",
        help = "Suffix of bamfile, [default \"%default\"].",
        metavar = "character"
        )
)

opts <- parse_args(
  OptionParser(
    option_list = opt_list,
    usage = "usage: %prog [options]",
    add_help_option = TRUE,
    prog = "GenomicAlimenstsQuant",
    description = "The Rscript to quantification gene and exon counts. Reads for SINGLE end, Fragments for PAIRED end.")
  )

if (is.null(opts$in_dir)) {
    stop("-i --in_dir not set")
} else if (is.null(opts$out_dir)) {
    stop("-o --out_dir not set")
} else if (is.null(opts$project_name)) {
    stop("-n --project_name not set")
} else if (is.null(opts$meta_file)) {
    stop("-m --meta_file not set")
} else if (is.null(opts$ref)) {
    stop("-r --ref not set")
} else if (is.null(opts$reftype)) {
    stop("-t --reftype not set")
} else if (is.null(opts$num_threads)) {
    stop("-p --num_threads not set")
}

if (!dir.exists(opts$out_dir)) {
        dir.create(opts$out_dir)
    if (!dir.exists(opts$out_dir)) {
        stop("Can't create output directory.")
    }
}

suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomicAlignments"))


#-----------------------------------------------------------------------------#
filterfiles <- function(meta, ids, path, bamfiles) {
  if (length(ids) > 0){
    sampleTable <- DataFrame(meta[ids,])
    bamfile <- NULL
    for(i in 1:length(ids)){
        bamfile[i] <- grep(paste0(ids[i], "[^0-9a-zA-Z]"), bamfiles, value=TRUE)
    }
    bamfile <- file.path(path, bamfile)
    return(list(sampleTable=sampleTable, bamfile=bamfile))
  }
}

merge_se <- function(x, y) {
  if (!is.null(x) & !is.null(y)) {
    rowranges <- rowRanges(x)
    x_cnt <- assay(x)
    y_cnt <- assay(y)
    x_coldata <- colData(x)
    y_coldata <- colData(y)
    se <- SummarizedExperiment(
      assays = SimpleList(counts = cbind(x_cnt, y_cnt)), 
      rowRanges = rowranges, 
      colData = rbind(x_coldata, y_coldata))
  } else if (!is.null(x)) {
      se <- x
  } else if (!is.null(y)) {
      se <- y
  } else {
      se <- NULL
  }
  return(se)
}

flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)

bamfiles <- list.files(path = file.path(opts$in_dir), pattern = paste0("*", opts$bam_suffix, "$"))
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
cat(date(), "make genome features objects.\n")

## meta information table
meta <- read.table(opts$meta_file, header=T, sep="\t", stringsAsFactors = FALSE)
if (all(c("Run", "R1", "R2", "Layout", "Strand_specificity") %in% colnames(meta))) {
    rownames(meta) <- meta$Run
} else {
    stop("colnames error in metainfo!")
}

refpath <- file.path(opts$ref)
txdb <- makeTxDbFromGFF(
  refpath,
  format = tolower(opts$reftype),
  dataSource = basename(opts$ref)
  )

ebg <- exonsBy(txdb, by="gene")
#-----------------------------------------------------------------------------#


cat(date(), "reads counting.\n")

if (TRUE) {
    options(mc.cores=opts$num_threads)
}

#########
## in summarizeOverlaps, ignore.strand has higher priority than strandMode
#########

if (opts$counts_type == "gene") {
  #-----------------------------------------------------------------------------#
  ## single end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  SE_0_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-unstranded")$Run

  if (length(SE_0_id) > 0){
    files <- filterfiles(meta, SE_0_id, opts$in_dir, bamfiles)
    SE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    ## gene reads count
    SE_0_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = SE_0_bamlist,
                                      ignore.strand = TRUE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_0_geneCnt) <- SE_0_id
    colData(SE_0_geneCnt) <- files$sampleTable
  } else {
    SE_0_geneCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  SE_1_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-secondstrand")$Run

  if (length(SE_1_id) > 0){
    files <- filterfiles(meta, SE_1_id, opts$in_dir, bamfiles)
    SE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    SE_1_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = SE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_1_geneCnt) <- SE_1_id
      colData(SE_1_geneCnt) <- files$sampleTable
  } else {
    SE_1_geneCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  SE_2_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-firststrand")$Run

  if (length(SE_2_id) > 0){
    files <- filterfiles(meta, SE_2_id, opts$in_dir, bamfiles)
    SE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    SE_2_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = SE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = invertStrand)
    colnames(SE_2_geneCnt) <- SE_2_id
    colData(SE_2_geneCnt) <- files$sampleTable
  } else {
    SE_2_geneCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## PAIRED end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  PE_0_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-unstranded")$Run

  if (length(PE_0_id) > 0){
    files <- filterfiles(meta, PE_0_id, opts$in_dir, bamfiles)
    PE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_0_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = PE_0_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 0,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_0_geneCnt) <- PE_0_id
    colData(PE_0_geneCnt) <- files$sampleTable
  } else {
    PE_0_geneCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  PE_1_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-secondstrand")$Run

  if (length(PE_1_id) > 0){
    files <- filterfiles(meta, PE_1_id, opts$in_dir, bamfiles)
    PE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_1_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = PE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 1,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_1_geneCnt) <- PE_1_id
    colData(PE_1_geneCnt) <- files$sampleTable
  } else {
    PE_1_geneCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  PE_2_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-firststrand")$Run

  if (length(PE_2_id) > 0){
    files <- filterfiles(meta, PE_2_id, opts$in_dir, bamfiles)
    PE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_2_geneCnt <- summarizeOverlaps(features = ebg,
                                      mode = "Union",
                                      reads = PE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 2,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_2_geneCnt) <- PE_2_id
    colData(PE_2_geneCnt) <- files$sampleTable
  } else {
    PE_2_geneCnt <- NULL
  }

  se_gene <- merge_se(merge_se(SE_0_geneCnt, SE_1_geneCnt), SE_2_geneCnt)
  pe_gene <- merge_se(merge_se(PE_0_geneCnt, PE_1_geneCnt), PE_2_geneCnt)
  gene_se <- merge_se(se_gene, pe_gene)

  saveRDS(gene_se, file = file.path(opts$out_dir,paste0(opts$project_name,".GeneCnt.Rds")))

} else if (opts$counts_type == "exonBin") {
  ep <- exonicParts(txdb, linked.to.single.gene.only=FALSE)
  #-----------------------------------------------------------------------------#
  ## single end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  SE_0_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-unstranded")$Run

  if (length(SE_0_id) > 0){
    files <- filterfiles(meta, SE_0_id, opts$in_dir, bamfiles)
    SE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    SE_0_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = SE_0_bamlist,
                                      ignore.strand = TRUE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_0_exonBinCnt) <- SE_0_id
    colData(SE_0_exonBinCnt) <- files$sampleTable
  } else {
    SE_0_exonBinCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  SE_1_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-secondstrand")$Run

  if (length(SE_1_id) > 0){
    files <- filterfiles(meta, SE_1_id, opts$in_dir, bamfiles)
    SE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    SE_1_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = SE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_1_exonBinCnt) <- SE_1_id
    colData(SE_1_exonBinCnt) <- files$sampleTable
  } else {
    SE_1_exonBinCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  SE_2_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-firststrand")$Run

  if (length(SE_2_id) > 0){
    files <- filterfiles(meta, SE_2_id, opts$in_dir, bamfiles)
    SE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    SE_2_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = SE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = invertStrand)
    colnames(SE_2_exonBinCnt) <- SE_2_id
    colData(SE_2_exonBinCnt) <- files$sampleTable
  } else {
    SE_2_exonBinCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## PAIRED end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  PE_0_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-unstranded")$Run

  if (length(PE_0_id) > 0){
    files <- filterfiles(meta, PE_0_id, opts$in_dir, bamfiles)
    PE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    PE_0_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = PE_0_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 0,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_0_exonBinCnt) <- PE_0_id
    colData(PE_0_exonBinCnt) <- files$sampleTable
  } else {
    PE_0_exonBinCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  PE_1_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-secondstrand")$Run

  if (length(PE_1_id) > 0){
    files <- filterfiles(meta, PE_1_id, opts$in_dir, bamfiles)
    PE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    PE_1_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = PE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 1,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_1_exonBinCnt) <- PE_1_id
    colData(PE_1_exonBinCnt) <- files$sampleTable
  } else {
    PE_1_exonBinCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  PE_2_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-firststrand")$Run

  if (length(PE_2_id) > 0){
    files <- filterfiles(meta, PE_2_id, opts$in_dir, bamfiles)
    PE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)
    ## exon reads count
    PE_2_exonBinCnt <- summarizeOverlaps(features = ep,
                                      mode = "Union",
                                      reads = PE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 2,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_2_exonBinCnt) <- PE_2_id
    colData(PE_2_exonBinCnt) <- files$sampleTable
  } else {
    PE_2_exonBinCnt <- NULL
  }

  se_exonBin <- merge_se(merge_se(SE_0_exonBinCnt, SE_1_exonBinCnt), SE_2_exonBinCnt)
  pe_exonBin <- merge_se(merge_se(PE_0_exonBinCnt, PE_1_exonBinCnt), PE_2_exonBinCnt)
  exonBin_se <- merge_se(se_exonBin, pe_exonBin)

  saveRDS(exonBin_se, file = file.path(opts$out_dir,paste0(opts$project_name,".ExonBinCnt.Rds")))

} else if (opts$counts_type == "exon") {
  exon <- exons(txdb, use.names = TRUE)
  #-----------------------------------------------------------------------------#
  ## single end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  SE_0_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-unstranded")$Run

  if (length(SE_0_id) > 0){
    files <- filterfiles(meta, SE_0_id, opts$in_dir, bamfiles)
    SE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    ## gene reads count
    SE_0_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = SE_0_bamlist,
                                      ignore.strand = TRUE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_0_exonCnt) <- SE_0_id
    colData(SE_0_exonCnt) <- files$sampleTable
  } else {
    SE_0_exonCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  SE_1_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-secondstrand")$Run

  if (length(SE_1_id) > 0){
    files <- filterfiles(meta, SE_1_id, opts$in_dir, bamfiles)
    SE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    SE_1_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = SE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(SE_1_exonCnt) <- SE_1_id
      colData(SE_1_exonCnt) <- files$sampleTable
  } else {
    SE_1_exonCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## single end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  SE_2_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-firststrand")$Run

  if (length(SE_2_id) > 0){
    files <- filterfiles(meta, SE_2_id, opts$in_dir, bamfiles)
    SE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    SE_2_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = SE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = TRUE,
                                      param = sbp,
                                      preprocess.reads = invertStrand)
    colnames(SE_2_exonCnt) <- SE_2_id
    colData(SE_2_exonCnt) <- files$sampleTable
  } else {
    SE_2_exonCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## PAIRED end & non strand-specific samples
  #-----------------------------------------------------------------------------#
  PE_0_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-unstranded")$Run

  if (length(PE_0_id) > 0){
    files <- filterfiles(meta, PE_0_id, opts$in_dir, bamfiles)
    PE_0_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_0_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = PE_0_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 0,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_0_exonCnt) <- PE_0_id
    colData(PE_0_exonCnt) <- files$sampleTable
  } else {
    PE_0_exonCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-secondstrand samples
  #-----------------------------------------------------------------------------#
  PE_1_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-secondstrand")$Run

  if (length(PE_1_id) > 0){
    files <- filterfiles(meta, PE_1_id, opts$in_dir, bamfiles)
    PE_1_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_1_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = PE_1_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 1,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_1_exonCnt) <- PE_1_id
    colData(PE_1_exonCnt) <- files$sampleTable
  } else {
    PE_1_exonCnt <- NULL
  }
  #-----------------------------------------------------------------------------#
  ## paired end & fr-firststrand samples
  #-----------------------------------------------------------------------------#
  PE_2_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-firststrand")$Run

  if (length(PE_2_id) > 0){
    files <- filterfiles(meta, PE_2_id, opts$in_dir, bamfiles)
    PE_2_bamlist <- BamFileList(files$bamfile,yieldSize=2000000)

    PE_2_exonCnt <- summarizeOverlaps(features = exon,
                                      mode = "Union",
                                      reads = PE_2_bamlist,
                                      ignore.strand = FALSE,
                                      inter.feature = FALSE,
                                      singleEnd = FALSE,
                                      fragments = FALSE,
                                      strandMode = 2,
                                      param = sbp,
                                      preprocess.reads = NULL)
    colnames(PE_2_exonCnt) <- PE_2_id
    colData(PE_2_exonCnt) <- files$sampleTable
  } else {
    PE_2_exonCnt <- NULL
  }

  se_exon <- merge_se(merge_se(SE_0_exonCnt, SE_1_exonCnt), SE_2_exonCnt)
  pe_exon <- merge_se(merge_se(PE_0_exonCnt, PE_1_exonCnt), PE_2_exonCnt)
  exon_se <- merge_se(se_exon, pe_exon)

  saveRDS(exon_se, file = file.path(opts$out_dir,paste0(opts$project_name,".exonCnt.Rds")))
} else {
  cat("Not ready!")
}

cat(date(), "all done!\n")

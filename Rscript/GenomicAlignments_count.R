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

opts <- parse_args(OptionParser(option_list = opt_list, 
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

if (TRUE) {
    options(mc.cores=opts$num_threads)
}

#-----------------------------------------------------------------------------#
cat(date(), ", make genome features objects.\n")

## meta information table
meta <- read.table(opts$meta_file, header=T, sep="\t", stringsAsFactors = FALSE)
if (all(c("Run", "R1", "R2", "Layout", "Strand_specificity") %in% colnames(meta))) {
    rownames(meta) <- meta$Run
} else {
    stop("colnames error in metainfo!")
}

suppressPackageStartupMessages(library("GenomicFeatures"))
suppressPackageStartupMessages(library("GenomicAlignments"))

refpath <- file.path(opts$ref)
bamfile <- list.files(path = file.path(opts$in_dir), pattern = paste0("*", opts$bam_suffix, "$"))
print(bamfile)

## parse gtf or gff
txdb <- makeTxDbFromGFF(refpath,
                        format = tolower(opts$reftype),
                        dataSource = basename(opts$ref)
                        )

## gene level annotation
ebg <- exonsBy(txdb, by="gene")
## exon level annotation
ep <- exonicParts(txdb, linked.to.single.gene.only=FALSE)
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isDuplicate = FALSE,
                    isUnmappedQuery = FALSE)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)
#-----------------------------------------------------------------------------#



cat(date(), ", reads counting ...\n")
#-----------------------------------------------------------------------------#
## single end & non strand-specific samples
#-----------------------------------------------------------------------------#
SE_0_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-unstranded")$Run

if (length(SE_0_id) > 0){
    sampleTable <- DataFrame(meta[SE_0_id,])
    SE_0_bamfile <- NULL
    for(i in 1:length(SE_0_id)){
        SE_0_bamfile[i] <- grep(paste0(SE_0_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    SE_0_bamfile <- file.path(opts$in_dir, SE_0_bamfile)
    print(SE_0_bamfile)
    ## bamfilelist
    SE_0_bamlist <- BamFileList(SE_0_bamfile,yieldSize=2000000)
    ## exon reads count
    SE_0_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = SE_0_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      # strandMode = 0, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_0_exonCnt) <- SE_0_id
    colData(SE_0_exonCnt) <- sampleTable

    ## gene reads count
    SE_0_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = SE_0_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      # strandMode = 0, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_0_geneCnt) <- SE_0_id
    colData(SE_0_geneCnt) <- sampleTable
} else {
    SE_0_exonCnt <- NULL
    SE_0_geneCnt <- NULL
}



#-----------------------------------------------------------------------------#
## single end & fr-secondstrand samples
#-----------------------------------------------------------------------------#
SE_1_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-secondstrand")$Run

if (length(SE_1_id) > 0){
    sampleTable <- DataFrame(meta[SE_1_id,])
    SE_1_bamfile <- NULL
    for(i in 1:length(SE_1_id)){
        SE_1_bamfile[i] <- grep(paste0(SE_1_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    SE_1_bamfile <- file.path(opts$in_dir,SE_1_bamfile)
    ## bamfilelist
    SE_1_bamlist <- BamFileList(SE_1_bamfile,yieldSize=2000000)
    ## exon reads count
    SE_1_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = SE_1_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 1, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_1_exonCnt) <- SE_1_id
    colData(SE_1_exonCnt) <- sampleTable

    ## gene reads count
    SE_1_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = SE_1_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 1, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_1_geneCnt) <- SE_1_id
    colData(SE_1_geneCnt) <- sampleTable
} else {
    SE_1_exonCnt <- NULL
    SE_1_geneCnt <- NULL
}



#-----------------------------------------------------------------------------#
## single end & fr-firststrand samples
#-----------------------------------------------------------------------------#
SE_2_id <- subset(meta, Layout=="SINGLE" & Strand_specificity == "fr-firststrand")$Run

if (length(SE_2_id) > 0){
    sampleTable <- DataFrame(meta[SE_2_id,])
    SE_2_bamfile <- NULL
    for(i in 1:length(SE_2_id)){
        SE_2_bamfile[i] <- grep(paste0(SE_2_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    SE_2_bamfile <- file.path(opts$in_dir,SE_2_bamfile)
    ## bamfilelist
    SE_2_bamlist <- BamFileList(SE_2_bamfile,yieldSize=2000000)
    ## exon reads count
    SE_2_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = SE_2_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 2, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_2_exonCnt) <- SE_2_id
    colData(SE_2_exonCnt) <- sampleTable

    ## gene reads count
    SE_2_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = SE_2_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 2, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(SE_2_geneCnt) <- SE_2_id
    colData(SE_2_geneCnt) <- sampleTable
} else {
    SE_2_exonCnt <- NULL
    SE_2_geneCnt <- NULL
}




#-----------------------------------------------------------------------------#
## PAIRED end & non strand-specific samples
#-----------------------------------------------------------------------------#
PE_0_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-unstranded")$Run

if (length(PE_0_id) > 0){
    sampleTable <- DataFrame(meta[PE_0_id,])
    PE_0_bamfile <- NULL
    for(i in 1:length(PE_0_id)){
        PE_0_bamfile[i] <- grep(paste0(PE_0_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    PE_0_bamfile <- file.path(opts$in_dir, PE_0_bamfile)
    ## bamfilelist
    PE_0_bamlist <- BamFileList(PE_0_bamfile,yieldSize=2000000)
    ## exon reads count
    PE_0_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = PE_0_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = FALSE,
                                      fragments = FALSE, 
                                      strandMode = 0, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_0_exonCnt) <- PE_0_id
    colData(PE_0_exonCnt) <- sampleTable

    ## gene reads count
    PE_0_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = PE_0_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 0, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_0_geneCnt) <- PE_0_id
    colData(PE_0_geneCnt) <- sampleTable
} else {
    PE_0_exonCnt <- NULL
    PE_0_geneCnt <- NULL
}



#-----------------------------------------------------------------------------#
## paired end & fr-secondstrand samples
#-----------------------------------------------------------------------------#
PE_1_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-secondstrand")$Run

if (length(PE_1_id) > 0){
    sampleTable <- DataFrame(meta[PE_1_id,])
    PE_1_bamfile <- NULL
    for(i in 1:length(PE_1_id)){
        PE_1_bamfile[i] <- grep(paste0(PE_1_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    PE_1_bamfile <- file.path(opts$in_dir,PE_1_bamfile)
    ## bamfilelist
    PE_1_bamlist <- BamFileList(PE_1_bamfile,yieldSize=2000000)
    ## exon reads count
    PE_1_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = PE_1_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = FALSE,
                                      fragments = FALSE, 
                                      strandMode = 1, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_1_exonCnt) <- PE_1_id
    colData(PE_1_exonCnt) <- sampleTable

    ## gene reads count
    PE_1_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = PE_1_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 1, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_1_geneCnt) <- PE_1_id
    colData(PE_1_geneCnt) <- sampleTable
} else {
    PE_1_exonCnt <- NULL
    PE_1_geneCnt <- NULL
}



#-----------------------------------------------------------------------------#
## single end & fr-firststrand samples
#-----------------------------------------------------------------------------#
PE_2_id <- subset(meta, Layout=="PAIRED" & Strand_specificity == "fr-firststrand")$Run

if (length(PE_2_id) > 0){
    sampleTable <- DataFrame(meta[PE_2_id,])
    PE_2_bamfile <- NULL
    for(i in 1:length(PE_2_id)){
        PE_2_bamfile[i] <- grep(paste0(PE_2_id[i], "[^0-9a-zA-Z]"), bamfile, value=T)
    }
    PE_2_bamfile <- file.path(opts$in_dir,PE_2_bamfile)
    ## bamfilelist
    PE_2_bamlist <- BamFileList(PE_2_bamfile,yieldSize=2000000)
    ## exon reads count
    PE_2_exonCnt <- summarizeOverlaps(features = ep, 
                                      mode = "IntersectionStrict",
                                      reads = PE_2_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 2, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_2_exonCnt) <- PE_2_id
    colData(PE_2_exonCnt) <- sampleTable

    ## gene reads count
    PE_2_geneCnt <- summarizeOverlaps(features = ebg, 
                                      mode = "Union",
                                      reads = PE_2_bamlist,
                                      ignore.strand = FALSE, 
                                      inter.feature = FALSE, 
                                      singleEnd = TRUE,
                                      fragments = FALSE, 
                                      strandMode = 2, 
                                      param = sbp, 
                                      preprocess.reads = NULL)
    colnames(PE_2_geneCnt) <- PE_2_id
    colData(PE_2_geneCnt) <- sampleTable
} else {
    PE_2_exonCnt <- NULL
    PE_2_geneCnt <- NULL
}


merge.se <- function(x, y) {
    if ((exists(x) & !is.null(get(x))) & (exists(y) & !is.null(get(y)))) {
        rowranges <- rowRanges(get(x))
        x_cnt <- assay(get(x))
        y_cnt <- assay(get(y))
        x_coldata <- colData(get(x))
        y_coldata <- colData(get(y))
        se <- SummarizedExperiment(assays = cbind(x_cnt, y_cnt), rowRanges = rowranges, colData = rbind(x_coldata, y_coldata))
    } else if (exists(x) & !is.null(get(x))) {
        # rowranges <- rowRanges(get(x))
        # x_cnt <- assay(get(x))
        # se <- SummarizedExperiment(assays = x_cnt, rowRanges = rowranges, colData = x_coldata)
        se <- get(x)
    } else if (exists(y) & !is.null(get(y))) {
        # rowranges <- rowRanges(get(y))
        # x_cnt <- assay(get(y))
        # se <- SummarizedExperiment(assays = y_cnt, rowRanges = rowranges, colData = y_coldata)
        se <- get(y)
    } else {
        se <- NULL
    }
    return(se)
}

se.gene <- merge.se(merge.se("SE_0_geneCnt", "SE_1_geneCnt"), "SE_2_geneCnt")
pe.gene <- merge.se(merge.se("PE_0_geneCnt", "PE_1_geneCnt"), "PE_2_geneCnt")
gene_se <- merge.se("se.gene", "pe.gene")

se.exon <- merge.se(merge.se("SE_0_exonCnt", "SE_1_exonCnt"), "SE_2_exonCnt")
pe.exon <- merge.se(merge.se("PE_0_exonCnt", "PE_1_exonCnt"), "PE_2_exonCnt")
exon_se <- merge.se("se.exon", "pe.exon")

saveRDS(list(gene_se=gene_se, exon_se=exon_se), file = file.path(opts$out_dir,paste0(opts$project_name,".gene_exon.RawCnt.Rds")))

cat(date(), ", all done!\n")
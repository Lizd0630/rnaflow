#!/usr/bin/env python
# -*- coding: utf-8 -*-


import json
import pathlib as pl

# ENV
__dir__ = str(pl.Path(__file__).resolve().parents[1])
VERSION = "2.0.1"
LABEL = "RNA seq pipeline"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
READ1 = r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[1]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz"
READ2 = r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[2]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz"

RSEM_STRAND = { "fr-unstranded": "none",
                "fr-firststrand": "reverse",
                "fr-secondstrand": "forward"}
PICARD_STRAND = {   "fr-unstranded": "NONE",
                    "fr-firststrand": "FIRST_READ_TRANSCRIPTION_STRAND",
                    "fr-secondstrand": "SECOND_READ_TRANSCRIPTION_STRAND"}
RESULTS_SUFFIX = {  "infer_expr": "inferExpr.log",
                    "CollectRnaSeqMetrics": "RNAmetrics.txt",
                    "geneBody_coverage": "geneBodyCoverage.txt"}

# softwares path
SOFTWARES = pl.Path(__dir__, "config/softwares.json")


# software parameters
FASTQC_CONFIG = pl.Path(__dir__, "config/fastqc.json")
FASTP_CONFIG = pl.Path(__dir__, "config/fastp.json")
TRIMMOMATIC_CONFIG = pl.Path(__dir__, "config/trimmomatic.json")
STAR_CONFIG = pl.Path(__dir__, "config/STAR_align.json")
RSEM_CONFIG = pl.Path(__dir__, "config/RSEM_quant.json")
BAMQC_CONFIG = pl.Path(__dir__, "config/bamqc.json")


# bash, R scripts
GA_SCRIPT = pl.Path(__dir__, "Rscript/GenomicAlignments_count.R")
STAR_STAT = pl.Path(__dir__, "py/StarStat.py")
RSEM_MERGE = pl.Path(__dir__, "Rscript/RSEM_merge_results.R")
INFEREXPR_MERGE = pl.Path(__dir__, "Rscript/inferExpr_merge.R")
RESULTS_RNAMETRICS = pl.Path(__dir__, "Rscript/collectRnaMetrics_plot.R")
RESULTS_GENEBODY = pl.Path(__dir__, "Rscript/geneBodyCoverage_plot.R")


class ParseDict:
    def __init__(self, json_file):
        with open(json_file, 'r') as f:
            temp = json.loads(f.read())
        self.dict = temp

    def make_config(self, sep=" "):
        param = ' '.join([f"{k}{sep}{v}" for k, v in self.dict.items()])
        return(param)

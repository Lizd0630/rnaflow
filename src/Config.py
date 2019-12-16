#!/usr/bin/env python
# -*- coding: utf-8 -*-


import json
import pathlib as pl


__dir__ = str(pl.Path(__file__).resolve().parents[1])
VERSION = "0.1.0"
LABEL = "RNA seq pipeline"
CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
SOFTWARES = pl.Path(__dir__, "config/software.json")
FASTQC_CONFIG = pl.Path(__dir__, "config/fastqc.json")
FASTP_CONFIG = pl.Path(__dir__, "config/fastp.json")
TRIMMOMATIC_CONFIG = pl.Path(__dir__, "config/trimmomatic.json")
STAR_CONFIG = pl.Path(__dir__, "config/STAR_align.json")
RSEM_CONFIG = pl.Path(__dir__, "config/RSEM_quant.json")
GA_CONFIG = pl.Path(__dir__, "config/GenomicAlignments.json")


RSEM_STRAND = {"fr-unstranded": "none", "fr-firststrand": "reverse", "fr-secondstrand": "forward"}


class Param:
    def __init__(self, json_file):
        self.file = json_file

    def make_config(self):
        with open(self.file, 'r') as f:
            temp = json.loads(f.read())
        param = " ".join([f"{k} {v}" for k, v in temp.items()])
        return(param)


class Soft:
    def __init__(self, json_file):
        self.file = json_file

    def make_config(self):
        with open(self.file, 'r') as f:
            soft = json.loads(f.read())
        return(soft)

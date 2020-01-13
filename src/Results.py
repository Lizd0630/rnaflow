#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-20 11:44:22
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

import pathlib as pl
import re
from src.Config import STAR_STAT, RSEM_MERGE, INFEREXPR_MERGE, RESULTS_RNAMETRICS, RESULTS_GENEBODY


class Results:
    def __init__(self, input_dir, output_dir, project_name):
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.proj = project_name

    # def star(self):
    #     cmds = []
    #     cmd = f"bash {STAR_STAT}\
    #             -i {self.input_dir} \
    #             > {self.output_dir}{self.proj}_STARfinal.tsv"
    #     cmd = re.sub(" {2,}", " ", cmd)
    #     cmds.append(cmd)
    #     return(cmds)

    def star(self):
        cmds = []
        cmd = f"python3 {STAR_STAT} \
                -i {self.input_dir} \
                -o {self.output_dir} \
                -p {self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def rsem(self):
        cmds = []
        cmd = f"Rscript {RSEM_MERGE} \
                -i {self.input_dir} \
                -o {self.output_dir} \
                --project_name {self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def infer_expr(self, suffix):
        cmds = []
        cmd = f"Rscript {INFEREXPR_MERGE} \
                -i {self.input_dir} \
                -o {self.output_dir} \
                --suffix {suffix} \
                --project_name {self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def RNAmetrics(self, suffix):
        cmds = []
        cmd = f"Rscript {RESULTS_RNAMETRICS} \
                -i {self.input_dir} \
                -o {self.output_dir} \
                --suffix {suffix} \
                --project_name {self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def geneBody(self, suffix):
        cmds = []
        cmd = f"Rscript {RESULTS_GENEBODY} \
                -i {self.input_dir} \
                -o {self.output_dir} \
                --suffix {suffix} \
                --project_name {self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

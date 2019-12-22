#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-20 11:44:22
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

import pathlib as pl
import re
from src.Config import STAR_STAT, RSEM_MERGE, INFEREXPR_MERGE


class Results:
    def __init__(self, input_dir, output_dir, project_name):
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.proj = project_name

    def star(self):
        cmds = []
        cmd = f"bash {STAR_STAT}\
                -i {self.input_dir} \
                > {self.output_dir}{self.proj}_STARfinal.tsv"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def rsem(self):
        cmds = []
        cmd = f"bash {RSEM_MERGE}\
                -i {self.input_dir} \
                -o {self.output_dir}{self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

    def infer_expr(self):
        cmds = []
        cmd = f"bash {INFEREXPR_MERGE}\
                -i {self.input_dir} \
                -p inferExpr.log \
                -o {self.output_dir}{self.proj}"
        cmd = re.sub(" {2,}", " ", cmd)
        cmds.append(cmd)
        return(cmds)

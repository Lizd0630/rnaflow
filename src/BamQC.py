#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-20 19:12:35
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

from src.Meta import Meta
from src.Config import ParseDict
from src.Config import PICARD_STRAND
import pathlib as pl
import re
import src.MyError as me


class BamQC:
    def __init__(self, meta_file, input_dir, output_dir, softwares, config, suffix, bed12, refFlat):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        self.meta_info = Meta(meta_file).get_info()
        self.input_dir = str(pl.Path(input_dir)) + r"/"
        self.output_dir = str(pl.Path(output_dir)) + r"/"
        self.soft_path = ParseDict(softwares).dict
        self.config = ParseDict(config).dict
        self.suffix = suffix
        self.bed12 = bed12
        self.refFlat = refFlat

    # ----------------------------------RSeQC--------------------------------- #
    def infer_expr(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*")]
        for meta in self.meta_info:
            pattern = re.compile(self.input_dir + meta["run"] + r".*?" + self.suffix + r"$")
            inbam = list(filter(lambda x: pattern.match(x) is not None, allfiles))
            if len(inbam) == 1:
                inbam = inbam[0]
            else:
                raise me.CntError(f"Bam file number error: {meta['run']}")
            cmd = f"infer_experiment.py \
                    {' '.join([f'{k} {v}' for k, v in self.config['infer_expr'].items()])} \
                    -i {inbam} \
                    -r {self.bed12} \
                    2>&1 | tee {self.output_dir}{meta['run']}_inferExpr.log"
            cmd = re.sub(" {2,}", " ", cmd)
            cmds.append(cmd)
        return(cmds)

    def genebody(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*")]
        for meta in self.meta_info:
            pattern = re.compile(self.input_dir + meta["run"] + r".*?" + self.suffix + r"$")
            inbam = list(filter(lambda x: pattern.match(x) is not None, allfiles))
            if len(inbam) == 1:
                inbam = inbam[0]
            else:
                raise me.CntError(f"Bam file number error: {meta['run']}")
            cmd = f"geneBody_coverage.py \
                    {' '.join([f'{k} {v}' for k, v in self.config['genebody'].items()])} \
                    -i {inbam} \
                    -r {self.bed12} \
                    -o {self.output_dir}{meta['run']} \
                    2>&1 | tee {self.output_dir}{meta['run']}_geneBody.log"
            cmd = re.sub(" {2,}", " ", cmd)
            cmds.append(cmd)
        return(cmds)

    # ----------------------------------picard-------------------------------- #
    def rnametrics(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*")]
        for meta in self.meta_info:
            pattern = re.compile(self.input_dir + meta["run"] + r".*?" + self.suffix + r"$")
            inbam = list(filter(lambda x: pattern.match(x) is not None, allfiles))
            if len(inbam) == 1:
                inbam = inbam[0]
            else:
                raise me.CntError(f"Bam file number error: {meta['run']}")
            cmd = f"java -jar {self.soft_path['picard']} \
                    CollectRnaSeqMetrics \
                    {' '.join([f'{k}={v}' for k, v in self.config['rnametrics'].items()])} \
                    I={inbam} \
                    O={self.output_dir}{meta['run']}_RNAmetrics.txt \
                    REF_FLAT={self.refFlat} \
                    STRAND={PICARD_STRAND[meta['strand']]}"
            cmd = re.sub(" {2,}", " ", cmd)
            cmds.append(cmd)
        return(cmds)

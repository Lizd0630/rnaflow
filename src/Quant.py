#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-16 01:03:33
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


from src.Meta import Meta
from src.Config import ParseDict
from src.Config import RSEM_STRAND
import pathlib as pl
import re
import src.MyError as me


class Quant:
    def __init__(self, meta_file, input_dir, output_dir, tool, softwares, config, suffix, ref):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        self.meta_info = Meta(meta_file).get_info()
        self.input_dir = str(pl.Path(input_dir)) + r"/"
        self.output_dir = str(pl.Path(output_dir)) + r"/"
        self.software = tool
        self.soft_path = ParseDict(softwares).dict
        self.config = ParseDict(config).make_config()
        self.suffix = suffix
        self.ref = ref

    def make_cmds(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*")]
        for meta in self.meta_info:
            if meta["layout"] == "PAIRED":
                pattern = re.compile(self.input_dir + meta["run"] + r".*?" + self.suffix)
                inbam = list(filter(lambda x: pattern.match(x) is not None, allfiles))
                if len(inbam) == 1:
                    inbam = inbam[0]
                else:
                    raise me.CntError(f"Bam file number error: {meta['run']}")
                if self.software == "RSEM":
                    cmd = f"{self.soft_path['RSEM']} \
                            --paired-end \
                            --strandedness {RSEM_STRAND[meta['strand']]} \
                            {self.config} \
                            {inbam} \
                            {self.ref} \
                            {self.output_dir}{meta['run']}."
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    pass
            elif meta["layout"] == "SINGLE":
                pattern = re.compile(self.input_dir + meta["run"] + r".*?" + self.suffix)
                inbam = list(filter(lambda x: pattern.match(x) is not None, allfiles))
                if len(inbam) == 1:
                    inbam = inbam[0]
                else:
                    raise me.CntError(f"Bam file number error: {meta['run']}")
                if self.software == "RSEM":
                    cmd = f"{self.soft_path['RSEM']} \
                            --strandedness {RSEM_STRAND[meta['strand']]} \
                            {self.config} \
                            {inbam} \
                            {self.ref} \
                            {self.output_dir}{meta['run']}."
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    pass
            else:
                print("'Layout' Error!")
        return(cmds)

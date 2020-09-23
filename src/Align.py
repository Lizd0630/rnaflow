#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-16 00:20:53
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


from src.Meta import Meta
from src.Config import ParseDict
from src.Config import READ1, READ2
import pathlib as pl
import re
import src.MyError as me


class Align:
    def __init__(self, meta_file, input_dir, output_dir, tool, ref, softwares, config):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        self.meta_info = Meta(meta_file).get_info()
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.software = tool
        self.ref = ref
        self.soft_path = ParseDict(softwares).dict
        self.config = ParseDict(config).make_config()

    def make_cmds(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*gz")]
        for meta in self.meta_info:
            if meta["layout"] == "PAIRED":
                pattern1 = re.compile(self.input_dir + meta["r1"] + READ1 + r"$")
                pattern2 = re.compile(self.input_dir + meta["r2"] + READ2 + r"$")
                R1 = list(filter(lambda x: pattern1.match(x) is not None, allfiles))
                if len(R1) == 1:
                    R1 = R1[0]
                else:
                    raise me.CntError(f"Fastq file number error: {meta['r1']}")
                    # try:
                    #     raise me.CntError(f"fastq file number error: {meta['r1']}")
                    # except me.CntError as ce:
                    #     print(ce.message)
                R2 = list(filter(lambda x: pattern2.match(x) is not None, allfiles))
                if len(R2) == 1:
                    R2 = R2[0]
                else:
                    raise me.CntError(f"Fastq file number error: {meta['r1']}")
                    # try:
                    #     raise me.CntError(f"fastq file number error: {meta['r2']}")
                    # except me.CntError as ce:
                    #     print(ce.message)
                if self.software == "STAR":
                    cmd = f"{self.soft_path['STAR']} \
                            --genomeDir {self.ref} {self.config} \
                            --readFilesIn {R1} {R2} \
                            --outFileNamePrefix {self.output_dir}{meta['run']}. \
                            2>&1 | tee {self.output_dir}{meta['run']}.STAR.log \
                            && {self.soft_path['samtools']} \
                            index {self.output_dir}{meta['run']}.Aligned.sortedByCoord.out.bam"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    pass
            elif meta["layout"] == "SINGLE":
                pattern1 = re.compile(self.input_dir + meta["r1"] + READ1 + r"$")
                R1 = list(filter(lambda x: pattern1.match(x) is not None, allfiles))
                if len(R1) == 1:
                    R1 = R1[0]
                else:
                    raise me.CntError(f"Fastq file number error: {meta['r1']}")
                    # try:
                    #     raise me.CntError(f"fastq file number error: {meta['r1']}")
                    # except me.CntError as ce:
                    #     print(ce.message)
                if self.software == "STAR":
                    cmd = f"{self.soft_path['STAR']} \
                            --genomeDir {self.ref} {self.config} \
                            --readFilesIn {R1} \
                            --outFileNamePrefix {self.output_dir}{meta['run']}. \
                            2>&1 | tee {self.output_dir}{meta['run']}.STAR.log \
                            && {self.soft_path['samtools']} \
                            index {self.output_dir}{meta['run']}.Aligned.sortedByCoord.out.bam"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    pass
            else:
                print("'Layout' Error!")
        return(cmds)

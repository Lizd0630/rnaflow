#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-09-03 00:00:41
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


from src.Meta import Meta
from src.Config import ParseDict
import pathlib as pl
import re
import src.MyError as me


class FastQC:
    def __init__(self, meta_file, input_dir, output_dir, config, softwares):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        self.meta_info = Meta(meta_file).get_info()
        self.input_dir = str(pl.Path(input_dir)) + r"/"
        self.output_dir = str(pl.Path(output_dir)) + r"/"
        self.fastqc = ParseDict(softwares).dict['fastqc']
        self.config = ParseDict(config).make_config()

    def make_cmds(self):
        cmds = []
        allfiles = [str(i) for i in pl.Path(self.input_dir).rglob(r"*")]
        for meta in self.meta_info:
            if meta["layout"] == "PAIRED":
                pattern1 = re.compile(self.input_dir + meta["r1"] + r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[1]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz" + r"$")
                pattern2 = re.compile(self.input_dir + meta["r2"] + r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[2]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz" + r"$")
                R1 = list(filter(lambda x: pattern1.match(x) is not None, allfiles))
                if len(R1) == 1:
                    R1 = R1[0]
                else:
                    raise me.CntError(f"fastq file number error: {meta['r1']}")
                R2 = list(filter(lambda x: pattern2.match(x) is not None, allfiles))
                if len(R2) == 1:
                    R2 = R2[0]
                else:
                    raise me.CntError(f"fastq file number error: {meta['r2']}")
                cmd1 = f"{self.fastqc} {self.config} -o {self.output_dir} {R1} 2>&1 | tee {self.output_dir}{meta['r1']}.fastqc.log"
                cmd2 = f"{self.fastqc} {self.config} -o {self.output_dir} {R2} 2>&1 | tee {self.output_dir}{meta['r2']}.fastqc.log"
                cmds.append(re.sub(" {2,}", " ", cmd1))
                cmds.append(re.sub(" {2,}", " ", cmd2))
            elif meta["layout"] == "SINGLE":
                pattern1 = re.compile(self.input_dir + meta["r1"] + r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[1]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz" + r"$")
                R1 = list(filter(lambda x: pattern1.match(x) is not None, allfiles))
                if len(R1) == 1:
                    R1 = R1[0]
                else:
                    raise me.CntError("fastq file number error: {meta['r1']}")
                cmd = f"{self.fastqc} {self.config} -o {self.output_dir} {R1} 2>&1 | tee {self.output_dir}{meta['r1']}.fastqc.log"
                cmds.append(re.sub(" {2,}", " ", cmd))
            else:
                print("'Layout' Error!")
        return(cmds)

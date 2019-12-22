#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-15 16:39:52
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


from src.Meta import Meta
from src.Config import ParseDict
import pathlib as pl
import re
import src.MyError as me


class Trim:
    def __init__(self, meta_file, input_dir, output_dir, tool, softwares, config):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        self.meta_info = Meta(meta_file).get_info()
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.software = tool
        self.soft_path = ParseDict(softwares).dict
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
                if self.software == "fastp":
                    out1 = f"{meta['r1']}_R1.cln.fq.gz"
                    out2 = f"{meta['r2']}_R2.cln.fq.gz"
                    up1 = f"{meta['r1']}_R1.UP.fq.gz"
                    up2 = f"{meta['r2']}_R2.UP.fq.gz"
                    failed = f"{meta['run']}_failed.fq.gz"
                    j_file = f"{meta['run']}_fastp.json"
                    h_file = f"{meta['run']}_fastp.html"
                    cmd = f"{self.soft_path['fastp']} {self.config} \
                            -i {R1} -I {R2} \
                            -o {self.output_dir}{out1} \
                            --unpaired1 {self.output_dir}{up1} \
                            -O {self.output_dir}{out2} \
                            --unpaired2 {self.output_dir}{up2} \
                            --failed_out {self.output_dir}{failed} \
                            -j {self.output_dir}{j_file} \
                            -h {self.output_dir}{h_file} \
                            2>&1 | tee {self.output_dir}{meta['run']}.fastp.log"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    trimmomatic = self.soft_path['Trimmomatic']
                    trimmomatic_dir = str(pl.Path(trimmomatic).parents[0])
                    out1 = f"{meta['r1']}_R1.cln.fq.gz"
                    out2 = f"{meta['r2']}_R2.cln.fq.gz"
                    up1 = f"{meta['r1']}_R1.UP.fq.gz"
                    up2 = f"{meta['r2']}_R2.UP.fq.gz"
                    cmd = f"java -jar {trimmomatic} PE \
                            {R1} {R2} \
                            {self.output_dir}{out1} {self.output_dir}{up1} \
                            {self.output_dir}{out2} {self.output_dir}{up2} \
                            ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-PE.fa:2:30:10 \
                            {self.config} \
                            2>&1 | tee {self.output_dir}{meta['run']}.trimmomatic.log"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
            elif meta["layout"] == "SINGLE":
                pattern1 = re.compile(self.input_dir + meta["r1"] + r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[1]?[._]?(cln|clean)?[._]?(fq|fastq|fa).gz" + r"$")
                R1 = list(filter(lambda x: pattern1.match(x) is not None, allfiles))
                if len(R1) == 1:
                    R1 = R1[0]
                else:
                    raise me.CntError(f"fastq file number error: {meta['r1']}")
                if self.software == "fastp":
                    out1 = f"{meta['r1']}_cln.fq.gz"
                    failed = f"{meta['r1']}_failed.fq.gz"
                    j_file = f"{meta['r1']}_fastp.json"
                    h_file = f"{meta['r1']}_fastp.html"
                    cmd = f"{self.soft_path['fastp']} {self.config} \
                            -i {R1} \
                            -o {self.output_dir}{out1} \
                            --unpaired1 {self.output_dir}{up1} \
                            --failed_out {self.output_dir}{failed} \
                            -j {self.output_dir}{j_file} \
                            -h {self.output_dir}{h_file} \
                            2>&1 | tee {self.output_dir}{meta['r1']}.fastp.log"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
                else:
                    trimmomatic = self.soft_path['Trimmomatic']
                    trimmomatic_dir = str(pl.Path(trimmomatic).parents[0])
                    out1 = f"{meta['r1']}_cln.fq.gz"
                    cmd = f"java -jar {trimmomatic} SE \
                            {R1} \
                            {self.output_dir}{out1} \
                            ILLUMINACLIP:{trimmomatic_dir}/adapters/TruSeq3-SE.fa:2:30:10 \
                            {self.config} \
                            2>&1 | tee {self.output_dir}{meta['r1']}.trimmomatic.log"
                    cmd = re.sub(" {2,}", " ", cmd)
                    cmds.append(cmd)
            else:
                print("'Layout' Error!")
        return(cmds)

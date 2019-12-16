#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-09-03 10:07:36
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


import pandas as pd
import pathlib as pl
import re
import src.MyError


class Meta:
    def __init__(self, meta_file):
        self.meta_info = pd.read_csv(meta_file, sep="\t")
        self.layouts = ['PAIRED', 'SINGLE']
        self.strands = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']

    def get_info(self):
        run = self.meta_info["Run"]
        r1 = self.meta_info["R1"]
        r2 = self.meta_info["R2"]
        layout = self.meta_info["Layout"]
        strand = self.meta_info["Strand_specificity"]
        if set(layout.unique()).issubset(set(self.layouts)):
            pass
        else:
            raise ValsError("'Layout' value error!")
        if set(strand.unique()).issubset(set(self.strands)):
            pass
        else:
            raise ValsError("'Strand_specificity' value error!")
        return([{"run": run[i],
                "r1": r1[i],
                "r2": r2[i],
                "layout": layout[i],
                "strand": strand[i]} for i in self.meta_info.index])

    def STAR_input(self, input_dir, i):
        '''
        search and match single end RNAseq files.
        '''
        files = []
        run = self.meta_info["run"][i]
        regexpr = re.compile(run + r"[^0-9a-zA-Z]([rR]|[Rr]ead)?[12]?[._]?(fq|fastq|fa).gz")
        all = [str(i) for i in pl.Path(input_dir).rglob(run + r"*")]
        for file in all:
            if not regexpr.search(file):
                files.append(file)
            else:
                pass
        files.sort()
        return files

    def STAR_output(self, output_dir):
        '''
        '''
        STAR_outprefix = output_dir + r"/" + self.meta_info["run"][i] + r"." 
        return STAR_outprefix




    def trim_meta():
        pass

    def align_meta():
        pass

    def bamqc_meta():
        pass

    def cnt_meta():
        pass

    def rsem_meta():
        pass

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2019-09-03 10:07:36
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


import pandas as pd
import src.MyError as me


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
            raise me.ValsError("'Layout' value error!")
        if set(strand.unique()).issubset(set(self.strands)):
            pass
        else:
            raise me.ValsError("'Strand_specificity' value error!")
        return([{"run": run[i],
                 "r1": r1[i],
                 "r2": r2[i],
                 "layout": layout[i],
                 "strand": strand[i]} for i in self.meta_info.index])

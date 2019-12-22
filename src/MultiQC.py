#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-15 18:34:55
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

import pathlib as pl


class MultiQC:
    def __init__(self, input_dir, output_dir, output_name):
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.output_name = output_name

    def make_cmds(self):
        cmds = []
        cmd = f"multiqc -n {self.output_dir}{self.output_name} {self.input_dir}"
        cmds.append(cmd)
        return(cmds)

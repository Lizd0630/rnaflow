#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-22 23:20:45
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


import pathlib as pl


class WriteCmds:
    def __init__(self, cmds, project_name):
        self.cmds = cmds
        self.proj = project_name

    def dump(self, output_dir, prog):
        out_dir = str(pl.Path(output_dir).resolve()) + r"/"
        outfile = out_dir + self.proj + r"." + prog + ".cmds.bash"
        with open(outfile, "w") as out:
            out.write("\n".join(self.cmds) + "\n")

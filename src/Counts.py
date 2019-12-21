#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-17 21:31:29
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-16 01:03:33
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$


from src.Meta import Meta
from src.Config import ParseDict
from src.Config import GA_SCRIPT
import pathlib as pl
import re


class Counts:
    def __init__(self, meta_file, input_dir, output_dir, tool, softwares,
                 config, suffix, ref, ref_type, project_name, n_jobs):
        # ["Run", "R1", "R2", "Layout", "Strand_specificity"]
        # {"run": run[i], "r1": r1[i], "r2": r2[i], "layout": layout[i], "strand": strand[i]}
        if tool == "GenomicAlignments":
            self.meta_info = meta_file
        else:
            self.meta_info = Meta(meta_file).get_info()
            self.config = ParseDict(config).make_config()
        self.input_dir = str(pl.Path(input_dir)) + r"/"
        self.output_dir = str(pl.Path(output_dir)) + r"/"
        self.software = tool
        self.soft_path = ParseDict(softwares).dict
        self.suffix = suffix
        self.ref = ref
        self.ref_type = ref_type
        self.proj = project_name
        self.n_jobs = n_jobs

    def make_cmds(self):
        cmds = []
        if self.software == "GenomicAlignments":
            cmd = f"Rscript {GA_SCRIPT} \
                    -i {self.input_dir} \
                    -o {self.output_dir} \
                    -n {self.proj} \
                    -m {self.meta_info} \
                    -g {self.ref} \
                    -t {self.ref_type} \
                    -p {self.n_jobs} \
                    -s {self.suffix}"
            cmd = re.sub(" {2,}", " ", cmd)
            cmds.append(cmd)
        else:
            pass
        return(cmds)

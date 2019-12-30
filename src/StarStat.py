#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-30 19:58:20
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

import pathlib as pl
import pandas as pd
from functools import reduce
import re

import argparse

parser = argparse.ArgumentParser(description="Merge STAR Log.final.out")
parser.add_argument('-i', '--input_dir', type=str, required=True, help="")
parser.add_argument('-o', '--output_dir', type=str, required=True, help="")
parser.add_argument('-p', '--project_name', type=str, required=True, help="")
args = parser.parse_args()


def read_starLog(file):
    df = pd.read_csv(file, sep="|", header=None)
    filename = re.sub("[^0-9A-Za-z]Log.final.out", "", pl.Path(file).name)
    df.columns = ["Run", filename]
    df = df[df[filename].notna()]
    df[filename] = df[filename].str.replace("\t", "")
    return df


class StarStat:
    def __init__(self, input_dir, output_dir, project_name):
        self.input_dir = str(pl.Path(input_dir).resolve()) + r"/"
        self.output_dir = str(pl.Path(output_dir).resolve()) + r"/"
        self.proj = project_name

    def multi_merge(self):
        allfiles = [str(i) for i in pl.Path("align/").rglob(r"*Log.final.out")]
        alltabs = [read_starLog(file) for file in allfiles]
        df_final = reduce(lambda left, right: pd.merge(left, right, on=['Run']), alltabs)
        df_final["Run"] = df_final["Run"].str.replace(" $", "")
        df_final["Run"] = [re.sub("^ *", "", x) for x in list(df_final["Run"])]
        return df_final.T


def main(input_dir, output_dir, project_name):
    starstat = StarStat(input_dir, output_dir, project_name)
    df = starstat.multi_merge()
    pd.DataFrame.to_csv(df, f'{starstat.output_dir}{starstat.proj}_STARfinal.tsv', sep='\t', na_rep='.', index=True, header=False)


if __name__ == '__main__':
    main(args.input_dir, args.output_dir, args.project_name)

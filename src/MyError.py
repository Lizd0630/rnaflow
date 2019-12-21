#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-12-14 21:59:05
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$

import pathlib as pl


class MyError(Exception):
    """Base class for exceptions in this module."""
    pass


class colsError(MyError):
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message


class ValsError(MyError):
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message


class CntError(MyError):
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message


# meta_info is a data.frame from pandas
def check_cols(self, meta_info):
    cols = ["Run", "R1", "R2", "Layout", "Strand_specificity"]
    if set(cols).issubset(set(self.meta_info.columns.tolist())):
        pass
    else:
        raise colsError(f"{' '.join(cols)} should be in meta cols.")


class existError(MyError):
    def __init__(self, message):
        Exception.__init__(self)
        self.message = message


# check path exist
def check_dir(path):
    path = pl.Path(path)
    if path.exists():
        raise existError(f"{path} doesn't exist.")

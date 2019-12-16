#!/usr/bin/env python3
# -*- coding:utf-8
# @Date    : 2019-12-14 12:31:20
# @Author  : lizd (lzd_hunan@126.com)
# @Link    : ---
# @Version : $Id$
u"""
Call another programs to perform QC process
Created by Zhang yimng at 2018.11.24
"""

import click
from src.FastQC import FastQC
from src.MultiQC import MultiQC
from src.RunCmds import RunCmds
import src.Config as Config
from src.Trim import Trim
from src.Align import Align
from src.Quant import Quant


# ----------------------------------FastQC------------------------------------ #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Run RNAseq QC, FastQC and MultiQC."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input fastq files."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to meta information of samples."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of running jobs at the same time. Finally CPU: n_jobs * threads/job."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    type=str,
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of MultiQC."
)
@click.option(
    "--config",
    default=Config.FASTQC_CONFIG,
    type=click.Path(exists=True),
    help="Path to JSON file(parameters)."
)
@click.option(
    "--silent",
    default=True,
    type=bool,
    help="Whether suppress information in process."
)
@click.version_option(Config.VERSION, message='%(version)s')
def qc(
        input_dir,
        output_dir,
        meta_file,
        n_jobs,
        softwares,
        project_name,
        config,
        silent
):
    u"""
    Performing fastqc of files in meta information.
    """
    soft = Config.Soft(softwares).make_config()
    cmds_fastqc = FastQC(meta_file=meta_file,
                         input_dir=input_dir,
                         output_dir=output_dir,
                         config=config,
                         soft_path=soft['fastqc']).make_cmds()
    cmds_multiqc = MultiQC(input_dir=output_dir,
                           output_dir=output_dir,
                           output_name=project_name).make_cmds()
    RunCmds(cmds=cmds_fastqc, silent=silent).running(n_jobs)
    RunCmds(cmds=cmds_multiqc, silent=silent).running(1)
# ----------------------------------FastQC------------------------------------ #


# ----------------------------------Trimming---------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Run fastq trimming, Trimmomatic or fastp."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input fastq files."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to meta information of samples."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of running jobs at the same time. Finally CPU: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    type=click.Choice(["Trimmomatic", "fastp"]),
    help="Software to use.",
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    type=click.Path(exists=True),
    help="Path to json file of softwares."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of MultiQC if use fastp."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters)."
)
@click.option(
    "--silent",
    default=True,
    type=bool,
    help="Whether suppress information in process."
)
@click.version_option(Config.VERSION, message="%(version)s")
def trim(
        input_dir,
        output_dir,
        meta_file,
        n_jobs,
        tool,
        softwares,
        project_name,
        config,
        silent
):
    u"""
    Performing triming to all input files.
    """
    if config:
        pass
    else:
        if tool == "Trimmomatic":
            config = Config.TRIMMOMATIC_CONFIG
        else:
            config = Config.FASTP_CONFIG
    cmds_trim = Trim(meta_file=meta_file,
                     input_dir=input_dir,
                     output_dir=output_dir,
                     tool=tool,
                     softwares=softwares,
                     config=config).make_cmds()
    cmds_multiqc = MultiQC(input_dir=output_dir,
                           output_dir=output_dir,
                           output_name=project_name).make_cmds()
    RunCmds(cmds=cmds_trim, silent=silent).running(n_jobs)
    if tool == "Trimmomatic":
        pass
    else:
        RunCmds(cmds=cmds_multiqc, silent=silent).running(1)
# ----------------------------------Trimming---------------------------------- #


# ----------------------------------Align------------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform alignment, STAR, GMAP. Need samtools for index."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input fastq files."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to meta information of samples."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of running jobs at the same time. Finally CPU: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="STAR",
    type=click.Choice(["STAR", "GMAP", "GSNAP"]),
    help="Software to use.",
)
@click.option(
    "--ref",
    type=click.Path(exists=True),
    required=True,
    help="Alignment reference, like STAR index."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    type=click.Path(exists=True),
    help="Path to json file of softwares."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of MultiQC if use fastp."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters)."
)
@click.option(
    "--silent",
    default=True,
    type=bool,
    help="Whether suppress information in process."
)
@click.version_option(Config.VERSION, message="%(version)s")
def align(
        input_dir,
        output_dir,
        meta_file,
        n_jobs,
        tool,
        ref,
        softwares,
        project_name,
        config,
        silent
):
    u"""
    Performing triming to all input files.
    """
    if config:
        pass
    else:
        if tool == "STAR":
            config = Config.STAR_CONFIG
            cmds_align = Align(meta_file=meta_file,
                               input_dir=input_dir,
                               output_dir=output_dir,
                               tool=tool,
                               ref=ref,
                               softwares=softwares,
                               config=config).make_cmds()
            cmds_multiqc = MultiQC(input_dir=output_dir,
                                   output_dir=output_dir,
                                   output_name=project_name).make_cmds()
            RunCmds(cmds=cmds_align, silent=silent).running(n_jobs)
            RunCmds(cmds=cmds_multiqc, silent=silent).running(1)
        else:
            print("Not ready!")
# ----------------------------------Align------------------------------------- #


# ----------------------------------Quant------------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform quantification, RSEM."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input bam files."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to meta information of samples."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of running jobs at the same time. Finally CPU: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="RSEM",
    type=click.Choice(["RSEM"]),
    help="Software to use.",
)
@click.option(
    "--ref",
    type=str,
    required=True,
    help="Alignment reference, like RSEM index."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    type=click.Path(exists=True),
    help="Path to json file of softwares."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of MultiQC if use fastp."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters)."
)
@click.option(
    "--suffix",
    type=str,
    default="Aligned.toTranscriptome.out.bam",
    help="Suffix of bam files, default is STAR's."
)
@click.option(
    "--silent",
    default=True,
    type=bool,
    help="Whether suppress information in process."
)
@click.version_option(Config.VERSION, message="%(version)s")
def quant(
        input_dir,
        output_dir,
        meta_file,
        n_jobs,
        tool,
        ref,
        softwares,
        project_name,
        config,
        suffix,
        silent
):
    u"""
    Performing triming to all input files.
    """
    if config:
        pass
    else:
        if tool == "RSEM":
            config = Config.RSEM_CONFIG
            cmds_rsem = Quant(meta_file=meta_file,
                              input_dir=input_dir,
                              output_dir=output_dir,
                              tool=tool,
                              softwares=softwares,
                              config=config,
                              suffix=suffix,
                              ref=ref).make_cmds()
            RunCmds(cmds=cmds_rsem, silent=silent).running(n_jobs)
        else:
            print("Not ready!")
# ----------------------------------Quant------------------------------------- #


# ----------------------------------Count------------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform read counting."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input bam files."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to meta information of samples."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of running jobs at the same time. Finally CPU: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="GenomicAlignmets",
    type=click.Choice(["GenomicAlignmets", "HTSeq"]),
    help="Software to use.",
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    type=click.Path(exists=True),
    help="Path to json file of softwares."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of MultiQC if use fastp."
)
@click.option(
    "--config",
    default=Config.GA_CONFIG,
    type=click.Path(exists=True),
    help="Path to config file(parameters)."
)
@click.option(
    "--suffix",
    type=str,
    default="Aligned.toTranscriptome.out.bam",
    help="Suffix of bam files, default is STAR's."
)
@click.option(
    "--index",
    type=str,
    required=True,
    help="RSEM index."
)
@click.option(
    "--silent",
    default=True,
    type=bool,
    help="Whether suppress information in process."
)
@click.version_option(Config.VERSION, message='%(version)s')
def count(
        input,
        mode,
        output,
        config,
        n_jobs,
        silent
):
    u"""
    Calculate gene, transcripts and exon expression value
    \f
    :param input:
    :param mode:
    :param output:
    :param config:
    :param n_jobs:
    :return:
    """
    click.echo(mode)
# ----------------------------------Count------------------------------------- #


@click.group(
    context_settings=Config.CONTEXT_SETTINGS,
    invoke_without_command=False
)
@click.version_option(Config.VERSION, message='%(version)s')
def main():
    u"""
    Welcome
    \f
    Created by Zhang yiming at 2018.11.14
    """
    pass


if __name__ == '__main__':
    main.add_command(qc)
    main.add_command(trim)
    main.add_command(align)
    main.add_command(quant)
    main.add_command(count)
    main()

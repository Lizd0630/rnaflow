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
import src.Config as Config
from src.WriteCmds import WriteCmds
from src.FastQC import FastQC
from src.MultiQC import MultiQC
from src.RunCmds import RunCmds
from src.Trim import Trim
from src.Align import Align
from src.BamQC import BamQC
from src.Quant import Quant
from src.Counts import Counts
from src.Results import Results


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
    help="Path to directory of input files(fq.gz)."
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
    type=click.Path(exists=True),
    help="Path to file of meta information."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    show_default=True,
    type=str,
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of command-file and MultiQC output."
)
@click.option(
    "--config",
    default=Config.FASTQC_CONFIG,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during processing."
)
@click.version_option(Config.VERSION, message='%(version)s')
def fastqc(
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
    Perform fastqc of samples in meta information.
    """
    cmds_fastqc = FastQC(meta_file=meta_file,
                         input_dir=input_dir,
                         output_dir=output_dir,
                         config=config,
                         softwares=softwares).make_cmds()
    cmds_multiqc = MultiQC(input_dir=output_dir,
                           output_dir=output_dir,
                           output_name=project_name).make_cmds()
    WriteCmds(cmds=cmds_fastqc, 
              project_name=project_name).dump(output_dir=output_dir, prog="fastqc")
    RunCmds(cmds=cmds_fastqc, silent=silent).running(n_jobs)
    RunCmds(cmds=cmds_multiqc, silent=silent).running(1)
# ----------------------------------FastQC------------------------------------ #


# ----------------------------------Trimming---------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Run fastq trimming: Trimmomatic or fastp."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of input files(fq.gz)."
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
    type=click.Path(exists=True),
    help="Path to file of meta information."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="fastp",
    show_default=True,
    type=click.Choice(["Trimmomatic", "fastp"]),
    help="Software to use.",
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of command-file, and MultiQC output if use fastp."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
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
    Perform triming to all files in meta file through specific tools: \n
        fastp \n
        Trimmomatic \n
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
    WriteCmds(cmds=cmds_trim, project_name=project_name).dump(output_dir=output_dir, prog="trim")
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
    type=click.Path(exists=True),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.Path(exists=True),
    help="Path to file of meta information."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="STAR",
    show_default=True,
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
    show_default=True,
    type=click.Path(exists=True),
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of command-file and MultiQC."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
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
    Perform alignments to all sample in meta file through specific tools: \n
        STAR \n
        GMAP \n
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
            WriteCmds(cmds=cmds_align, project_name=project_name).dump(output_dir=output_dir, prog="align")
            RunCmds(cmds=cmds_align, silent=silent).running(n_jobs)
            RunCmds(cmds=cmds_multiqc, silent=silent).running(1)
        else:
            print("Not ready!")
# ----------------------------------Align------------------------------------- #


# ----------------------------------bam QC------------------------------------ #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform bam QC."
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
    type=click.Path(exists=True),
    help="Path to file of meta information."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "--suffix",
    type=str,
    default="Aligned.sortedByCoord.out.bam",
    show_default=True,
    help="Suffix of bam files, default is STAR's."
)
@click.option(
    "--bed12",
    type=click.Path(exists=True),
    required=True,
    help="bed12 file, for RSeQC: infer_experiment, genebody_coverage."
)
@click.option(
    "--refflat",
    type=click.Path(exists=True),
    required=True,
    help="refFlat, for picard: CollectRnaSeqMetrics."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of command-file."
)
@click.option(
    "--config",
    default=Config.BAMQC_CONFIG,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
)
@click.version_option(Config.VERSION, message="%(version)s")
def bamqc(
        input_dir,
        output_dir,
        meta_file,
        n_jobs,
        suffix,
        bed12,
        refflat,
        softwares,
        project_name,
        config,
        silent
):
    u"""
    Perform bam QC: \n
        Picard(CollectRnaSeqMetrics) \n
        RSeQC(infer_experiment, genebody_coverage)
    """
    config = Config.BAMQC_CONFIG
    cmds_infer = BamQC(meta_file=meta_file,
                       input_dir=input_dir,
                       output_dir=output_dir,
                       bed12=bed12,
                       refFlat=refflat,
                       softwares=softwares,
                       config=config,
                       suffix=suffix).infer_expr()
    WriteCmds(cmds=cmds_infer, project_name=project_name).dump(output_dir=output_dir, prog="inferExpr")
    RunCmds(cmds=cmds_infer, silent=silent).running(n_jobs)
    cmds_rnametrics = BamQC(meta_file=meta_file,
                            input_dir=input_dir,
                            output_dir=output_dir,
                            bed12=bed12,
                            refFlat=refflat,
                            softwares=softwares,
                            config=config,
                            suffix=suffix).rnametrics()
    WriteCmds(cmds=cmds_rnametrics, project_name=project_name).dump(output_dir=output_dir, prog="RNAmetrics")
    RunCmds(cmds=cmds_rnametrics, silent=silent).running(n_jobs)
    cmds_genebody = BamQC(meta_file=meta_file,
                          input_dir=input_dir,
                          output_dir=output_dir,
                          bed12=bed12,
                          refFlat=refflat,
                          softwares=softwares,
                          config=config,
                          suffix=suffix).genebody()
    WriteCmds(cmds=cmds_genebody, project_name=project_name).dump(output_dir=output_dir, prog="geneBody")
    RunCmds(cmds=cmds_genebody, silent=silent).running(n_jobs)
# ----------------------------------bam QC------------------------------------ #


# ----------------------------------Quant------------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform quantification."
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
    type=click.Path(exists=True),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.File("r"),
    help="Path to file of meta information."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="RSEM",
    show_default=True,
    type=click.Choice(["RSEM"]),
    help="Software to use.",
)
@click.option(
    "--ref",
    type=str,
    required=True,
    help="Quantification reference, like RSEM index."
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of command-file and results if generate."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task."
)
@click.option(
    "--suffix",
    type=str,
    default="Aligned.toTranscriptome.out.bam",
    show_default=True,
    help="Suffix of input bam files, default is STAR's."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
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
    Perform quantification of samples in meta file though specific tools: \n
        RSEM \n
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
            WriteCmds(cmds=cmds_rsem, project_name=project_name).dump(output_dir=output_dir, prog="RSEM")
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
    type=click.Path(exists=True),
    help="Path to directory of output files."
)
@click.option(
    "-m",
    "--meta_file",
    required=True,
    type=click.Path(exists=True),
    help="Path to file of meta information."
)
@click.option(
    "-c",
    "--counts_type",
    default="gene",
    show_default=True,
    type=click.Choice(["gene", "exon", "exonbin"]),
    help="Counts type."
)
@click.option(
    "-n",
    "--n_jobs",
    default=1,
    show_default=True,
    type=click.IntRange(1, 12, clamp=True),
    help="Number of jobs to run in parallel. Finally CPU usage: n_jobs * threads/job."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="GenomicAlignments",
    show_default=True,
    type=click.Choice(["GenomicAlignments", "HTSeq"]),
    help="Software to use.",
)
@click.option(
    "--softwares",
    default=Config.SOFTWARES,
    show_default=True,
    type=click.Path(exists=True),
    help="Path to JSON file which contains software path."
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of output."
)
@click.option(
    "--config",
    type=click.Path(exists=True),
    help="Path to config file(parameters) of software used in current task. Useless for GenomicAlignmets."
)
@click.option(
    "--suffix",
    type=str,
    default="Aligned.sortedByCoord.out.bam",
    show_default=True,
    help="Suffix of input bam files, default is STAR's."
)
@click.option(
    "--ref",
    type=click.Path(exists=True),
    required=True,
    help="Reference, gtf or gff."
)
@click.option(
    "--ref_type",
    required=True,
    default="gtf",
    show_default=True,
    type=click.Choice(["gtf", "gff"]),
    help="Reference format, 'gtf' or 'gff'."
)
@click.option(
    "--silent",
    default=True,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
)
@click.version_option(Config.VERSION, message='%(version)s')
def count(
        input_dir,
        output_dir,
        meta_file,
        counts_type,
        n_jobs,
        tool,
        softwares,
        project_name,
        config,
        suffix,
        ref,
        ref_type,
        silent
):
    u"""
    Perform reads counting of samples in meta file through specific tools: \n
        GenomicAlignments \n
        HTSeq \n
    """
    if config:
        pass
    else:
        if tool == "GenomicAlignments":
            cmds_ga = Counts(meta_file=meta_file,
                             counts_type=counts_type,
                             input_dir=input_dir,
                             output_dir=output_dir,
                             tool=tool,
                             softwares=softwares,
                             config=config,
                             suffix=suffix,
                             ref=ref,
                             ref_type=ref_type,
                             project_name=project_name,
                             n_jobs=n_jobs).make_cmds()
            RunCmds(cmds=cmds_ga, silent=silent).running(1)
        else:
            print("Not ready!")
# ----------------------------------Count------------------------------------- #


# ------------------------------reuslts process------------------------------- #
@click.command(
    context_settings=Config.CONTEXT_SETTINGS,
    short_help="Perform results files manipulation: STAR final merge, RSEM results merge, etc."
)
@click.option(
    "-i",
    "--input_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of results."
)
@click.option(
    "-o",
    "--output_dir",
    required=True,
    type=click.Path(exists=True),
    help="Path to directory of output files."
)
@click.option(
    "-t",
    "--tool",
    required=True,
    default="STAR",
    show_default=True,
    type=click.Choice(["STAR", "RSEM", "infer_expr", "CollectRnaSeqMetrics", "geneBody_coverage"]),
    help="Software which results produced by.",
)
@click.option(
    "--suffix",
    type=str,
    help="Suffix of individual files of samples.",
)
@click.option(
    "--project_name",
    required=True,
    type=str,
    help="Out prefix of summaried files."
)
@click.option(
    "--silent",
    default=False,
    show_default=True,
    type=bool,
    help="Whether suppress information produced during process."
)
@click.version_option(Config.VERSION, message='%(version)s')
def results(
        input_dir,
        output_dir,
        tool,
        suffix,
        project_name,
        silent
):
    u"""
    Perform results files manipulation: STAR final merge, RSEM results merge, etc. \n
        STAR final logs \n
        RSEM results: TPM/FPKM/reads \n
        RSeQC: infer_experiment \n
        Picard: CollectRNAseqMetrics \n
    """
    if tool == "STAR":
        cmds_star = Results(input_dir=input_dir,
                            output_dir=output_dir,
                            project_name=project_name).star()
        RunCmds(cmds=cmds_star, silent=silent).running(1)
    elif tool == "RSEM":
        cmds_rsem = Results(input_dir=input_dir,
                            output_dir=output_dir,
                            project_name=project_name).rsem()
        RunCmds(cmds=cmds_rsem, silent=silent).running(1)
    elif tool == "infer_expr":
        suffix = Config.RESULTS_SUFFIX["infer_expr"]
        cmds_infer = Results(input_dir=input_dir,
                             output_dir=output_dir,
                             project_name=project_name).infer_expr(suffix)
        RunCmds(cmds=cmds_infer, silent=silent).running(1)
    elif tool == "CollectRnaSeqMetrics":
        suffix = Config.RESULTS_SUFFIX["CollectRnaSeqMetrics"]
        cmds_metrics = Results(input_dir=input_dir,
                               output_dir=output_dir,
                               project_name=project_name).RNAmetrics(suffix)
        RunCmds(cmds=cmds_metrics, silent=silent).running(1)
    elif tool == "geneBody_coverage":
        suffix = Config.RESULTS_SUFFIX["geneBody_coverage"]
        cmds_genebody = Results(input_dir=input_dir,
                                output_dir=output_dir,
                                project_name=project_name).geneBody(suffix)
        RunCmds(cmds=cmds_genebody, silent=silent).running(1)
    else:
        print("Not ready!")
# ------------------------------reuslts process------------------------------- #


@click.group(
    context_settings=Config.CONTEXT_SETTINGS,
    invoke_without_command=False
)
@click.version_option(Config.VERSION, message='%(version)s')
def main():
    u"""
    Welcome to Lizdʻs RNAseq studio! Here we can perform Analysis below: \n
        fastqc: FastQC\n
        trim: Trimming\n
        align: Alignments\n
        bamqc: bam QC\n
        quant: Quantification\n
        count: Reads counting\n
        results: Summarization or Plotting\n
    \n
                                                            Li zhidan, 2019.12
    """
    pass


if __name__ == '__main__':
    main.add_command(fastqc)
    main.add_command(trim)
    main.add_command(align)
    main.add_command(bamqc)
    main.add_command(quant)
    main.add_command(count)
    main.add_command(results)
    main()

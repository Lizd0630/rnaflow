## RNA-seq pipeline
#### lizhidan, lzd_hunan@126.com
rnaflow cover basic steps of RNAseq analyses, but differential analysis do not included.

## For what?
1. fastqc - FastQC: FastQC, multiqc
2. trim - Trimming: fastp, Trimmomatic
3. align - Alignments: STAR
4. bamqc - bam QC: RSeQC(infer_experiment.py, geneBody_coverage.py), picard(collectRnaSeqMetrics)
5. quant - Quantification: RSEM
6. count - Reads counting: GenomicAlignments
7. results - Summarization and plotting

#### Before using this workflow, you should have basic coneption about RNAseq, like: 
- RNA sequencing technologies
- strand specificity in RNAseq library
- Quality control about RNAseq
- Counts and differential analysis
- etc.

## Basic ideas
### Package dependency
- python3: pandas, click, multiqc, RSeQC
- R packages: optparse, GenomicAlignments, GenomicFeatures, ggplot2
- other softwares: FastQC, fastp(0.20.0), Trimmomatic, STAR, picard, RSEM, gtfToGenePred, genePredToBed

### Basic pipe
1. Parse metainformation
2. Parse software path and parameters
3. Combine software, parameters and I/O files to linux commands.
4. Launch commands(parallel)

### File dependency
1. JSON file contains path to softwares be used (dictionary), which should be modified in different environment. [softwares.json](./config/softwares.json)
2. JSON file contains parameters used in software, except for I/O files which will be determinated by metafile (dictionary, leave null if parameter has no values), which should be modified for different version. For STAR2.7.1a, [STAR_align.json](./config/STAR_align.json)
3. Meta information, tsv file(Table seperated file) contain at least 5 columns: 
- Run: prefix of sample name, like alignment output, bamQC, etc. 
- R1: prefix of reads 1 fastq file (raw/clean)
- R2: prefix of reads 2 fastq file, leave NULL if PAIRED (raw/clean)
- Layout: SINGLE/PAIRED
- Strand_specificity: unstranded/fr-firststrand/fr-secondstrand

## Notes
### Software parameter JSON file
#### Input, output, strandness, reference, log file, prefix, suffix, etc, should not in JSON files.
#### Exact performed commands will be recorded in output directory.
#### For details,
0. All input/output files will be parsed from meta file, you just only offer --project_name as prefix or not.
1. Trimmomatic: adapter will be parsed via softpath(Do not support Illumina GA series; Modified Trim Class if you have to)
2. STAR align: --genomeDir be specified via command --ref
3. RSEM quantification: --strandness, --paired-end will be parsed from meta info, RSEM index be specified via command --ref
4. GenomicAlignments count:
```r
flag <- scanBamFlag(isSecondaryAlignment = FALSE,
                    isNotPassingQualityControls = FALSE,
                    isUnmappedQuery = FALSE)
sbp <- ScanBamParam(flag=flag, mapqFilter = 255)
```
##### for single end, unstranded library
```r
summarizeOverlaps(features = ebg,
                  mode = "Union", # exonic reads counting using mode "IntersectionStrict"
                  reads = SE_0_bamlist,
                  ignore.strand = TRUE, ## can be change according to strand-specificity
                  inter.feature = FALSE,
                  singleEnd = TRUE,
                  param = sbp,
                  preprocess.reads = NULL) ## can be change according to strand-specificity
```
##### for paired end, unstranded library
```r
summarizeOverlaps(features = ebg,
                  mode = "Union", # exonic reads counting using mode "IntersectionStrict"
                  reads = PE_0_bamlist,
                  ignore.strand = FALSE,
                  inter.feature = FALSE,
                  singleEnd = FALSE,
                  fragments = FALSE,
                  strandMode = 0, ## can be change according to strand-specificity
                  param = sbp,
                  preprocess.reads = NULL)
```

## Example
### directory tree
```bash
.
|-- align
|-- bam_QC
|-- fastq
|   |-- clean
|   |   `-- QC
|   `-- raw
|       `-- QC
|-- counts
|-- meta
|   `-- meta.tsv
|-- ref
|   |-- chr22.bed12
|   |-- chr22.fa
|   |-- chr22.genePred
|   |-- chr22.gtf
|   |-- chr22.refFlat
|   |-- chr22.sorted.bed12
|   |-- chr22.splicesites
|   |-- chr22.splicesites.fmt
|   |-- RSEM
|   `-- STAR
|-- RSEM
`-- summary
```

### create work space
```bash
mkdir -p align bam_QC fastq/{clean/QC,raw/QC} counts meta ref RSEM summary
```

### RNAseq simulation
```r
setwd("fastq/raw/")
library(polyester)
library(Biostrings)
fasta = readDNAStringSet("RSEM/chr22.transcripts.fa")
readspertx = round(20 * width(fasta) / 100)
fold_changes = matrix(runif(4578 * 3, 1, 5), nrow=4578)
simulate_experiment('RSEM/chr22.transcripts.fa', 
                    reads_per_transcript=readspertx, 
                    num_reps=c(1,1,1), 
                    fold_changes=fold_changes, 
                    outdir='simulated_reads2', 
                    paired=TRUE, 
                    strand_specific=TRUE)
simulate_experiment('RSEM/chr22.transcripts.fa', 
                    reads_per_transcript=readspertx, 
                    num_reps=c(1,1,1), 
                    fold_changes=fold_changes, 
                    outdir='simulated_reads1', 
                    paired=FALSE, 
                    strand_specific=TRUE)
```

### data manipulation
```bash
## fr-secondstrand
seqtk seq -r sample_03.fasta > sample_03r.fasta

## fasta2fastq https://github.com/ekg/fasta-to-fastq.git
perl fasta-to-fastq/fasta_to_fastq.pl fasta > fastq

gzip
```

### examples related software
1. STAR index version: 2.7.1a
2. RSEM index version: v1.3.1
```bash
rsem-prepare-reference --gtf chr22.gtf chr22.fa ./RSEM/chr22
STAR --runMode genomeGenerate \
    --runThreadN 8 \
    --genomeDir ./ \
    --genomeFastaFiles ../chr22.fa \
    --sjdbGTFfile ../chr22.gtf
bash gtf_2_genePred_refFlat_bed12.bash -i chr22.gtf ./
```

### 01 FastQC
```bash
python3 path/to/main.py fastqc -i fastq/raw/ -o fastq/raw/QC -m meta/meta.tsv -n 4 --project_name exam
```
### 02 trim
```bash
python3 path/to/main.py trim -i fastq/raw/ -o fastq/clean/ -m meta/meta.tsv -n 4 --project_name exam -t fastp
```
### 03 align
```bash
python3 path/to/main.py align -i fastq/clean/ -o align/ -m meta/meta.tsv -n 4 --project_name exam -t STAR --ref ref/STAR
## STAR final log
python3 path/to/main.py results -t STAR -i align/ -o summary/ --project_name exam
```
### 04 bamqc
```bash
python3 path/to/main.py bamqc -i align/ -o bam_QC/ -m meta/meta.tsv -n 4 --bed12 ref/chr22.bed12 --refflat ref/chr22.refFlat --project_name exam
## infer experiment
python3 path/to/main.py results -t infer_expr -i bam_QC/ -o summary/ --project_name exam
## RNA metrics
python3 path/to/main.py results -t CollectRnaSeqMetrics -i bam_QC/ -o summary/ --project_name exam
## genebody coverage
python3 path/to/main.py results -t geneBody_coverage -i bam_QC/ -o summary/ --project_name exam
```
### 05 RSEM
```bash
python3 path/to/main.py quant -i align/ -o RSEM/ -m meta/meta.tsv -n 4 --ref ref/RSEM/chr22 --project_name exam
## merge TPM
python3 path/to/main.py results -t RSEM -i RSEM/ -o summary/ --project_name exam
```
### 06 reads counting
```bash
python3 path/to/main.py count -i align/ -o counts/ -m meta/meta.tsv -n 4 --ref ref/chr22.gtf -c gene --project_name exam
```




## To Do
1. remove duplicates and multi-mapped reads
2. Check commands finished or not, in order to not do it again.
3. Error: suffix file exist?





## Learning

### Data manipulation tips
http://bioconductor.jp/packages/3.6/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf



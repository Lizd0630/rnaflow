## RNA-seq pipeline
#### lizhidan

## Basic ideas
### Package dependency
pandas, click

### File dependency
1. JSON file contains path to software be used (dictionary)
2. JSON file contains parameters used in software, except for I/O files which will be determinated by metafile (dictionary, leave null if parameter has no values)
3. Table seperated file contain at least 5 columns: Run, R1, R2, Layout, strand_specificity.


### Basic pipe
1. Parse metainformation
2. Parse software path and parameters
3. Combine software, parameters and I/O files to linux commands.
4. Launch


### example
##### directory tree
```bash
.
|-- align
|-- bam_QC
|-- clean
|   |-- fastp
|   `-- trimmomatic
|-- counts
|-- meta
|   `-- meta.tsv
|-- raw
|   |-- sample_A.fastq.gz
|   |-- sample_B.fastq.gz
|   |-- sample_C.fastq.gz
|   |-- sample_D_1.fastq.gz
|   |-- sample_D_2.fastq.gz
|   |-- sample_E_1.fastq.gz
|   |-- sample_E_2.fastq.gz
|   |-- sample_F_1.fastq.gz
|   `-- sample_F_2.fastq.gz
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

```r
library(polyester)
library(Biostrings)
fasta = readDNAStringSet("RSEM/chr22.transcripts.fa")
readspertx = round(20 * width(fasta) / 100)
fold_changes = matrix(runif(4578 * 3, 1, 5), nrow=4578)
 simulate_experiment('RSEM/chr22.transcripts.fa', reads_per_transcript=readspertx, num_reps=c(1,1,1), fold_changes=fold_changes, outdir='simulated_reads2', paired=TRUE, strand_specific=TRUE)
 simulate_experiment('RSEM/chr22.transcripts.fa', reads_per_transcript=readspertx, num_reps=c(1,1,1), fold_changes=fold_changes, outdir='simulated_reads1', paired=FALSE, strand_specific=TRUE)
```
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
STAR   --runMode genomeGenerate \
       --runThreadN 8 \
       --genomeDir ./ \
       --genomeFastaFiles ../chr22.fa \
       --sjdbGTFfile ../chr22.gtf

bash gtf_2_genePred_refFlat_bed12.bash -i chr22.gtf ./
```

### Notes: software parameter JSON file
#### Input, output, strandness, reference, log file, prefix, suffix, etc, should not in JSON files.
#### For details,
1. fastp: input/output, -j, -h, log file will be parsed from meta info
2. Trimmomatic: adapter will be parsed via softpath
3. STAR align: --genomeDir specified via command, --readFilesIn, --outFileNamePrefix will be parsed from meta info
4. RSEM quantification: --strandness, --paired-end will be parsed from meta info, ref index specified via command


## To Do
1. bam QC(rnaMetrics summary), remove duplicates and multi-mapped reads
2. Check commands finished in order not do it again.
3. Error: parameter with suffix file exist?
4. counts single end strand-specificity

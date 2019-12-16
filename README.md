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


### examples related software
1. STAR index version: 2.7.1a
2. RSEM index version: v1.3.1


### Notes: software parameter JSON file
#### Input, output, strandness, reference, log file, prefix, suffix, etc, should not in JSON files.
#### For details,
1. fastp: input/output, -j, -h, log file will be parsed from meta info
2. Trimmomatic: adapter will be parsed via softpath
3. STAR align: --genomeDir specified via command, --readFilesIn, --outFileNamePrefix will be parsed from meta info
4. RSEM quantification: --strandness, --paired-end will be parsed from meta info, ref index specified via command


## To Do
1. bam QC, raw counts, remove duplicates and multi-mapped reads
2. Check commands finished in order not do it again.
3. Error

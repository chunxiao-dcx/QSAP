# QSAP
QSAP version 1.0  
##Introduction

QSAP is an automatic annotation pipeline for fast annotation and classification of QS-like sequences from sequencing data.  
Author: DAI Chunxiao.  
Email: 2446378845\@qq.com.  
## Prepare compulsory command  
1. Diamond Buchfink B, Reuter K, Drost HG, \"Sensitive protein alignments at tree-of-life
    scale using DIAMOND\", *Nature Methods* **18**, 366–368 (2021).
    [doi:10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)
2. HUMMR 
## Usage  
```
$ perl QSAP.pl -i \<Input files list\> -o \<Output dir\> -s \[sub|union|hmmscan|diamond\] -r \<raw reads list\> -n \<2.0\> -p \[conda|/your/install/path/\] -d \<0.50\> -G -k \<21\> -A \<Gene abundance table list\> -h
```

**General options:**  
 --input(-i)      Input files list,necessary.  
 --output(-o)     Output files directory, default current directory.  
 --strategy(-s)   Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.  
 --rawread(-r)    The metagenome raw reads list for salmon quant, if -r, Input files should be nucleotide sequences.  
 --thread(-n)     Number of threads, default value is 2.  
 --path(-p)       Software install directory: default installed by conda, or set to /your/install/path/. EXAMPLE: [-p /home/soft_for_qsg/bin/].  

**Diamond parameters:**  
 --id(-d)         The identity value, default 50.  
 --diamondopt  Other options in Diamond blastp, the "-" need to be "/-". Inalterable options: --in;--db;--out;--outfmt;--max-target-seqs. EXAMPLE: [--diaoptions /--query-cover 50 /--fast /-e 0.00001].  

**Hmmscan parameters:**  
 -G               Default option for model-specific thresholding. if -G , use profile's GA gathering cutoffs to set all thresholding.
 --hmmscanopt     Other options in Hmmscan, the "-" need to be "/-". if --hmmscanopt, -G will not be allowed. Inalterable options: -o;--tblout;--noali. Alterable options: -E;-T;--domE;--domT. EXAMPLE: [--hmmscanopt /-E 0.01 /--domE 0.01];[--hmmscanopt /-T 30 /--domT 30].

**Salmon parameters:**  
 --Kmers(-k)      Kmers for building salmon index, default value is 31 .
 --salmonopt      Other options.in salmon quant, the "-" need to be "/-". Inalterable options: --validateMappings;-i;-l;--meta;-1;-2;-o. EXAMPLE: [--saloptions "/--hardFilter /--recoverOrphans"].

**Others:**  
 -A               Gene abundance table list. if -A, --rawread will not be allowed. The target columns should be named as "Name" and "Abundance".
 -h               Print help information.
 
 
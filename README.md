# QSAP
QSAP version 1.0  
Author: DAI Chunxiao.  
Email: 2446378845\@qq.com.

# Introduction
Quorum sensing (QS) is thought as an important process in bacterial cell-to-cell communication which coordinated population behaviors in the bacterial world. **QSAP is an automatic annotation pipeline for fast annotation and classification of QS-related sequences from sequencing data.**  


# Usage  

**` perl QSAP.pl -i <Input files list> -o <Output dir> -s [sub|union|hmmscan|diamond] -r <raw reads list> -n <2.0> -p [/your/install/path/] -d <0.50> -G -k <21> -A <Gene abundance table list> -h`**
   
**General options:**  
&emsp;`--input(-i)` &emsp;&emsp;Input files list,necessary.  
&emsp;`--output(-o)` &emsp;&ensp;Output files directory, default current directory.  
&emsp;`--strategy(-s)`&ensp; Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.  
&emsp;`--rawread(-r)`&ensp;&ensp; The metagenome raw reads list for salmon quant, if `-r`, Input files should be nucleotide sequences.   
&emsp;`--thread(-n)`&emsp;&ensp; Number of threads, default value is 2.  
&emsp;`--path(-p)`&emsp;&ensp;&emsp; Software install directory: Default installed by conda, or set to /your/install/path/.  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;EXAMPLE: `-p /home/soft_for_qsap/bin/`.  

**Diamond parameters:**  
 &emsp;`--id(-d)`&emsp;&emsp;&emsp;&ensp;The identity value, default 50.  
 &emsp;`--diamondopt`  &ensp;&ensp; Other options in Diamond blastp, the "-" need to be "/-".   
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Inalterable options: `--in`; `--db`; `--out`; `--outfmt`; `--max-target-seqs`.   
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;EXAMPLE: `--diaoptions /--query-cover 50 /--fast /-e 0.00001`.  

**Hmmscan parameters:**  
 &emsp; `-G` &emsp;&emsp;&emsp;&emsp;&ensp;&ensp;&ensp; Default option for model-specific thresholding. if `-G`, use profile's GA gathering cutoffs to set all thresholding.  
 &emsp;`--hmmscanopt`&emsp;&ensp;Other options in Hmmscan, the "-" need to be "/-". if `--hmmscanopt`, `-G` will not be allowed.  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Inalterable options: `-o`; `--tblout`; `--noali`.  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Alterable options: `-E`; `-T`; `--domE`; `--domT`.  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;EXAMPLE: `--hmmscanopt /-E 0.01 /--domE 0.01`;`--hmmscanopt /-T 30 /--domT 30`.

**Salmon parameters:**  
  &emsp;`--Kmers(-k)` &emsp;&ensp;Kmers for building salmon index, default value is 31.  
  &emsp;`--salmonopt`&emsp;&ensp; Other options.in salmon quant, the "-" need to be "/-".  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Inalterable options: `--validateMappings`;`-i`;`-l`;`--meta`;`-1`;`-2`;`-o`.  
 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;EXAMPLE: `--saloptions /--hardFilter /--recoverOrphans`.

**Others:**  
 &emsp;`-A` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Gene abundance table list. if `-A`, `--rawread` will not be allowed. The target columns should be named as "Name" and "Abundance".  
 &emsp;`-h`&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;Print help information.
 

# Before start
## Software involved in this pipeline  
1. **Diamond v2.0.15.153**: Buchfink B, Reuter K, Drost HG, \"Sensitive protein alignments at tree-of-life scale using DIAMOND\", *Nature Methods* 18, 366–368 (2021). [DOI: 10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)&emsp;[Git-hub: https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond).  
2. HMMER 3.3.2: [DOI:](http://hmmer.org/)&emsp;[Git-hub:https://github.com/EddyRivasLab/hmmer](https://github.com/EddyRivasLab/hmmer).  
3. Seqkit 2.2.0: [DOI:](http://hmmer.org/)&emsp;[Git-hub:](http://hmmer.org/).
### Installation: 
default: installed by conda.  
if `-p` \[your/install/path/\].  
### notice:

## Input files list

## Metagenome raw reads list

## Gene abundance table list



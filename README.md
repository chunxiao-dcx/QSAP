# QSAP

# Introduction
**QSAP version 1.0**    
&ensp;Quorum sensing (QS) is thought as an important process in bacterial cell-to-cell communication which coordinated population behaviors in the bacterial world. **QSAP** is an automatic annotation pipeline for fast annotation and classification of **QS-related sequences** from sequencing data. This pipeline is constructed based on 809721 protein sequences from the [QSP database](https://github.com/chunxiao-dcx/QSP).  
Author: DAI Chunxiao.  
Email: 2446378845\@qq.com.  

# Workflow
![image](https://github.com/chunxiao-dcx/QSAP/blob/main/QSAPpipeline.png)

# Before start
To run QSAP, users should download the QSAP source code into local computer system (Unix/Linux), and installed the software involved in QSAP.
### clone source code  
`git clone https://github.com/chunxiao-dcx/QSAP.git`  

## involved Software 
### Software:   
1. **Diamond v2.0.15.153**: Buchfink B, Reuter K, Drost HG, Sensitive protein alignments at tree-of-life scale using DIAMOND, *Nature Methods*. [DOI: 10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)&emsp;[Git-hub: https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond).  
2. **HMMER 3.3.2**: Potter SC, Luciani A, Eddy SR, Park Y, 
Lopez R, Finn RD,HMMER web server: 2018 update,*Nucleic Acids Res*.[DOI: 10.1093/nar/gky448](http://doi.org/10.1093/nar/gky448)&emsp;[Git-hub:https://github.com/EddyRivasLab/hmmer](https://github.com/EddyRivasLab/hmmer).  
3. **SeqKit 2.2.0**: 
W Shen, S Le, Y Li\*, F Hu\*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. *PLOS ONE*. [DOI: 10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)&emsp;[Git-hub:https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit).
### notice:
The copyrights of the involved software belong to their original developers.  
### Installation:  
**Default:**
     Installed by conda, and the script will Using these software directly under conda environment.The user of this pipeline should follow the guideline and regulations of these software.  
**if `-p \[your/install/path/\]`:**  

##  Laguage used in QSAP
QSAP is writed using Perl language, and used the following Module from [CPAN (Comprehensive Perl Archive Network)](www. cpan.org).
1. 
## Prepare your input files list [necessary]
Input file list is nessasery for QSAP. The list contained 

File|path| Character
---------|-----------------------|------
TM1.fasta|~/QSAP/example/test/nuc| nuc
TM2.fasta|~/QSAP/example/test/nuc| nuc    
TP.fasta|~/QSAP/example/test/pro| pro

### Tips:

## Prepare Metagenome raw reads list for gene abundance caculatation
Metagenome raw reads list

File name |Forward file|Forward file|Path
---|-----------|-----------|-----------------------
TM1|TM1_1.fq.gz|TM1_2.fq.gz|~/QSAP/example/testpair
TM2|TM2_1.fq.gz|TM2_2.fq.gz|~/QSAP/example/testpair  

## Prepare your own gene abundance table list 

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
 &emsp;`-A` &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;Gene abundance table list. if `-A`, `--rawread` will not be allowed.  
&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;The target columns should be named as "Name" and "Abundance".  
 &emsp;`-h`&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&ensp;Print help information.

## Note
&ensp;This pipeline is distributed in the hope to achieve the aim of management of QS-related sequences in envrionment, but WITHOUT ANY WARRANTY. No other warranty is expressed or implied including warranties or merchantability and fitness for any particular purpose. This pipeline is only allowed to be used for non-commercial and academic purpose.

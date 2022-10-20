# QSAP

# Introduction
**QSAP version 1.0**    
&emsp;Quorum sensing (QS) is thought as an important process in bacterial cell-to-cell communication which coordinated population behaviors in the bacterial world. **QSAP** is an automatic annotation pipeline for fast annotation and classification of **QS-related sequences** from sequencing data. This pipeline is constructed based on 809721 protein sequences from the [QSP database](https://github.com/chunxiao-dcx/QSP). This pipeline provided two homolog alignment strategies offered by Diamond (Blastp) or HMMER (Hmmscan), as well as a data cleansing function for a subset or union set of the hits. The establishment of this pipeline, which provided a useful tool for QS-related sequence annotation in a wide range of projects, will increase our understanding of QS communication in the bacteria world.  
&emsp;Author: DAI Chunxiao.  
&emsp;Email: 2446378845\@qq.com.  

# Workflow
![image](https://github.com/chunxiao-dcx/QSAP/blob/main/QSAPpipeline.png)

# Before start
## Clone source code  
&emsp;To run QSAP, users should download the QSAP source code into local computer system (Unix/Linux), and installed the software involved in QSAP.  
&emsp;**`git clone https://github.com/chunxiao-dcx/QSAP.git`** 

## Involved Software 
1. **Diamond v2.0.15**: 
Buchfink B, Reuter K, Drost HG, Sensitive protein alignments at tree-of-life scale using DIAMOND, *Nature Methods*. [DOI: 10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x)&emsp;[Git-hub: https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond).  
2. **HMMER 3.3.2**: 
Potter SC, Luciani A, Eddy SR, Park Y, Lopez R, Finn RD,HMMER web server: 2018 update,*Nucleic Acids Res*.[DOI: 10.1093/nar/gky448](http://doi.org/10.1093/nar/gky448)&emsp;[Git-hub:https://github.com/EddyRivasLab/hmmer](https://github.com/EddyRivasLab/hmmer).  
3. **SeqKit 2.2.0**: 
Shen W, Le S, Li Y, Hu F, SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation, *PLOS ONE*. [DOI: 10.1371/journal.pone.0163962](https://doi.org/10.1371/journal.pone.0163962)&emsp;[Git-hub:https://github.com/shenwei356/seqkit](https://github.com/shenwei356/seqkit).  
4. **Salmon 1.8.0**:
Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C, Salmon provides fast and bias-aware quantification of transcript expression, *Nature Methods*. [DOI: 10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)&emsp;[Git-hub:https://github.com/COMBINE-lab/salmon](https://github.com/COMBINE-lab/salmon)

### notice:
&emsp;The copyrights of the involved software belong to their original developers, please cite orrectly if they are integral to your study.
  
### Installation:  
**1. Default:**  
&emsp;Install by conda. and the script will Using these software directly under conda environment.The user of this pipeline should follow the guideline and regulations of these software.  
    
    conda create -n software  
    conda activate software
    conda install diamond=2.0.15
    conda install HMMER=3.3.2
    conda install -c bioconda seqkit=2.3.1
    conda install Salmon=1.8.0

**2. If `-p \[your/install/path/\]`:**  
&emsp;To install in the specified folder: 

    mkdir -p /your/install/path
    cd /your/install/path
    # Install Diamond
    wget https://github.com/bbuchfink/diamond/releases/download/v2.0.15/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
    # Install HMMER
    wget http://eddylab.org/software/hmmer/hmmer.tar.gz
    tar zxf hmmer.tar.gz
    cd hmmer-3.3.2 
    ./configure --prefix /your/install/path/hmmer
    make
    make install
    # Install SeqKit
    cd /your/install/path
    wget https://github.com/shenwei356/seqkit/releases/download/v2.3.1/seqkit_linux_amd64.tar.gz
    tar zxvf seqkit_linux_amd64.tar.gz
    # Install Salmon
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz
    tar zxvf salmon-1.8.0_linux_x86_64.tar.gz

##  Language & Modules
&emsp;QSAP is writed using Perl language, and used the following Modules from [CPAN \(Comprehensive Perl Archive Network\)](www.cpan.org). Before start, the Perl and modules need to be installed correctly following the guideline.  
1. [Getopt::Long](https://metacpan.org/pod/Getopt::Long)  
2. [FindBin](https://metacpan.org/pod/FindBin)  
3. [Data::Dumper](https://metacpan.org/pod/Data::Dumper)  
4. [Array::Utils](https://metacpan.org/pod/Array::Utils)  
5. [List::Util](https://metacpan.org/pod/List::Util) 

## Prepare your input files list [Necessary]
&emsp;A input file list is nessasery for QSAP which contain the name of orignal sequences file, its pathway and its character. The files need to be FASTA/Q format containing nucleotide sequences ("nuc") or protein sequences("pro"). The input file list looks like this:

File|Pathway| nuc/pro
---------|-----------------------|------
TM1.fasta|~/QSAP/example/test/nuc| nuc                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
TM2.fasta|~/QSAP/example/test/nuc| nuc    
TP.fasta|~/QSAP/example/test/pro| pro

### Tips:
&emsp;The pepiline will report error when exist extra points with no meaning of format in file names. For example: "TM1.fasta" is ok, "TM1.1.fasta" is wrong.

## Prepare pair-end reads list for gene abundance estimation [Optional]
&emsp;QSAP also provides a method of gene abundance estimation using [Salmon](https://github.com/COMBINE-lab/salmon). if `-r`, nucleotide sequences and pair-end reads are needed. And the pair-end reads list looks like this: 

File name |Left reads|Right reads|Path
---|-----------|-----------|-----------------------
TM1|TM1_1.fq.gz|TM1_2.fq.gz|~/QSAP/example/testpair
TM2|TM2_1.fq.gz|TM2_2.fq.gz|~/QSAP/example/testpair  
### Tips:
"File name" in this table need to be corresponded to the files in the input files list
for example: 
Column name|content 
---------------------------|--------
File \<input files list\>|TM1.fasta
File name \<pair-end reads list\>|TM1

## Prepare your own gene abundance tables list [Optional]
You can also provide own gene abundance tables for merge the gene abundance . The gene abundance tables list looks like this:
File name|Gene abundance table
----|-----------------------
YM1|~/QSAP/example/GeneAbundance/YM1.sf
YM2|~/QSAP/example/GeneAbundance/YM2.sf

### Tips:
Same as the pair-end reads list, file names in gene abundance tables list need to be corresponded to the files in the input files list. Gene abundance table for each sample need contain the column "Name" and "abundance", and QSAP will automatic merge the gene abunance and the alignment result. For example:

Name|Length|Abundance
-----------|---|------
k141_8534_1|330|66.23
k141_502_1|123|18.22
k141_1506_1|312|27.77

# Usage  
**` perl QSAP.pl -i <Input files list> -o <Output dir> -s [sub|union|hmmscan|diamond] -r <Pair-end reads list> -n <2.0> -p [/your/install/path] -d <0.50> -G -k <21> -A <Gene abundance tables> -h`**
   
**General options:**  
&emsp;`--input(-i)` &emsp;&emsp;Input files list,necessary.  
&emsp;`--output(-o)` &emsp;&ensp;Output files directory, default current directory.  
&emsp;`--strategy(-s)`&ensp; Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.  
&emsp;`--rawread(-r)`&ensp;&ensp; The metagenome raw reads list for salmon quant, if `-r`, Input files should be nucleotide sequences.   
&emsp;`--thread(-n)`&emsp;&ensp; Number of threads, default value is 2.  
&emsp;`--path(-p)`&emsp;&ensp;&emsp; Software install directory: Default installed by conda, or set to /your/install/path.  
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



# Note
&emsp;This pipeline is distributed in the hope to achieve the aim of management of QS-related sequences in envrionment, but WITHOUT ANY WARRANTY. No other warranty is expressed or implied including warranties or merchantability and fitness for any particular purpose. This pipeline is only allowed to be used for non-commercial and academic purpose.

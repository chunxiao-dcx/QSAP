# QSAP
QSAP_version_1.0  
An automatic annotation pipeline for fast annotation and classification of QS-like sequences from sequencing data.

Author: DAI Chunxiao.

Email: 2446378845\@qq.com.

# Options
```sh
$ perl QSAP.pl -i <Input files list> -o <Output dir> -s [sub|union|hmmscan|diamond] -r <raw reads list> -n <2.0> -p [conda|/your/install/path/] -d <0.50> -G -k <21> -A <Gene abundance table list> -h
```

General options:  
 --input(-i)	Input files list,necessary.  
 --output(-o)	Output files directory, default current directory.  
 --strategy(-s)	Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.  
 --rawread(-r) The metagenome raw reads list for salmon quant, if -r, Input files should be nucleotide sequences.  
	--thread(-n)	Number of threads, default value is 2.  
	--path(-p)	Software install directory: default installed by conda, or set to /your/install/path/.  
				EXAMPLE: [-p /home/soft_for_qsg/bin/].  

Diamond parameters:  
	--id(-d)	The identity value, default 50.  
	--diamondopt	Other options in Diamond blastp, the "-" need to be "/-".  
				Inalterable options: --in;--db;--out;--outfmt;--max-target-seqs.  
				EXAMPLE: [--diaoptions /--query-cover 50 /--fast /-e 0.00001].  

Hmmscan parameters:
	-G		default option for model-specific thresholding. if -G , use profile's GA gathering cutoffs to set all thresholding.
	--hmmscanopt	Other options in Hmmscan, the "-" need to be "/-". if --hmmscanopt, -G will not be allowed.
				Inalterable options: -o;--tblout;--noali.
				Alterable options: -E;-T;--domE;--domT.
				EXAMPLE: [--hmmscanopt /-E 0.01 /--domE 0.01];[--hmmscanopt /-T 30 /--domT 30].

Salmon parameters:
	--Kmers(-k)	Kmers for building salmon index, default value is 31.
	--salmonopt	other options.in salmon quant, the "-" need to be "/-".
				\033[34mInalterable options:\033[0m --validateMappings;-i;-l;--meta;-1;-2;-o.
				\033[34mEXAMPLE:\033[0m [--saloptions "/--hardFilter /--recoverOrphans"].

Others:
	-A		Gene abundance table list. if -A, --rawread will not be allowed. The target columns should be named as "Name" and "Abundance".
	-h		Print help information.
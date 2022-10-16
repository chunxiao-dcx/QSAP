# QSAP
An automatic annotation pipeline, named QSAP, was developed for fast annotation and classification of QS-like sequences from sequencing data.
QSAP_version_1.0
This program is used for quorum sensing related genes prediction.
Author: DAI Chunxiao.
Email: 2446378845\@qq.com

$0 -i <Input files list> -o <Output dir> -s [sub|union|hmmscan|diamond] -r <raw reads list> -n <2.0> -p [conda|/your/install/path/] -d <0.50> -G -k <21> -A <Gene abundance table list> -h

General options:
	--input(-i)	Input files list,necessary.
	--output(-o)	Output files directory, default current directory.
	--strategy(-s)	Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.
	--rawread(-r) The metagenome raw reads list for salmon quant, if -r, Input files should be nucleotide sequences.
	--thread(-n)	Number of threads, default value is 2.
	--path(-p)	Software install directory: default installed by conda, or set to /your/install/path/.
				\033[34mEXAMPLE:\033[0m [-p /home/soft_for_qsg/bin/].

Diamond parameters:
	--id(-d)	The identity value, default 50.
	--diamondopt	Other options in Diamond blastp, the "-" need to be "/-".
				\033[34mInalterable options:\033[0m --in;--db;--out;--outfmt;--max-target-seqs.
				\033[34mEXAMPLE:\033[0m [--diaoptions /--query-cover 50 /--fast /-e 0.00001].

Hmmscan parameters:
	-G		default option for model-specific thresholding. if -G , use profile's GA gathering cutoffs to set all thresholding.
	--hmmscanopt	Other options in Hmmscan, the "-" need to be "/-". if --hmmscanopt, -G will not be allowed.
				\033[34mInalterable options:\033[0m -o;--tblout;--noali.
				\033[34mAlterable options:\033[0m -E;-T;--domE;--domT.
				\033[34mEXAMPLE:\033[0m [--hmmscanopt /-E 0.01 /--domE 0.01];[--hmmscanopt /-T 30 /--domT 30].

Salmon parameters:
	--Kmers(-k)	Kmers for building salmon index, default value is 31.
	--salmonopt	other options.in salmon quant, the "-" need to be "/-".
				\033[34mInalterable options:\033[0m --validateMappings;-i;-l;--meta;-1;-2;-o.
				\033[34mEXAMPLE:\033[0m [--saloptions "/--hardFilter /--recoverOrphans"].

Others:
	-A		Gene abundance table list. if -A, --rawread will not be allowed. The target columns should be named as "Name" and "Abundance".
	-h		Print help information.
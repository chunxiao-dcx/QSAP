#!/usr/bin/perl -w
use strict;

#QSG_1.0_version  by dcx
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#Usage information
my $Usage = <<USAGE;
########################
QSAP_1.0_version
This program is used for quorum sensing related genes prediction.
Author: DAI Chunxiao.
Email: 2446378845\@qq.com

$0 -i <Input files list> -o <Output dir> -s [sub|union|hmmscan|diamond] -r <raw reads list> -n <2.0> -p [conda|/your/install/path/] -d <0.50> -G -k <21> -A <Gene abundance table list> -h

General options:
	--input(-i)	Input files list,necessary.
	--output(-o)	Output files directory, default current directory.
	--strategy(-s)	Softwares used for extrated QSGs, default use both hmmscan and diamond blastp, and keep the subset.
	--rawread(-r)	The metagenome pair-end reads list for salmon quant, if -r, Input files should be nucleotide sequences.
	--thread(-n)	Number of threads, default value is 2.
	--path(-p)	Software install directory: default installed by conda, or set to /your/install/path/. The final "/" is nessasery.
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
	-A		Gene abundance tables list. if -A, --rawread will not be allowed. The target columns should be named as "Name" and "Abundance".
	-h		Print help information.
USAGE

##########################################################################
#1 Neccessary input information and read the file list and generalize dir for this program;
my $date = localtime;
my $root_path = $Bin;
print "########################\nLog on $date.\n########################\n\033[37;1m1 Read the file list.\033[0m\n";
$date =~ s/\s/\_/g;
my ($help,$GA)=(0,1);
my ($hmmE,$hmmT,$hmmdomE,$hmmdomT)=(20,1e-5,20,1e-5);
my $thread = 2;
my $identity = 0.50;
my $input= "none";
my $rawread = "none";
my $output = "$root_path/QSG_output_$date";
my $strategy = "sub";
my $bin_path = "$root_path/bin";
my $soft_path = " ";
my $Kmers = 31;
my $Abundance = "none";
my (@Diaoptions,@Hmmoptions,@Saloption);
my %options;
GetOptions(
	'help|h' => \$help,
	'input|i=s' => \$input,
	'output|o=s' => \$output,
	'strategy|s=s' => \$strategy,
	'rawread|r=s' => \$rawread,
	'thread|n=i' => \$thread,
	'path|p=s' => \$soft_path,
	'id|d=o' => \$identity,
	'G' => \$GA,
	'A=s' => \$Abundance,
	'Kmers|k=i' => \$Kmers,
	'diamondopt=s{,}' => \@Diaoptions,
	'hmmscanopt=s{,}' => \@Hmmoptions,
	'saloption=s{,}' => \@Saloption,
);
#print "$GA*******\n";
print "File List: $input.\n";
if(($input eq "none")||($help==1)){
	die "$Usage\n";
}else{
	open(IN,"<$input")|| die "\033[31;1mERROR:\033[0m Can not open the files list: $!\n";
	readline IN;
}
unless($strategy eq "sub" || ($strategy eq "hmmscan") || ($strategy eq "diamond")){
	die "\033[31;1mERROR:\033[0m"."-u [sub|union|hmmscan|diamond] Strategy used for aligned QSGs, default keep the subset.\n";
}
if((($GA==1)&&($Hmmoptions[0] eq ""))||(($GA==1)&&(@Hmmoptions))){
	die "\033[31;1mERROR:\033[0m"."if --hmmscanopt, -G will not be allowed.";
}
print "@Hmmoptions******\n";
unless(-d $output){
	`mkdir $output`;
}
print "The results will be restored in $output.\n";

#read the raw reads list;
our %list;

print "The abundance list:\n";
if(($rawread eq "none")&&($Abundance eq "none")){
	print "No metagenome raw reads list or Gene abundance table list was detected.\n";
}elsif((!($rawread eq "none"))&& (!($Abundance eq "none"))){
	die "\033[31;1mERROR:\033[0m if --rawread|-r, -A will not be allowed.";
}elsif(!($Abundance eq "none")){
	(open IN_AB,"<$Abundance")||die "\033[31;1mERROR:\033[0m Can not open Gene abundance table list: $!\n";
	readline IN_AB;
	while(<IN_AB>){
		chomp;
		my($name,$other) = (split /\t/,$_,2)[0,1];
		print "$name\t$other\n";
		open INSS,"<$other";
		our $firstline = readline INSS;
		my @st=(split /\t/,$firstline);
		unless((grep /^Name$/, @st)&&(grep /^Abundance$/, @st)){
			die "\033[31;1mERROR:\033[0m No target columns named as \"Name\" and \"Abundance\" in $other.";
		}
		close INSS;
	}
	close IN_AB;
}elsif(!($rawread eq "none")){
	open (RAW,"<$rawread") || die "\033[31;1mERROR:\033[0m Can not open raw reads list: $!\n\n";
	readline RAW;
	while(<RAW>){
		chomp;
		$_=~ s/\r//g;
		my($name,$other) = (split /\t/,$_,2)[0,1];
		$list{$name}=$other;
		print "$name\t$list{$name}\n";
		}
	close RAW;
}

#Output dir
my $database_path = "$root_path/DB";
my $output_sh = "$output/QSG_test_version_$date.sh";
my $output_Log = "$output/QSG_test_version_$date.log";
open (SH_file, ">>$output_sh"); 
open (STDOUT, "| tee -ai $output_Log");

##########################################################################
#2 Make the sh file;
print "########################\n\033[37;1m2 Make the sh file.\033[0m\n";
print SH_file "#!/bin/bash\n";
our @filess;
while(<IN>){
	chomp;
	my($name,$pathway,$character) = (split /\t/,$_,3)[0,1,2];
	print "$name | ";
	print SH_file "echo -e \"\\033[33mDealing with the $name file\.\\033[0m\"\n";
	my ($file_name_noch,$fa,$gz) = (split /\./,$name,3)[0,1,2];
	push @filess, $file_name_noch;
	$character =~ s/\r//g;
	my $pathway_ori = "";
	my $name_ori = "";
	#pro or nuc
	if($character eq "pro"){
		#print "";
	}elsif($character eq "nuc"){
		my ($RESULT_trans) = &trans_to_nuc($name,$pathway);
		print SH_file "$RESULT_trans\n";
		$pathway_ori = "$pathway";
		$name_ori= "$name";
		$pathway = "$output/Trans_file/$file_name_noch";
		$name = "$file_name_noch.pro.fasta";
	}else{
		print "\033[31;1mERROR:\033[0m"."The character of $name file should be 'pro' or 'nuc'.\n";
		die "\033[31;1mERROR\033[0m";
	}
	
	my $RESULT_hmmscan = &predict_QSG_hmmscan($name,$pathway);
	my $RESULT_blastp = &predict_QSG_diamond($name,$pathway);
	my $RESULT_salmon;
	if ($strategy eq "sub"||$strategy eq "union"){
		print SH_file "$RESULT_hmmscan\n";
		print SH_file "$RESULT_blastp\n";
	}elsif($strategy eq "hmmscan"){
		print SH_file "$RESULT_hmmscan\n";
	}elsif($strategy eq "diamond"){
		print SH_file "$RESULT_blastp\n";
	}else{
		print "\033[31;1mERROR:\033[0m"."-u [sub|union|hmmscan|diamond] Softwares used for extrated QSGs, default use both hmmscan and diamond blastp.\n";
		die;
	}
	if ((exists $list{$file_name_noch})&&$character eq "nuc"){
		$RESULT_salmon = &salmon_abundance($name_ori,$pathway_ori,$list{$file_name_noch});
		print SH_file "$RESULT_salmon";
	}else{
		#print "";
	}
}

print "\nOk\n";
close IN;

##transform to protein
sub trans_to_nuc{
	my $file_name = shift @_;
	my $file_path = shift @_;
	my $protein = ();
	my ($file_name_noch,$fa,$gz) = (split /\./,$file_name,3)[0,1,2];
	my $format_path = "$output\/Trans_file";
	my $file_pro = "$format_path/$file_name_noch/$file_name_noch";
	$gz ||="none";
	unless(-d $format_path){
		`mkdir $format_path`;
	}
	
	if($fa eq "fasta"||$fa eq "fa"||$fa eq "ffn"||$fa eq "faa"||$fa eq "fna"){
		$protein = "$soft_path"."seqkit translate $file_path/$file_name -o $file_pro.pro.fasta";
	}elsif($fa eq "fastq"||$fa eq "fq"){
	my $trans_fq_fa= "$format_path/$file_name_noch/$file_name_noch.fasta";
	$protein = "$soft_path"."seqkit translate $trans_fq_fa -o $file_pro.pro.fasta\n";
	}else{
		print "\033[31;1mERROR:\033[0m"."The file formate should be \".fa/fasta/fq/fastq\", \".gz\" format files will be identified automatic";
		die "\033[31;1mERROR: $date\033[0m";
	}
	return $protein;
}

##hmmscan
sub predict_QSG_hmmscan{
	my $file_name = shift @_;
	my $file_path = shift @_;
	my $RESULT_hmmdir = "$output/Hmmscanresults";
	unless(-d $RESULT_hmmdir || $strategy eq "diamond"){
		`mkdir $RESULT_hmmdir`;
	}
	my ($file_name_noch,$fa,$gz) = (split /\./,$file_name,3)[0,1,2];
	my $RESULT_hmm = "$RESULT_hmmdir/$file_name_noch";
	my $DB_hmm = "$root_path/DB/QSG.hmms";
	my $HMM = ();
	my $HMMOPT;
	if($GA==1){
		$HMMOPT= "--cut_ga";
	}elsif(defined $Hmmoptions[0]){
		$HMMOPT = join(" ",@Hmmoptions);
		$HMMOPT=~s/\/-/-/g;
	}
	if($soft_path eq " "){
		$HMM = "hmmscan -o $RESULT_hmm.txt --tblout $RESULT_hmm.tbl --noali --cpu $thread $HMMOPT $DB_hmm $file_path/$file_name";
	}else{
		$HMM = "$soft_path"."hmmer/bin/hmmscan -o $RESULT_hmm.txt --tblout $RESULT_hmm.tbl --cpu $thread --noali $HMMOPT $DB_hmm $file_path/$file_name";
	}
	return $HMM;
}

##diamond
sub predict_QSG_diamond{
	my $file_name = shift @_;
	my $file_path = shift @_;
	my ($file_name_noch,$fa,$gz) = (split /\./,$file_name,3)[0,1,2];
	my $RESULT_diadir = "$output/diamondresults";
	unless(-d $RESULT_diadir || $strategy eq "hmmscan"){
		`mkdir $RESULT_diadir`;
	}
	my $RESULT_dia = "$RESULT_diadir/$file_name_noch.txt";
	my $DB_dia = "$root_path/DB/QSG";
	my $DIAOPT= join(" ",@Diaoptions);
	$DIAOPT=~s/\/-/-/g;
	my $DIA =  "$soft_path"."diamond blastp --db $DB_dia -q $file_path/$file_name -p $thread -o $RESULT_dia -f 6 --max-target-seqs 1 --id $identity $DIAOPT";
	return $DIA;
}

#salmon
sub salmon_abundance{
	my $file_name = shift @_;
	my ($file_name_noch,$fa,$gz) = (split /\./,$file_name,3)[0,1,2];
	my $file_path = shift @_;
	my $rawfile = shift @_;
	my ($f1,$f2,$rawpath)=(split /\t/,$rawfile,3)[0,1,2];
	my $RESULT_saldir = "$output/salmonresults";
	unless(-d $RESULT_saldir){
		`mkdir $RESULT_saldir`;
	}
	my $soft_path2;
	if($soft_path eq " "){
		$soft_path2="$soft_path";
	}else{
		$soft_path2="$soft_path/salmon-1.8.0_linux_x86_64/bin/";
	}
	my $index="$soft_path2"."salmon index -t $file_path/$file_name -p $thread -k $Kmers -i $RESULT_saldir/index_$file_name_noch";
	my $SALOPT= join(" ",@Saloption);
	$SALOPT=~s/\/-/-/g;
	my $quant= "$soft_path"."salmon quant --validateMappings -i $RESULT_saldir/index_$file_name_noch -l A -p $thread --meta -1 $rawpath/$f1 -2 $rawpath/$f2 -o $RESULT_saldir/$file_name_noch.quant $SALOPT";
	my $SAL= "$index\n$quant\n";
	return $SAL;
}
close SH_file;

##########################################################################
#3 run the sh file;
print "########################\n\033[37;1m3 run the sh file.\033[0m\n";
`chmod 764 $output_sh`;
system ($output_sh);

##########################################################################
#4 Cleansing the sample results;
print "########################\n\033[37;1m4 Cleansing the sample results.\033[0m\n";
#Result cleansing sub|union|diamond|hmmscan
foreach my $na(@filess){
	my $cleansing= "perl $bin_path/cleansing.pl -u $strategy -H $output/Hmmscanresults/$na.tbl -D $output/diamondresults/$na.txt -o $output -f $na -S $root_path/DB/QSP_subtype_information.txt";
	system ($cleansing);
}

my $abun_merge;
if(($rawread eq "none")&&($Abundance eq "none")){
	print "No output for gene abundance.";
}elsif(!($Abundance eq "none")){
	$abun_merge = "perl $bin_path/abundance.pl -u $strategy -o $output -A $Abundance";
	system ($abun_merge);
}elsif(!($rawread eq "none")){
	$abun_merge= "perl $bin_path/abundance.pl -u $strategy -o $output -R $rawread";
	system ($abun_merge);
}
#print "$abun_merge";
#merge samples 
##########################################################################
#5. Finsh;
print "\n########################\n\033[37;1m5. Finsh qsg_dcx: $date\033[0m\n";
close STDOUT;


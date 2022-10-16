#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use Array::Utils qw(:all);

my $Usage = <<USAGE;
This program is used for quorum sensing genes prediction.
Author: DAI Chunxiao.
Email: 2446378845\@qq.com
perl $0 -u <Strategy used> -H <Hmmscan result> -D <Diamond result> -o <outputdir> -f <file name> -S <QSP_subtyoe_information> -h";
USAGE

our %input=();
getopts('u:H:D:o:f:S:h',\%input);
if($input{h}){
	die "$Usage\n";
}
my $results= "$input{o}/outputs";
unless(-d $results){
	`mkdir $results`;
}

our %hmmscan = ();
our %diamond = ();
our %subsub = ();
#######################################
open (OUT_SUB,">$results/$input{f}_$input{u}.txt");
open SUBTYPE,"<$input{S}";
while(<SUBTYPE>){
	chomp;
	my @LI = (split /\t/,$_);
	$subsub{$LI[0]}=$LI[1];
	#print "$_\n$subsub{$LI[0]}\n";
}
close SUBTYPE;
#######################################
if($input{u} eq "sub"){
	print OUT_SUB "GeneName\tHmmscanResults\tFullSequenceE-value\tDiamondResults\tSubject\tSubtype\tPident\n";
	&hmmread;
	&diaread;
	my @namedia = keys %diamond;
	my @namehmm = keys %hmmscan;
	my @names = intersect(@namedia,@namehmm);
	for(@names){
		print OUT_SUB "$_\t$hmmscan{$_}\t$diamond{$_}\n";
	}
}elsif($input{u} eq "union"){
	print OUT_SUB "GeneName\tHmmscanResults\tFullSequenceE-value\tDiamondResults\tSubject\tSubtype\tPident\n";
	&hmmread;
	&diaread;
	my @namedia = keys %diamond;
	my @namehmm = keys %hmmscan;
	my @names = unique(@namedia,@namehmm);
	for(@names){
		print OUT_SUB "$_\t$hmmscan{$_}\t$diamond{$_}\n";
	}
}elsif($input{u} eq "diamond"){
	print OUT_SUB "GeneName\tDiamondResults\tSubtype\tPident\n";
	&diaread;
	my @namedia = keys %diamond;
	for(@namedia){
		print OUT_SUB "$_\t$diamond{$_}\n";
	}
}elsif($input{u} eq "hmmscan"){
	print OUT_SUB "GeneName\tHmmscanResults\tFullSequenceE-value\n";
	&hmmread;
	my @namehmm = keys %hmmscan;
	for(@namehmm){
		print OUT_SUB "$_\t$hmmscan{$_}\n";
	}
}
#####################################################
sub hmmread{
	open(MMM,"<$input{H}")||die "$!\n";
	while (<MMM>){
		chomp;
		if($_=~/^#.*/){
			#Line;
		}else{
		$_ =~ s/\s+/\t/g;
		my @hmmscan_re=(split /\t/,$_);
		$hmmscan{$hmmscan_re[2]} = "$hmmscan_re[0]\t$hmmscan_re[4]";
		}
	}
	close MMM;
}

sub diaread{
	open(DIA,"<$input{D}")||die "$!\n";
	while (<DIA>){
		chomp;
		if($_=~/^#\s*/){
			
		}else{
			my @diamond_re=(split /\t/,$_);
			my @type = (split /%/,$diamond_re[1]);
			my $subs;
			if(exists $subsub{$diamond_re[1]}){
				$subs = "$subsub{$diamond_re[1]}";
			}else{
				$subs= "";
			}
			$diamond{$diamond_re[0]} = "$type[0]\t$type[1]\t$subs\t$diamond_re[2]";
		}
	}
	close DIA;
}

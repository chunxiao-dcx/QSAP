#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw/sum/;

my $Usage = <<USAGE;
This program is used for quorum sensing genes prediction.
Author: DAI Chunxiao.
Email: 2446378845\@qq.com
perl $0 -u <Strategy used> -o <outputdir> -R <Rawread list> -A <Gene abundance table> -h";
USAGE

our %input=();
getopts('u:o:R:A:h',\%input);
if($input{h}){
	die "$Usage\n";
}
my $path_cleasing="$input{o}/outputs";
our @list;
our %ABUN;
our %ABUNtotal;
my @qsg = ("AqdC","Clp","CqsA","CqsS","dCACHE","DGC","HdtS","Lactonase","LsrA","LsrB","LsrC","LsrD","LsrG","LsrK","LsrR","LuxI","LuxN","LuxQ","LuxR","LuxS","PDE","PqsD","PqsH","RpfB","RpfC","RpfF","RpfG","LuxM","Acylase","LuxP","PqsR","PhnA","PhnB","PqsA","PqsB","PqsC","PqsE","PqsL");

if($input{R}){
	open IN_R,"<$input{R}";
	while (<IN_R>){
		###list
		our $firstline2th = readline IN_R;
		chomp;
		my($name,$others) = (split /\t/,$_,2)[0,1];
		push @list,$name;
		foreach my $qs(@qsg){
			$ABUN{$name}{$qs}=[];
		}
		###read quant.sf;
		my $path_sf= "$input{o}/salmonresults/$name.quant/quant.sf";
		open SF,"<$path_sf";
		my %salmon;
		while(<SF>){
			chomp;
			my @salmon_re=(split /\t/,$_); 
			$salmon{$salmon_re[0]} = $salmon_re[3];
			#print "$salmon_re[0]\t$salmon_re[3]\n";
		}
		close SF;
		###merge;
		open OUTPUT,"<$path_cleasing/$name\_$input{u}.txt";
		print "$path_cleasing/$name\_$input{u}.txt";
		while(<OUTPUT>){
			chomp;
			my @lines =(split /\t/,$_);
			push @{$ABUN{$name}{$lines[1]}},$salmon{$lines[0]};
		}
		close OUTPUT;
		###caculate
		foreach my $ke(@qsg){
			#print "@{$ABUN{$name}{$ke}}";
			my $sums = sum @{$ABUN{$name}{$ke}};
			$ABUNtotal{$name}{$ke}=$sums;
		}
	}
}elsif($input{A}){
	open IN_AB,"<$input{A}";
	while(<IN_AB>){
		###list;
		chomp;
		my($name,$other) = (split /\t/,$_,2)[0,1];
		print "";
		push @list,$name;
		foreach my $qs(@qsg){
			$ABUN{$name}{$qs}=[];
		}
		###read gene abundance.table;
		open INSS,"<$other";
		our $firstline = readline INSS;
		my @st=(split /\t/,$firstline);
		my($colGene,$colAbundance);
		my $i = 0;
		my $long = @st;
		until($i==$long){
			if($st[$i] eq "Name"){
				$colGene=$i;
				my $r=$i+1;
				print "The column $r contains Gene name in $other\n";
			}elsif($st[$i] eq "Abundance"){
				$colAbundance=$i;
				my $r=$i+1;
				print "The column $r contains Gene abundance in $other\n";
			}
			$i+=1;
		}
		my %table;
		while(<INSS>){
			chomp;
			my @lines =(split /\t/,$_);
			$table{$lines[$colGene]} = $lines[$colAbundance];
		}
		close INSS;
		###merge;
		open OUTPUT,"<$path_cleasing/$name\_$input{u}.txt";
		while(<OUTPUT>){
			chomp;
			my @lines =(split /\t/,$_);
			push @{$ABUN{$name}{$lines[1]}},$table{$lines[0]};
			#print "$table{$lines[0]}";
		}
		###caculate
		foreach my $ke(@qsg){
			my $sums = sum @{$ABUN{$name}{$ke}};
			$ABUNtotal{$name}{$ke}=$sums;
		}
	}
	close IN_AB;
}

open OUTS,">$input{o}/merge_abundance.txt";
print OUTS "QSP";
foreach(@list){
	print OUTS "\t$_";
}
print OUTS "\n";
foreach my $qss(@qsg){
	print OUTS "$qss";
	foreach my $namea(@list){
		print OUTS "\t";
		unless(defined $ABUNtotal{$namea}{$qss}){
			$ABUNtotal{$namea}{$qss}=0;
		}
		print OUTS "$ABUNtotal{$namea}{$qss}";
	}
	print OUTS "\n";
}
close OUTS;
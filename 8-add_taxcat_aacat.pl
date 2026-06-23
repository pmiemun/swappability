#!/usr/bin/perl
use strict;
use Bio::SeqIO;  

#Add a column with species comparison and a column with aa pair

my $input = "desglose_onetoone2.txt";	

open (OUT,">>desglose_onetoone_taxcat.txt");
my %aa = ("A"=>1,"C"=>2,"D"=>3,"E"=>4,"F"=>5,"G"=>6,"H"=>7,"I"=>8,"K"=>9,"L"=>10,"M"=>11,"N"=>12,"P"=>13,"Q"=>14,"R"=>15,"S"=>16,"T"=>17,"V"=>18,"W"=>19,"Y"=>20);
my %cat = ("PROT"=>1,"OBAC"=>2,"ARCH"=>3,"ALVE"=>4,"AMOE"=>5,"EUGL"=>6,"OEUK"=>7,"VIRI"=>8,"FUNG"=>9,"META"=>10);

#Save abbreviations in memory
my %tax = '';	open(IN1,"<../0-db/tax_info.csv");	
while (<IN1>) {	chomp $_;	my @ii = split(/\,/,$_);	$tax{$ii[1]} = $ii[0];	}	close IN1;	

#Go over orthologs
open (IN,"<$input");	my $index = 0;
while (<IN>) { chomp $_;	$index++;		my $line = $_;
	if ($index == 1)	{		print OUT "$line\taapair\ttaxcat\n";	}		next if ($index == 1);			#Only for first line
	
	my @info = split(/\t/,$_);	my $pair = "-";			#Process line

	if ($aa{$info[2]} <= $aa{$info[5]}) { $pair = "$info[2]$info[5]"; } 
	else 																{	$pair = "$info[5]$info[2]"; }			#Add aa pair info
	
	my @sp1 = split(/\|/,$info[0]);	my @sp2 = split(/\|/,$info[3]);		my $taxcat = "-";
	
	if ($cat{$tax{$sp1[0]}} <= $cat{$tax{$sp2[0]}}) { $taxcat = "$tax{$sp1[0]}-$tax{$sp2[0]}"; } 
	else 																						{	$taxcat = "$tax{$sp2[0]}-$tax{$sp1[0]}"; }			#Add taxonomy_category info
	print OUT "$line\t$pair\t$taxcat\n";	
} close IN; close OUT;
exit;

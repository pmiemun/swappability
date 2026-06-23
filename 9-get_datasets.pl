#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;       																																				#BioPerl module to handle the input fasta file    

my $aa1 = "E";	my $aa2 = "E";			my $pair = "EE";

my %seqs;		open my $fa_fh, '<', "filtered_seqs.fasta";	my $curr_id;
while (<$fa_fh>) { chomp; if (/^>(\S+)/) { $curr_id = $1;  $seqs{$curr_id} = '';  } elsif (defined $curr_id) {  $seqs{$curr_id} .= $_; }	}	close $fa_fh;
my $out_fasta = "$aa1\-$pair.fasta";
my $raw1 = "desglose_onetoone.txt";	my $raw2 = "overlaps.txt";
my $ya = "";

#Raw1
open (RAW1,"<$raw1"); open (OUT,">>$out_fasta");
while (<RAW1>) {  chomp $_;  	my @info = split(/\t/,$_);
	if (($info[8] eq "A") || ($info[8] eq "C")) {																				#Only for overlapping polyX
		my $pp1 = "$info[2]$info[5]";		my $pp2 = "$info[5]$info[2]";	
		if (($pp1 eq $pair) || ($pp2 eq $pair)) {																					#Only for current pair
			if ($info[2] eq $aa1) {	next if (index($ya, $info[0]) != -1);	print OUT ">$info[0]\n$seqs{$info[0]}\n";	$ya .= "|$info[0]|"; }
			if ($info[5] eq $aa1) {	next if (index($ya, $info[3]) != -1);	print OUT ">$info[3]\n$seqs{$info[3]}\n";	$ya .= "|$info[3]|"; }	
		}
	}	
}	close RAW1;
#Raw2
open (RAW2,"<$raw2"); open (OUT,">>$out_fasta");
while (<RAW2>) {  chomp $_;  	my @info = split(/\t/,$_);
	if (($info[8] eq "A") || ($info[8] eq "C")) {																				#Only for overlapping polyX
		my $pp1 = "$info[2]$info[5]";		my $pp2 = "$info[5]$info[2]";	
		if (($pp1 eq $pair) || ($pp2 eq $pair)) {																					#Only for current pair
			if ($info[2] eq $aa1) {	next if (index($ya, $info[0]) != -1); print OUT ">$info[0]\n$seqs{$info[0]}\n";	$ya .= "|$info[0]|";	}
			if ($info[5] eq $aa1) {	next if (index($ya, $info[3]) != -1);	print OUT ">$info[3]\n$seqs{$info[3]}\n";	$ya .= "|$info[3]|";	}	
		}
	}
}	close RAW2;
close OUT;
exit;

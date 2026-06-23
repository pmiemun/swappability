#!/usr/bin/perl
use strict;
use Bio::SeqIO;

my %hsa = '';	my $mf_hsa = Bio::SeqIO->new( -file  => "<all.fasta" , '-format' => 'fasta');
while (my $seqq = $mf_hsa->next_seq()) {	my $id = $seqq->display_id; my $qq = $seqq->seq;	$hsa{$id} = $qq; }	

my $input = "relative2.tsv"; 
open (OUT2, ">>polyx_disorder.txt"); print OUT2 "polyx\tid\tposition\tAAinDisorder\tprotlength\tratioINdisorder\n";
open (INx, "<$input");	my $line = 0;
while (<INx>) {	chomp $_;	$line++; next if ($line==1);	print "$line\n";
	my @set = split(/\t/,$_);
	open (OUT,">>aux.fasta"); 	print OUT ">$set[0]\n$hsa{$set[0]}\n";	close OUT;
	my $length = length($hsa{$set[0]});
	my $output = `python3 ./IUPred/iupred3.py aux.fasta short >> aux2.txt`;	unlink "aux.fasta";
	`cut -f3 aux2.txt >> aux3.txt`;		unlink "aux2.txt";
	open (IN,"<aux3.txt"); my $ya = 0; my $disorder = 0; my $init = "IUPRED2";
	while (<IN>) {	chomp $_;
		if (index($_, $init) == 0) { $ya = 1; }
		if ($ya == 1) {	next if (index($_, $init) == 0);	if ($_ > 0.5) { $disorder++;}	}		
	} close IN;	unlink "aux3.txt";
	my $perce = $disorder/$length;
	print OUT2 "$set[2]\t$set[0]\t$set[3]\t$disorder\t$length\t$perce\n";
} close INx; close OUT2;
exit;

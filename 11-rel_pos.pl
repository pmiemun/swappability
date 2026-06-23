#!/usr/bin/perl
use strict;
use Bio::SeqIO;       																											

my $input    = 'polyEE.txt';
my $fasta    = 'aa_seqs_OrthoMCL-CURRENT.fasta';	my $output   = 'relative.tsv';

# read FASTA and save each ID’s length
my %length = '';			
my $mf_hsa = Bio::SeqIO->new( -file  => "<$fasta" , '-format' => 'fasta');
while (my $seqq = $mf_hsa->next_seq()) {	
	my $idx = $seqq->display_id; 	    (my $id = $idx) =~ s/-old//g;
	my $qq = $seqq->seq;	if (length($id) > 0) {	$length{$id} = length($qq);	 }
}

# process input and write output
open my $fh_in,  '<', $input;	open my $fh_out, '>>', $output;	 
#printf $fh_out "id\trelative\tpolyx\tposition\ttaxa\n";

while (<$fh_in>) {	  chomp;
    my @f = split "\t";  my ($id1, $r1, $aa1, $id2, $r2, $aa2, @rest) = @f;   
    my $pair = $rest[-2];  my $taxa = $rest[-1];
    # first polyX
    my ($s1, $e1) = split /-/, $r1;  my $mid1  = ($s1 + $e1) / 2;
    my $rel1  = $length{$id1} > 1
              ? ($mid1 - 1) / ($length{$id1} - 1)
              : 0;
    printf $fh_out "$id1\t$rel1\t$aa1-$pair\t$r1\t$taxa\n";

    # second polyx
    my ($s2, $e2) = split /-/, $r2;  my $mid2  = ($s2 + $e2) / 2;
    my $rel2  = $length{$id2} > 1
              ? ($mid2 - 1) / ($length{$id2} - 1)
              : 0;
    printf $fh_out "$id2\t$rel2\t$aa2-$pair\t$r2\t$taxa\n";
}
close $fh_in;	close $fh_out;

exit;

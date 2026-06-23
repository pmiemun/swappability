#!/usr/bin/perl
use strict;
use warnings;

my $id_file   = 'unique_id.txt';
my $in_fasta  = 'aa_seqs_OrthoMCL-CURRENT.fasta';
my $out_fasta = 'filtered_seqs.fasta';

open my $ID, '<', $id_file;

my %ids;	while (<$ID>) {  chomp;  $ids{$_} = 1 if length;	}	close $ID;

open my $IN,  '<', $in_fasta;	open my $OUT, '>', $out_fasta;

my $keep = 0;
while (my $line = <$IN>) {
    if ($line =~ /^>(\S+)/) {
        my $key = $1;
        if (exists $ids{$key}) { 		$keep = 1;   print $OUT $line; } 
        else {					           	$keep = 0; 							       }
    }
    elsif ($keep) {  print $OUT $line;	   }
}
close $IN;	close $OUT;

exit;

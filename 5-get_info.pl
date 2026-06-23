#!/usr/bin/perl
use strict;
use warnings;

my $pairs_file = "11polyx.txt";
my $polyx_file = "../1-filter/polyx.txt";
my $out_file = "polyx_pair.txt";

# Load polyx data into a hash
open my $px_fh, '<', $polyx_file;		my %polyx;
while (<$px_fh>) {
    chomp;	 next if /^\s*$/;
    my ($id, $start, $end, $aa, @rest) = split /\t/;
    $polyx{$id} = { aa => $aa, start => $start, end => $end };
}
close $px_fh;

# Process each line in 11polyx.txt, look up both IDs, and write merged fields
open my $in_fh,  '<', $pairs_file;		open my $out_fh, '>', $out_file;

print $out_fh "ID1\tAa1\tStart1\tEnd1\tID2\tAa2\tStart2\tEnd2\n";
    
while (<$in_fh>) {
    chomp;	 next if /^\s*$/;   my ($id1, $id2, @rest) = split /\t/;
    # Retrieve fields for ID1 & ID2
    my ($aa1, $s1, $e1) = ('NA','NA','NA');
    if (exists $polyx{$id1}) {  $aa1 = $polyx{$id1}{aa};    $s1  = $polyx{$id1}{start};   $e1  = $polyx{$id1}{end};  }
    my ($aa2, $s2, $e2) = ('NA','NA','NA');
    if (exists $polyx{$id2}) {  $aa2 = $polyx{$id2}{aa};    $s2  = $polyx{$id2}{start};   $e2  = $polyx{$id2}{end};  }

    # Print: ID1, Aa1, Start1, End1, ID2, Aa2, Start2, End2
    print $out_fh join("\t", $id1, $aa1, $s1, $e1, $id2, $aa2, $s2, $e2), "\n";
}
close $in_fh;	close $out_fh;

exit;

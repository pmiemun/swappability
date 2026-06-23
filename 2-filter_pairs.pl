#!/usr/bin/perl
use strict;
use warnings;

my $infile  = 'orthologs.txt';
my $outfile = 'filtered_pairs.txt';

open my $in, '<', $infile;
my @lines; my @scores;

while (my $line = <$in>) {
    chomp $line;
    next unless length $line;
    my @fields = split /\t/, $line;
    push @scores, $fields[2];
    push @lines,  $line;
}
close $in;

my @sorted = sort { $a <=> $b } @scores;
my $n      = @sorted;
my $median;
if ($n % 2) {
    $median = $sorted[int($n / 2)];
} else {
    $median = ($sorted[$n / 2 - 1] + $sorted[$n / 2]) / 2;
}

open my $out, '>', $outfile;
for my $line (@lines) {
    my @fields = split /\t/, $line;
    print $out "$line\n" if $fields[2] > $median;
}
close $out;
print "Median normalized similarity score: $median\n";

exit;

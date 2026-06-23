#!/usr/bin/perl
use strict;
use warnings;

# Count polyX per ID
my %count;	open my $px_fh, '<', '../1-filter/polyx.txt';
while (my $line = <$px_fh>) {  
	chomp $line;	 next unless length $line;  my ($id) = split /\t/, $line;   $count{$id}++;
}
close $px_fh;

open my $out1, '>', '00polyx.txt';	open my $out2, '>', '01polyx.txt';	open my $out3, '>', '11polyx.txt';
open my $out4, '>', 'm1polyx.txt';	open my $out5, '>', 'mmpolyx.txt';	open my $out6, '>', 'm0polyx.txt';

# Process each pair and save in the correct file
open my $pairs_fh, '<', '../1-filter/filtered_pairs.txt';
while (my $line = <$pairs_fh>) {
    chomp $line;   next unless length $line;
    my ($id1, $id2, @rest) = split /\t/, $line;
    my $c1 = $count{$id1} // 0;	 my $c2 = $count{$id2} // 0;  my $suffix = "[$c1|$c2]";
    if ($c1 == 0 and $c2 == 0) {														        print $out1 "$line\t$suffix\n";    }
    elsif (($c1 == 0 and $c2 == 1) or ($c1 == 1 and $c2 == 0)) {    print $out2 "$line\t$suffix\n";	   }
    elsif ($c1 == 1 and $c2 == 1) {													        print $out3 "$line\t$suffix\n";    }
    elsif (($c1 > 1 and $c2 == 1) or ($c1 == 1 and $c2 > 1)) {	    print $out4 "$line\t$suffix\n";    }
    elsif ($c1 > 1 and $c2 > 1) {														        print $out5 "$line\t$suffix\n";    }
    elsif (($c1 == 0 and $c2 > 1) or ($c1 > 1 and $c2 == 0)) {	    print $out6 "$line\t$suffix\n";    }
}
close $pairs_fh;
close $_ for ($out1, $out2, $out3, $out4, $out5, $out6);

exit;

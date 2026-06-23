#!/usr/bin/perl
use strict;
use Bio::SeqIO;  
use File::Temp qw/tempfile/;

				#	A = same_overlap			→ same polyX type in overlapping regions	
				#	C = diff_aligned			→ different polyX type in overlapping regions

# Define files
my $mult_file = "multiple.txt";							my $fasta_file = "../0-db/aa_seqs_OrthoMCL-CURRENT.fasta";	
my $polyx_file = "../1-filter/polyx.txt";		my $out_file = "overlaps.txt";

# Load polyX data
my %polyx; open my $px_fh, '<', $polyx_file;
while (<$px_fh>) {  
		chomp;  next if /^\s*$/;  my ($id, $start, $end, $aa, $seq) = split /\t/;  my $len = $end - $start + 1;
   push @{ $polyx{$id} }, { start => $start, end => $end, aa => $aa, len => $len };
}	close $px_fh;

# Load all sequences into memory
my %seq;		open my $fa_fh, '<', $fasta_file;	my $curr_id;
while (<$fa_fh>) { chomp; if (/^>(\S+)/) { $curr_id = $1;  $seq{$curr_id} = '';  } elsif (defined $curr_id) {  $seq{$curr_id} .= $_; }	}	close $fa_fh;

open my $out_fh, '>', $out_file;	open my $mult_fh, '<', $mult_file;
my $indd = 0;
while (<$mult_fh>) {		$indd++; print "$indd\n";
    chomp;
    next if /^\s*$/;
    my ($id1_raw, $id2_raw, undef, undef) = split /\t/;

    # Fetch sequences
    my $seq1 = $seq{$id1_raw} // '';
    my $seq2 = $seq{$id2_raw} // '';
    next unless $seq1 && $seq2;

    # Fetch polyX lists
    my $px1 = $polyx{$id1_raw} // [];
    my $px2 = $polyx{$id2_raw} // [];
    next unless @$px1 && @$px2;

    # Write combined FASTA to temp
    my ($fh_tmp, $tmp_file) = tempfile( SUFFIX => '.fa' );
    print $fh_tmp ">$id1_raw\n$seq1\n>$id2_raw\n$seq2\n";
    close $fh_tmp;

    # Run MAFFT
    my $align_out = `mafft --quiet $tmp_file 2>/dev/null`;
    unlink $tmp_file;
    chomp $align_out;

    # Parse alignment
    my ($a1, $a2);
    for my $block (split />/, $align_out) {
        next unless $block =~ /\S/;
        my ($hdr, @rows) = split /\n/, $block;
        my $aligned = join('', @rows);
        if ($hdr eq $id1_raw) {
            $a1 = $aligned;
        }
        elsif ($hdr eq $id2_raw) {
            $a2 = $aligned;
        }
    }
    next unless defined $a1 && defined $a2;
    my $alen = length $a1;

    # Build maps: original position -> alignment column
    my @map1;  # index 1..len1 => alignment pos (0-based)
    my @map2;
    {
        my $pos = 0;
        for my $i (0 .. $alen - 1) {
            if (substr($a1, $i, 1) ne '-') {
                $pos++;
                $map1[$pos] = $i;
            }
        }
    }
    {
        my $pos = 0;
        for my $i (0 .. $alen - 1) {
            if (substr($a2, $i, 1) ne '-') {
                $pos++;
                $map2[$pos] = $i;
            }
        }
    }

    # Check each polyX pair for overlap
    for my $r1 (@$px1) {
        for my $r2 (@$px2) {
            my $len1 = $r1->{len};
            my $len2 = $r2->{len};
            my $min_len = $len1 < $len2 ? $len1 : $len2;

            # Gather alignment columns covered by each polyX
            my %cols1;
            for my $pos ($r1->{start} .. $r1->{end}) {
                my $col = $map1[$pos] // next;
                $cols1{$col} = 1;
            }
            my %cols2;
            for my $pos ($r2->{start} .. $r2->{end}) {
                my $col = $map2[$pos] // next;
                $cols2{$col} = 1;
            }

            # Count intersection
            my $overlap = 0;
            for my $c (keys %cols1) {
                $overlap++ if $cols2{$c};
            }
            next unless $overlap * 2 >= $min_len;  # at least 50%

            my $frac  = sprintf("%.2f", $overlap / $min_len);
            my $meanl = sprintf("%.1f", ($len1 + $len2) / 2);
            my $type  = $r1->{aa} eq $r2->{aa} ? 'A' : 'C';

            print $out_fh join("\t",
                $id1_raw,
                "$r1->{start}-$r1->{end}",
                $r1->{aa},
                $id2_raw,
                "$r2->{start}-$r2->{end}",
                $r2->{aa},
                $meanl,
                $frac,
                $type
            ), "\n";
        }
    }
}

close $mult_fh;	close $out_fh;

exit;

#!/usr/bin/perl
use strict;
use Bio::SeqIO;  

#		one_to_one	→ each prot has one polyX
				#	A = same_overlap			→ same polyX type in overlapping regions	
				#	B = same_notaligned		→ same polyX type in not overlapping regions
				#	C = diff_aligned			→ different polyX type in overlapping regions
				#	D = diff_notaligned		→ different polyX type in not overlapping regions

#Load all sequences into memory
my %seqs = ''; my $mf_hsa = Bio::SeqIO->new( -file  => "<aa_seqs_OrthoMCL-CURRENT.fasta" , '-format' => 'fasta');	
while (my $seqq = $mf_hsa->next_seq()) {	my $id = $seqq->display_id;  my $qq = $seqq->seq;	$seqs{$id} = $qq; }	

# Process each pair
my $pair_file = "polyx_pair.txt";		open my $pf, '<', $pair_file;	my $count = 0; my $total = 0;
open (OUTPUT,">>desglose_onetoone.txt");
while (<$pf>) {   chomp;	 next unless length;	$count++; next if ($count==1); $total++;
    my ($id1, $aa1, $s1, $e1, $id2, $aa2, $s2, $e2) = split /\t/;
    open(OUT,">>aux.fasta"); print OUT ">$id1\n$seqs{$id1}\n>$id2\n$seqs{$id2}\n";  close OUT;	 # Write the sequences to a FASTA file
	  my $outfn = "out.fasta"; `mafft --quiet aux.fasta > out.fasta`;	unlink "aux.fasta";			   	 # Run MAFFT
		# Read alignment
    open my $aln_fh, '<', $outfn;  local $/ = "\n>";   my %aln;
    while (<$aln_fh>) {
        s/^>//;  chomp;   my ($hdr, @lines) = split /\n/, $_; my ($id) = split /\s+/, $hdr;  my $aseq = join '', @lines;
        $aseq =~ s/\s+//g;   $aln{$id} = $aseq;
    }  close $aln_fh; unlink $outfn;
		my $aligned1 = $aln{$id1};   my $aligned2 = $aln{$id2};   my $alen = length $aligned1;
	  # For each alignment column, record original position (0 if gap)
    my @map1; my %rev1;   my $pos = 0;
    for my $i (0 .. $alen-1) {
    	my $c = substr($aligned1, $i, 1);
      if ($c ne '-') {   $pos++;   $map1[$i] = $pos;   $rev1{$pos} = $i;   } 
      else {  $map1[$i] = 0;   }
    }
    my @map2; my %rev2;  $pos = 0;
    for my $i (0 .. $alen-1) {
    	my $c = substr($aligned2, $i, 1);
   		if ($c ne '-') { $pos++;   $map2[$i] = $pos;   $rev2{$pos} = $i;   } 
     	else {$map2[$i] = 0;   }
    }
    my $len1 = $e1 - $s1 + 1;  my $len2 = $e2 - $s2 + 1;   	 # PolyX lengths
    # Count overlapping residues based on the shortest polyX
    my $overlap_count = 0;
    if ($len1 <= $len2) {
    	# Iterate each original position in seq1's polyX
      for my $orig ($s1 .. $e1) {
      	next unless exists $rev1{$orig};  my $i = $rev1{$orig};  my $m2 = $map2[$i];
        if ($m2 >= $s2 && $m2 <= $e2) {	 $overlap_count++;	 }
       }
      $overlap_count = $len1 ? $overlap_count / $len1 : 0;   # Fraction of seq1's run that overlaps seq2's run
    } else {
    	for my $orig ($s2 .. $e2) {
    	  next unless exists $rev2{$orig}; my $i = $rev2{$orig};  my $m1 = $map1[$i];
        if ($m1 >= $s1 && $m1 <= $e1) {  $overlap_count++;   }
      }
      $overlap_count = $len2 ? $overlap_count / $len2 : 0;
    }
    # Classify category
    my $same_type = ($aa1 eq $aa2);
    my $overlap   = ($overlap_count >= 0.5);
    my $category;
    if ($same_type) {   $category = $overlap ? 'A' : 'B';  } 
    else {			        $category = $overlap ? 'C' : 'D';  }
    my $mean_len = sprintf("%.1f", ($len1 + $len2) / 2);
    # Print results: ID1, Start1-End1, Aa1, ID2, Start2-End2, Aa2, Mean_length, Overlap, Category
    print OUTPUT "$id1\t$s1-$e1\t$aa1\t$id2\t$s2-$e2\t$aa2\t$mean_len\t$overlap_count\t$category\n";
  	if ($count == 11) { print "$total\n"; $count = 1; }
}
close $pf; close OUTPUT;

exit;

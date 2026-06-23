#!/usr/bin/perl
use strict;
use Bio::SeqIO; 

# Inputs
my $a = "A";
my $aa = "A-AA";
my $input = "polyAA.txt";

# Read FASTA and store the sequence
my $fasta    = 'filtered_seqs.fasta';
my $mf_hsa = Bio::SeqIO->new( -file  => "<$fasta" , '-format' => 'fasta');	my %hsa = '';			
while (my $seqq = $mf_hsa->next_seq()) {	
	my $idx = $seqq->display_id; 	    (my $id = $idx) =~ s/-old//g;
	my $qq = $seqq->seq;	if (length($id) > 0) {	$hsa{$id} = $qq;	 }
}

# Find context per polyX
my $out = "$aa.csv";	open (OUT,">>$out");	print OUT "ac,$aa,length,ppi,M10,M9,M8,M7,M6,M5,M4,M3,M2,M1,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10\n";	close OUT;
open (IN,"<$input");	my $line = 0;	

while(<IN>) { chomp $_;	my @info = split(/\t/,$_);	$line++;	print "$line\n";
	my @info2 = split(/\-/,$info[1]);		my $start1 = $info2[0]; my $end1 = $info2[1];	my $length1 = $end1-$start1+1;
	my @info3 = split(/\-/,$info[4]);		my $start2 = $info3[0]; my $end2 = $info3[1];	my $length2 = $end2-$start2+1;	
	if ($info[2] eq $a) {	&context($info[0],$hsa{$info[0]},$start1,$end1,$length1,$out);	}
	if ($info[5] eq $a) {	&context($info[3],$hsa{$info[3]},$start2,$end2,$length2,$out);	}
} close IN; 
exit;
##############################################
sub context () {		
	my ($id,$fullseq,$start,$end,$length,$out) = @_;
	open (OUT,">>$out");	

	#Get aa from -1 to -10
	my $m1 = "-"; 	my $coord1 = $start-2; if ($coord1 < 0) {} else {	$m1 = substr($fullseq,$coord1,1);		if ($m1 eq "X") { $m1 = "-";	}}
	my $m2 = "-"; 	my $coord2 = $start-3; if ($coord2 < 0) {} else {	$m2 = substr($fullseq,$coord2,1);	  if ($m2 eq "X") { $m2 = "-";	}}
	my $m3 = "-"; 	my $coord3 = $start-4; if ($coord3 < 0) {} else {	$m3 = substr($fullseq,$coord3,1); 	if ($m3 eq "X") { $m3 = "-";	}}
	my $m4 = "-"; 	my $coord4 = $start-5; if ($coord4 < 0) {} else {	$m4 = substr($fullseq,$coord4,1); 	if ($m4 eq "X") { $m4 = "-";	}}
	my $m5 = "-"; 	my $coord5 = $start-6; if ($coord5 < 0) {} else {	$m5 = substr($fullseq,$coord5,1); 	if ($m5 eq "X") { $m5 = "-";	}}
	my $m6 = "-"; 	my $coord6 = $start-7; if ($coord6 < 0) {} else {	$m6 = substr($fullseq,$coord6,1); 	if ($m6 eq "X") { $m6 = "-";	}}
	my $m7 = "-"; 	my $coord7 = $start-8; if ($coord7 < 0) {} else {	$m7 = substr($fullseq,$coord7,1);		if ($m7 eq "X") { $m7 = "-";	}}
	my $m8 = "-"; 	my $coord8 = $start-9; if ($coord8 < 0) {} else {	$m8 = substr($fullseq,$coord8,1);		if ($m8 eq "X") { $m8 = "-";	}}
	my $m9 = "-"; 	my $coord9 = $start-10;if ($coord9 < 0) {} else {	$m9 = substr($fullseq,$coord9,1);		if ($m9 eq "X") { $m9 = "-";	}}
	my $m10 = "-"; 	my $coord10 = $start-11;if($coord10 < 0){} else {	$m10 = substr($fullseq,$coord10,1);	if ($m10 eq "X") { $m10 = "-";	}}
	
	#Get aa from +1 to +10
	my $leng = length($fullseq);
	my $p1 = "-"; 	my $coordp1 = $end;  	 if ($coordp1 > ($leng-1)) {} else {  $p1 = substr($fullseq,$coordp1,1); 	if ($p1 eq "X") { $p1 = "-";	}}
	my $p2 = "-"; 	my $coordp2 = $end+1;  if ($coordp2 > ($leng-1)) {} else {	$p2 = substr($fullseq,$coordp2,1); 	if ($p2 eq "X") { $p2 = "-";	}}
	my $p3 = "-"; 	my $coordp3 = $end+2;  if ($coordp3 > ($leng-1)) {} else {	$p3 = substr($fullseq,$coordp3,1); 	if ($p3 eq "X") { $p3 = "-";	}}
	my $p4 = "-"; 	my $coordp4 = $end+3;  if ($coordp4 > ($leng-1)) {} else {	$p4 = substr($fullseq,$coordp4,1); 	if ($p4 eq "X") { $p4 = "-";	}}
	my $p5 = "-"; 	my $coordp5 = $end+4;  if ($coordp5 > ($leng-1)) {} else {	$p5 = substr($fullseq,$coordp5,1); 	if ($p5 eq "X") { $p5 = "-";	}}
	my $p6 = "-"; 	my $coordp6 = $end+5;  if ($coordp6 > ($leng-1)) {} else {	$p6 = substr($fullseq,$coordp6,1); 	if ($p6 eq "X") { $p6 = "-";	}}
	my $p7 = "-"; 	my $coordp7 = $end+6;  if ($coordp7 > ($leng-1)) {} else {	$p7 = substr($fullseq,$coordp7,1);	if ($p7 eq "X") { $p7 = "-";	}}
	my $p8 = "-"; 	my $coordp8 = $end+7;  if ($coordp8 > ($leng-1)) {} else {	$p8 = substr($fullseq,$coordp8,1);	if ($p8 eq "X") { $p8 = "-";	}}
	my $p9 = "-"; 	my $coordp9 = $end+8;  if ($coordp9 > ($leng-1)) {} else {	$p9 = substr($fullseq,$coordp9,1);	if ($p9 eq "X") { $p9 = "-";	}}
	my $p10 = "-"; 	my $coordp10 = $end+9; if($coordp10 > ($leng-1)) {} else {	$p10 =substr($fullseq,$coordp10,1);	if ($p10 eq "X") { $p10 = "-";	}}

	print OUT "$id,$start-$end,$length,0,$m10,$m9,$m8,$m7,$m6,$m5,$m4,$m3,$m2,$m1,$p1,$p2,$p3,$p4,$p5,$p6,$p7,$p8,$p9,$p10\n";	
	close OUT;
}





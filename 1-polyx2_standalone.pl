#!/usr/bin/perl
use strict;
use Bio::SeqIO;       																																				#BioPerl module to handle the input fasta file    
#################################################################################################################################################################################
##  PolyX Scanner		//	Created by Pablo MIER on 2020		//	Copyright © 2020 Pablo MIER. All rights reserved.
##  Description:
##			The program looks for all polyX in the protein sequence from a given file using a threshold in the form of Y/Z, where Y is the minimum number of same amino 
##			acid in window Z. Overlapping polyX are joined; i.e. AAxxAAAAAAxxAA is one polyA using threshold 8/10, although not all of its windows match the threshold.
##			The default threshold value is 8/10.
##	Usage		--> 		perl polyx2_standalone.pl INPUT MINIMUM MAXIMUM
##	Input 	--> 		INPUT = Fasta file
##	Input2 	--> 		INPUT2 = Minimum number of same amino acid
##	Input3 	--> 		INPUT3 = Window length
##	Output 	--> 		output_polyx.txt		--> List of polyX found in input fasta file, with protein ID and coordinates
#################################################################################################################################################################################
my ($file,$minimum,$window) = @ARGV;																													#Input parameters
if (!$file) 	 { print "A fasta file is needed to start the execution.\n"; exit; }						#Input file is needed
if (!$minimum) { $minimum = 8; }																															#Minimum number of same amino acid. Default = 8
if (!$window)  { $window = 10; }																															#Window length. Window = 10
my $impurities = $window-$minimum;																														#Thresholds, and maximum number of impurities allowed in the polyX
my $polyxregion = "polyx.txt";																																#Output file
my $allresults = '';  	my $already = 0;																											#Variable to temporary save the results, and to print them
open (OUT, ">>$polyxregion");	print OUT "ID\tStart\tEnd\tAa\tpolyX\n"; 												#Print headline in output file
my $multifasta = Bio::SeqIO->new( -file  => "<$file" , '-format' => 'fasta');									#Open input file
while (my $seq = $multifasta->next_seq()) {  																									#Go over the sequences in the input file
	my $na = $seq->display_id;		my $sequence = $seq->seq;																			#Get info for the current protein
	my ($previous_polyx, $previous_na, $previous_start, $previous_end, $previous_aa,$position) = ('','','','','',-1);   #Initialize variables
	while ($position <= length($sequence)) {	$position++;																			#Go over the complete length of the current sequence
		my $win6 = substr($sequence,$position,$window);																						#Retrieve substring with maximum window length
		next if (length($win6) < $minimum);																												#Window must be at least as long as the minimum required amino acids		
		my (%count, $prevalence, $aa, $different, $pos, $out) = ('','','',0,0,0);									#Initialize variables		
		while ($pos < $window) {																																	#Count ocurrences per amino acid forming the current window (fastest way)
			my $t = substr($win6,$pos,1);			$count{$t}++; 	$pos++;																#Get the amino acid in current position and count it; move on to the next position
			$different++ if $count{$t} == 1; 																												#New different amino acid if it is the first time it is counted
			($pos,$out) = ($window,1) if $different > ($impurities+1);															#Finish counting if there are more different amino acids than allowed (number... 	
		} 																																													# ... of allowed impurities + amino acid forming the polyX)
		next if ($out == 1);																																			#Next window if the maximum number of impurities in current window was reached
		my @key_list = keys %count;	my $max_key = pop @key_list;																	#Prepare hash of occurrences
		foreach (@key_list) {  $max_key = $_ if $count{$_} > $count{$max_key};	}									#Locate the most ocurrent amino acid in current window
		($aa, $prevalence) = ($max_key, $count{$max_key});																				#Save the most ocurrent amino acid in current window, and its prevalence
		if ($prevalence >= $minimum) {																														#Consider polyX if the most ocurrent amino acid surpasses the minimum threshold
			my ($done, $initial) = (0, $position);																									#Start coordinate of current polyX
			while ($done == 0) {																																		#Keep on going over the sequence to look for the end of the polyX
				$position++;	%count = ''; my ($prevalence2,$aa2) = ('','');													#Initialize variables
				$win6 = substr($sequence,$position,$window);																					#Retrieve substring with maximum window length
				foreach my $str (split //, $win6) {	$count{$str}++;	}																	#Count ocurrences per amino acid forming the current window
				@key_list = keys %count;	$max_key = pop @key_list;																		#Prepare hash of occurrences
				foreach (@key_list) {  $max_key = $_ if $count{$_} > $count{$max_key};	}							#Locate the most ocurrent amino acid in current window
				($aa2, $prevalence2) = ($max_key, $count{$max_key});																	#Save the most ocurrent amino acid in current window, and its prevalence
				if (($aa2 ne $aa) || ($prevalence2 < $minimum))  { $done = 1;			}	#End of polyX when the most prevalent amino acid is not the same, or if it does not reach the minimum
			}
			my $polyx = substr($sequence, $initial, $position-$initial-1+$window);									#Retrieve polyX region, now that we know when it starts and ends
			if (length($previous_polyx) == 0) 				{																							#Save current polyX, if none has been found yet
				$previous_na = $na;$previous_polyx = $polyx;$previous_aa = $aa;												#Protein ID, polyX sequence, and amino acid forming the polyX
				$previous_start = $initial;$previous_end = $previous_start-1+length($previous_polyx);	#PolyX start and end
			} else {																																								#There is a previous polyX																				
				if (($na eq $previous_na) && ($aa eq $previous_aa) && ($previous_end >= $initial)) {	#Join current with last polyX if: same ID, same amino acid in polyX, they overlap
					$previous_polyx = substr($sequence, $previous_start, $position-$previous_start-1+$window);#Retrieve new joined polyX
					$previous_na = $na; $previous_aa = $aa;																							#Save current polyX details for later: Protein ID and amino acid in polyX
					$previous_end = $previous_start-1+length($previous_polyx);  												#Save current polyX details for later: polyX end
				} else {																																							#If at least one condition fails: print previous polyX, and save current one
					my $chop1 = 0; ($previous_polyx,$chop1) = &chop_polyx($previous_polyx,$previous_aa);#Chop off unspecific beginning and end
					my $start = $previous_start+1+$chop1; my $end = $start-1+length($previous_polyx);		#Calculate polyX start and end coordinates
					$allresults .= "$na\t$start\t$end\t$previous_aa\t$previous_polyx\n"; $already++;		#Keep results in a variable, to reduce time printing them
					if ($already == 10000) { print OUT "$allresults"; $allresults = ''; $already = 0; }	#Print results after 10000 polyX, to liberate memory
					($previous_polyx,$previous_na,$previous_aa,$previous_start) = ($polyx,$na,$aa,$initial); #Save current polyX details for later
					$previous_end = $previous_start-1+length($previous_polyx);													#Save current polyX details for later: polyX end
				}	
			}	
		}	
	}
	if (length($previous_polyx) > 0) {																													#Before going into the next protein, print the polyX that was saved but not printed
		$previous_polyx = substr($sequence, $previous_start, $previous_end-$previous_start+1);		#Retrieve new joined polyX
		my $chop1 = 0; ($previous_polyx,$chop1) = &chop_polyx($previous_polyx,$previous_aa);			#Chop off unspecific beginning and end
		my $start = $previous_start+1+$chop1; my $end = $start-1+length($previous_polyx);					#Calculate polyX start and end coordinates
		$allresults .= "$na\t$start\t$end\t$previous_aa\t$previous_polyx\n"; $already++;					#Keep results in a variable, to reduce time printing them
		if ($already == 10000) { print OUT "$allresults"; $allresults = ''; $already = 0; }				#Print results after 10000 polyX, to liberate memory
		($previous_polyx, $previous_na, $previous_start, $previous_end, $previous_aa) = ('','','','','');		#Reinitialize variables for next protein
	}
} 
print OUT "$allresults";																																			#Print the rest of the results
close OUT;																																										#Close output file
exit;
################################################################################################################################################################################
sub chop_polyx () {				#Chop off unspecific beginning (AQAQQQQQQQ -> QAQQQQQQQ) and end (QAQQQQQQQA -> QAQQQQQQQ)
	my ($polyx,$aa) = @_;																																				#IN -> PolyX region, Amino acid
	my $wq = "-"; 	my $chop1 = 0; 																															#Variables to account for chopped amino acid when calculating the coord
	while ($wq ne $aa) {																																				#Until first letter of polyX is not an X
		my $bf = substr($polyx, 0, 1);																														#Get first character of polyX
		if ($bf eq $aa) { $wq = $aa; }																														#PolyX beginning is an X; stop the loop
		else {	$polyx = substr($polyx, 1); 	$chop1++; }																					#Chop off first character of polyX
	}	$wq = "-";																																								#Reinitialize variable
	while ($wq ne $aa) {																																				#Until last letter of polyX is not an X
		my $leng = length($polyx);	my $oie = $leng - 1;	my $bf2 = substr($polyx, $oie, 1);	 		#Get last character of polyX
		if ($bf2 eq $aa) { $wq = $aa; }																														#PolyX end is an X; stop the loop
		else {	$polyx = substr($polyx,0,$oie); 	 }																							#Chop off last character of polyX
	}
	return($polyx,$chop1);																																			#OUT -> Chopped polyX region, Number of chopped residues at the beginning
}
################################################################################################################################################################################

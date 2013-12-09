#!/usr/bin/env perl -w
#use warnings
# 2010-9-8
# }elsif($Line=~ /Expect\(\d\)\s+=\s+(\S+)/){ => }elsif($Line=~ /Expect\(\d+\)\s+=\s+(\S+)/){

# 2/2/2010

# output homolog sequnces even the given flank length is 0

# 5/28/2009

# filter out alignments whose query loc is between previous alignment query locs

# 5/9/2009

# Set minimal intron = 60 bps

# remove $Min_Gap from parameter list and let it = 10;

# 4/20/2009 version 1.8

# (1) solve the problems caused by missing a repetitive part in the sbjct

# (2) get rid of small alignments that are in big ones for blastn mode

# (3) change $Max_Overlap from 50 to 100 for blastn mode

# (4) devide %B_E into $B_E_For{$Sbject_Begin} and $B_E_Rev{$Sbjcet_End}

# 4/8/2009 version 1.7

# Add a new function: with given genbank file, it can automatically check whether a candidate homolog's accession ID that will be used as the name.

# 3/31/2009 version 1.6

# Add database name to the fasta name of each homologs to seperate homologs from different databases

# 3/24/2009 version 1.5

# Modified the grouping step to make PHI can handle elements with TIRs (for BLASTN)

# 1/19/2009 version 1.4

# Add a parameter to filter out putative pseudogenes

# Change copy_len to copy_score to find the best matching between queries and homologs

# 12/11/2008 version 1.3

# forbid redundancy matches

# 10/17/2008 version 1.2

# 1) if(!$Location_Info{$Name}) -> {if(!(defined($Location_Info{$Name}))) {

# 2) -C $Composition is not supported in linux, delete it

# because the first piece of sequence begins with 0

# Doesn't use genewise as second alignment

# ------------------

# add frameshift detector

# ------------------

# add the DNA copy finding code. NOTE: due to frame shift, the DNA location in the genome may be not very acqurate (+- several bps)

# 5-11/08

# ------------------

# Add e value function

# 5-01/08

# ------------------

# Written by vehell

# 6/8/107

#------

# update to make it can find data for both protein and DNA alignment

#------

# write a single vertion to get rid of long introns, improve it later

#

#-----------------------------------------------------

use Getopt::Std;

#-----------------------------------------------------

getopts("i:t:q:D:P:e:M:d:g:n:R:c:f:p:G:o:h:");



$Input       = defined $opt_i ? $opt_i : "";

$Matrix      = defined $opt_t ? $opt_t : "";

$Query_File  = defined $opt_q ? $opt_q : "";

$Sbjct_File  = defined $opt_D ? $opt_D : "";

$Protein     = defined $opt_P ? $opt_P : 1;

$Max_Evalue  = defined $opt_e ? $opt_e : 0.01; 

$Min_Pro     = defined $opt_M ? $opt_M : 0.7;

$Max_Intron  = defined $opt_d ? $opt_d : 8000;

$Min_Intron  = defined $opt_g ? $opt_g : 60;

#$ID_Tag      = defined $opt_I ? $opt_I : "";

$Max_Output  = defined $opt_n ? $opt_n : 100;

$Realignment = defined $opt_R ? $opt_R : 0;

$Composition = defined $opt_c ? $opt_c : T;

$Flank       = defined $opt_f ? $opt_f : 0;

$Pseu_Filter = defined $opt_p ? $opt_p : 0;

$GenBank     = defined $opt_G ? $opt_G : 0;

$Output      = defined $opt_o ? $opt_o : "protein_copies";

$Help        = defined $opt_h ? $opt_h : "";


usuage() if((!$Input)||($Help));



$Min_Gap = 10;		# a mysterious parameter ...

$Min_Intron = int($Min_Intron/3);

#--------------------------------------------------------

if($Protein == 1) {

	$Max_Overlap = 50;

}else{

	$Max_Overlap = 100;

}





if(($GenBank == 1)&&(-e "$Sbjct_File.gb")) {

	$Locus_Format = 0;

	open(GF, "$Sbjct_File.gb")||die"$!\n";

	while(<GF>) {

		chomp;

		$Line = $_;

		if(/ACCESSION\s+(\S+)\s/) {

			$Accession = $1;

		}

		if(/\s+gene\s+/) {

			if(/\.\./) {

				($Left, $Right) = split(/\.\./, $Line);

				if($Left =~ /(\d+)/) {

					$Gene_Start = $1;

				}else{

					print "$Line can not be parsed\n";

				}



				if($Right =~ /(\d+)/) {

					$Gene_Stop = $1;

				}else{

					print "$Line can not be parsed\n";

				}

			}

		}

		if(/\s+\/gene=\"(\w+\d+)\"/) {

			$Gene_ID = $1;

			$GenBank_Info{$Gene_ID} = $Accession." ".$Gene_Start." ".$Gene_Stop;

		}

		if(/\s+\/locus_tag=\"(\w+\d+)\"/) {

			$Gene_ID = $1;

			$Locus_Info{$Gene_ID} = $Accession." ".$Gene_Start." ".$Gene_Stop;

			$Locus_Format = 1;

		}

	}

	close(GF);



	%GenBank_Info = %Locus_Info if($Locus_Format == 1);



}else{

	$GenBank = 0;

}



#-------------- to seperate homologs from different databases----------------

@Dirs = split(/\//, $Sbjct_File);

$Database_ID = $Dirs[$#Dirs];



#-------------- load query sequences --------------------

open(QF, "$Query_File")||die "$!\n";

while(<QF>) {

	$Line = $_;
    
    #print "$Line\n";

	if(/^>(\S+)/) {

		$Name = $1;
        print "$Line\n";

	}else{

		$Line =~ s/\s+$//g;

		$Query_Name_Seq{$Name} .= $Line;

	}

}

close(QF);



foreach(keys(%Query_Name_Seq)) {

	$Query_Length{$_} = length($Query_Name_Seq{$_});
    
    print "Query length: $Query_Length{$_}\nQuery_Name_Seq: $Query_Name_Seq{$_}\n";
}



#-------------- Loading matrix file ---------------------

if($Protein == 1) {

	open(MF, "$Matrix")||die"$!\n";

	$Line_Num = 0;

	while(<MF>) {

		chomp;

		if($Line_Num == 0) {

			@AAs = split(/\s+/, $_);

		}else{

			@Contents = split(/\s+/, $_);

			for($i = 0; $i < @Contents; $i ++) {

				if($i == 0) {

					$AA = $Contents[0];

				}else{

					$AA_Score{$AA." ".$AAs[$i-1]} = $Contents[$i];

				}

			}

		}

		$Line_Num ++;

	}	

	close(MF);

}else{

	$AA_Score{"A A"} = 1;

	$AA_Score{"A T"} = -1;

	$AA_Score{"A C"} = -1;

	$AA_Score{"A G"} = -1;

	$AA_Score{"A N"} = 0;

	$AA_Score{"T A"} = -1;

	$AA_Score{"T T"} = 1;

	$AA_Score{"T C"} = -1;

	$AA_Score{"T G"} = -1;

	$AA_Score{"T N"} = 0;

	$AA_Score{"C A"} = -1;

	$AA_Score{"C T"} = -1;

	$AA_Score{"C C"} = 1;

	$AA_Score{"C G"} = -1;

	$AA_Score{"C N"} = 0;

	$AA_Score{"G A"} = -1;

	$AA_Score{"G T"} = -1;

	$AA_Score{"G C"} = -1;

	$AA_Score{"G G"} = 1;

	$AA_Score{"G N"} = 0;

	$AA_Score{"N A"} = 0;

	$AA_Score{"N T"} = 0;

	$AA_Score{"N C"} = 0;

	$AA_Score{"N G"} = 0;

	$AA_Score{"N N"} = 0;



	$AA_Score{"a a"} = 1;

	$AA_Score{"a t"} = -1;

	$AA_Score{"a c"} = -1;

	$AA_Score{"a g"} = -1;

	$AA_Score{"a n"} = 0;

	$AA_Score{"t a"} = -1;

	$AA_Score{"t t"} = 1;

	$AA_Score{"t c"} = -1;

	$AA_Score{"t g"} = -1;

	$AA_Score{"t n"} = 0;

	$AA_Score{"c a"} = -1;

	$AA_Score{"c t"} = -1;

	$AA_Score{"c c"} = 1;

	$AA_Score{"c g"} = -1;

	$AA_Score{"c n"} = 0;

	$AA_Score{"g a"} = -1;

	$AA_Score{"g t"} = -1;

	$AA_Score{"g c"} = -1;

	$AA_Score{"g g"} = 1;

	$AA_Score{"g n"} = 0;

	$AA_Score{"n a"} = 0;

	$AA_Score{"n t"} = 0;

	$AA_Score{"n c"} = 0;

	$AA_Score{"n g"} = 0;

	$AA_Score{"n n"} = 0;

}



#-----------------Load blast result---------------------------------

$Database = tblastn_loader($Input);



#---------------- filter out highly possible false homologous alignment (1) small in big (sbject) ------- 

foreach(keys(%Query_Sbjct_Matches)) {

	last if($Protein == 1);				# ---------- this part currently only work for DNA query



	($Query, $Sbjct) = split(/ /, $_);

	@Match_Pairs = split(/   /, $Query_Sbjct_Matches{$Query." ".$Sbjct});

	$Total_Matched_Num = 0;

	foreach(@Match_Pairs) {

		$Total_Matched_Num ++;

	}

	for($i = 0; $i < $Total_Matched_Num; $i ++) {

		($Query_Info, $Sbjct_Info)  = split(/  /, $Match_Pairs[$i]);    

		($S_Bi, $S_M, $S_Ei, $Fi, $Si) = split(/ /, $Sbjct_Info);		# (sbject_Begin, sbject_Match, Sbject_End, Frame, Score)

		for($j = 0; $j < $Total_Matched_Num; $j ++) {

			next if($i == $j);

			next if($Match_Pairs[$j] eq "Filtered");

			

			($Query_Info, $Sbjct_Info)  = split(/  /, $Match_Pairs[$j]);    

			($S_Bj, $S_M, $S_Ej, $Fj, $Sj) = split(/ /, $Sbjct_Info);



			next if($Fi ne $Fj);



			if($Fi == 1) {

				if(($S_Bi >= $S_Bj - $Min_Gap)&&($S_Ei <= $S_Ej + $Min_Gap)) {

					if($Si < $Sj) {

						$Match_Pairs[$i] = "Filtered";

					}

				}

			}else{

				if(($S_Bi <= $S_Bj + $Min_Gap)&&($S_Ei >= $S_Ej - $Min_Gap)) {

					if($Si < $Sj) {

						$Match_Pairs[$i] = "Filtered";

					}

				}

			}

		}

	}



	$Query_Sbjct_Matches{$Query." ".$Sbjct} = ();

	foreach(@Match_Pairs) {

#		print "$_\n";

		$Query_Sbjct_Matches{$Query." ".$Sbjct} .= $_."   " if($_ ne "Filtered");

	}

}



#------------------------------------ find the putative connections among blast hits -------------

foreach(keys(%Query_Sbjcts)) {  

	$Query = $_;
    #print "Query print 2: $Query\n";
    
    $Query_Len = $Query_Length{$Query};
    #if nucleotide search, set max overlap to 25% of query if less than 100
    if($Protein == 0) {
        $Max_overlap_temp = $Query_Len * 0.25;
        if ($Max_overlap_temp < 100) {
            $Max_overlap = $Max_overlap_temp;
        }
        else{
            $Max_overlap = 100;
        }
    }

	@Matched_Sbjcts = split(/ /, $Query_Sbjcts{$_});

	foreach(@Matched_Sbjcts) {

		$Sbjct           = $_;



		%B_E_For         = ();

		%B_E_Rev         = ();

		@S_Begins_Plus   = ();

		@S_Begins_Minus  = ();

		@Begins_Sorted   = ();

		%First_Score     = ();

		#------------------------------------ Use this to check a start loc has been used or not ---------

		%Copy_Start_Loc = ();

		

		$Plus_Match_Num  = 0;

		$Minus_Match_Num = 0;



		@Match_Pairs = split(/   /, $Query_Sbjct_Matches{$Query." ".$Sbjct});



		foreach(@Match_Pairs) {

			($Query_Info, $Sbjct_Info)  = split(/  /, $_);             # $Match_Info = $Q_Begin." ".$Q_Match." ".$Q_End."  ".$S_Begin." ".$S_Match." ".$S_End."   ";

			($Q_Begin, $Q_Match, $Q_End) = split(/ /, $Query_Info);    # $Q_Begin." ".$Q_Match." ".$Q_End

			($S_Begin, $S_Match, $S_End, $Frame, $Score) = split(/ /, $Sbjct_Info);    # $Q_Begin." ".$Q_Match." ".$Q_End $Frame



			if($Frame > 0) {

				$B_E_For{$S_Begin} = $Q_Begin." ".$Q_Match." ".$Q_End." ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame." ".$Score;

				push(@S_Begins_Plus, $S_Begin);

				$Plus_Match_Num ++;

			}else{

				$B_E_Rev{$S_End} = $Q_Begin." ".$Q_Match." ".$Q_End." ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame." ".$Score;

				push(@S_Begins_Minus, $S_End); 

				$Minus_Match_Num ++;

			}

		}



		# Plus Direction: sort the matched segments in small to big order, using sbjct matched locations as hash keys

		@Begins_Sorted = sort {$a <=> $b} @S_Begins_Plus;   

		$Copy_Start = 0;

		for($i = 0; $i < $Plus_Match_Num; $i ++) {

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Score) = split(/ /, $B_E_For{$Begins_Sorted[$i]});

#			print "$Begins_Sorted[$i]: $Q_Begin, $Q_End, $S_Begin, $S_End, $Frame, $Score\n";



			if($Copy_Start == 0) {

				#-------------------------------- find the beginning part of each copy 

				$New_Start  = 1;

				$Copy_Start = 1;

			}else{

				if($S_Begin - $Pre_S_End > $Max_Intron) {

					# ---------------------------------- end of the putative copy, reason: sbjct locs too far away

					$New_Start = 1;

				}elsif($Q_Begin + $Max_Overlap > $Pre_Q_End) {

					# ---------------------------------- connect plus seperated matches(exons)

					if($Protein == 1) {

						$New_Start = 0;

					}else{
                        #check if end of query is > end of previous query plus 10 and start of query is >10 before previous query start
						if($Q_End > $Pre_Q_End + $Min_Gap and ($Q_Begin - $Pre_Q_Begin) > $Min_Gap) {
                            #check if non-overlapping begining of query is longer than the overlapping sequence
                            if (($Q_Begin - $Pre_Q_Begin) > ($Pre_Q_End - $Q_Begin)){

                                $New_Start = 0;
                            }
                            else{
                                $New_Start = 1;
                            }

						}else{

							$New_Start = 1;

						}

					}

				 }else{

					# ---------------------------------- end of the putative copy, reason: new copy start begins

					$New_Start = 1;

				}

			 }

			

			if($New_Start == 1) {

				$Start_Loc = $S_Begin;

				#-----------(begin) check whether the start loc has been used or not ------

				if(defined($Copy_Start_Loc{$Start_Loc})) {

					$Start_Loc_Used = 1;

					while($Start_Loc_Used == 1) {

						$Start_Loc += 0.01;

						if(defined($Copy_Start_Loc{$Start_Loc})) {

							$Start_Loc_Used = 1;

						}else{

							$Start_Loc_Used = 0;

							$Copy_Start_Loc{$Start_Loc} = 1;

						}

					}

				}else{

					$Copy_Start_Loc{$Start_Loc} = 1;

				}

				#-----------(end) check whether the start loc has been used or not ------



				if(defined($Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc})) {

					if($Score > $First_Score{$Start_Loc}) {

						$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} = $B_E_For{$S_Begin}."  ";

						$First_Score{$Start_Loc} = $Score;

					}

				}else{

					if(!(defined($B_E_For{$S_Begin}))) {

						print "not defined: $S_Begin\t$Begins_Sorted[$i]\t$S_Begin\t$S_Match\t$S_End\t$Frame\n";

					}

					$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} = $B_E_For{$S_Begin}."  ";

					$First_Score{$Start_Loc} = $Score;

				}

			}else{

				$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} .= $B_E_For{$S_Begin}."  ";

			}



			$Pre_Q_Begin = $Q_Begin;

			$Pre_S_Begin = $S_Begin;

			$Pre_Q_End   = $Q_End;

			$Pre_S_End   = $S_End;

			$Pre_Dir     = $Dir;

		}



		# Minus Direction: sort the matched segments in small to big order, using sbjct matched locations as hash keys

		@Begins_Sorted = sort {$a <=> $b} @S_Begins_Minus;   

		$Copy_Start = 0;

		for($i = 0; $i < $Minus_Match_Num; $i ++) {

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Score) = split(/ /, $B_E_Rev{$Begins_Sorted[$i]});

#			print "$Begins_Sorted[$i]: $Q_Begin, $Q_End, $S_Begin, $S_End, $Frame, $Score\n";

			

			if($Copy_Start == 0) {

				#-------------------------------- find the beginning part of each copy 

				$New_Start  = 1;

				$Copy_Start = 1;

			}else{

				if($S_End - $Pre_S_Begin > $Max_Intron) {

					# ---------------------------------- end of the putative copy, reason: sbjct locs too far away

					$New_Start = 1;

				}elsif($Pre_Q_Begin + $Max_Overlap > $Q_End) {

					# ---------------------------------- connect minus seperated matches(exons)

					if($Protein == 1) {

						$New_Start = 0;

					}else{

						if($Pre_Q_End > $Q_End + $Min_Gap and ($Pre_Q_Begin - $Q_Begin) > $Min_Gap) {
                            if (($Pre_Q_Begin - $Q_Begin) > ($Q_End - $Pre_Q_Begin)) {
                                
                                $New_Start = 0;
                            }
                            else{
                                $New_Start = 1;
                            }

						}else{

							$New_Start = 1;

						}

					}

				 }else{

					# ---------------------------------- end of the putative copy, reason: new copy start begins

					$New_Start = 1;

				  }

			 }

			

			if($New_Start == 1) {

				$Start_Loc = $S_End;

				#-----------(begin) check whether the start loc has been used or not ------

				if(defined($Copy_Start_Loc{$Start_Loc})) {

					$Start_Loc_Used = 1;

					while($Start_Loc_Used == 1) {

						$Start_Loc += 0.01;

						if(defined($Copy_Start_Loc{$Start_Loc})) {

							$Start_Loc_Used = 1;

						}else{

							$Start_Loc_Used = 0;

							$Copy_Start_Loc{$Start_Loc} = 1;

						}

					}

				}else{

					$Copy_Start_Loc{$Start_Loc} = 1;

				}

				#-----------(end) check whether the start loc has been used or not ------



				if(defined($Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc})) {

					if($Score > $First_Score{$Start_Loc}) {

						$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} = $B_E_Rev{$S_End}."  ";

						$First_Score{$Start_Loc} = $Score;

					}

				}else{

					$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} = $B_E_Rev{$S_End}."  ";

					$First_Score{$Start_Loc} = $Score;

				}

			}else{

				$Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc} .= $B_E_Rev{$S_End}."  ";

			}



			$Pre_Q_Begin = $Q_Begin;

			$Pre_S_Begin = $S_Begin;

			$Pre_Q_End   = $Q_End;

			$Pre_S_End   = $S_End;

			$Pre_Dir     = $Dir;

		}

	}

}



#----------------------- read the location information for reading DNA sequences later ------------------------------

open(SF, "$Sbjct_File")||die "Can not open database!\n";

open(IF, $Sbjct_File.".index")||die "The sbject fasta file has not been indexed?!\n";



while(<IF>) {

      chomp;

      ($Name, $Loc) = split(/ /, $_);

      $Location_Info{$Name} = $Loc;

}      

close(IF);   



#----------------------find the longest/best homolog for each loc in genome (when there are mutiple queries) ...

foreach(keys(%Query_Sbjct_Copy_Infor)) {

	($Query, $Sbjct, $Start_Loc) = split(/ /, $_);

	$Query_Len = $Query_Length{$Query};

	$Copy_Score = 0;

	$Copy_Len   = 0;

	@Matches = split(/  /, $Query_Sbjct_Copy_Infor{$_});

	

	@Query_Match_Region = ();

	for($i = 0; $i < $Query_Len; $i++) {

		$Query_Match_Region[$i] = 0;

	}



	foreach(@Matches) {

		($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Score) = split(/ /, $_);



		$Copy_Score += $Score;



		for($i = $Q_Begin; $i < $Q_End; $i ++) {

			$Query_Match_Region[$i-1] = 1;

		}

	}

	$Stop_Loc = $S_Begin < $S_End ? $S_End : $S_Begin;

	

	foreach(@Query_Match_Region) {

		$Copy_Len ++ if($_ == 1);

	}

	#------------------- check whether there is another homolog copy at the same locus --------------------

	if(defined($Sbjct_Copies{$Sbjct})) {

		@Copies = split(/  /, $Sbjct_Copies{$Sbjct});

		$Same_Num = 0;

		foreach(@Copies) {

			# ------------------  $Query." ".$Copy_Len." ".$Start_Loc." ".$Stop_Loc." ".$Copy_Score 

			($Q, $L, $B, $E, $S) = split(/ /, $_);  

			if(($Start_Loc <= $B)&&($Stop_Loc > $B)) {

				$Same_Loc = 1;

				$Same_Num ++;

			}elsif(($Stop_Loc >= $E)&&($Start_Loc < $E)) {

				$Same_Loc = 1;

				$Same_Num ++;

			}elsif(($Start_Loc >= $B)&&($Stop_Loc <= $E)) {

				$Same_Loc = 1;

				$Same_Num ++;

			}elsif(($Start_Loc <= $B)&&($Stop_Loc >= $E)) {

				$Same_Loc = 1;

				$Same_Num ++;

			}else{

				$Same_Loc = 0;

			}

			

			if($Same_Loc == 1) {

				#------------------------- replace the old copy with the new (longer) copy -----------------

				if(($L <= $Copy_Len)&&($Copy_Len/$Query_Len >= $Min_Pro - 0.2)) {

					if($Same_Num == 1) {

						$Sbjct_Copies{$Sbjct} =~ s/$Q $L $B $E $S/$Query $Copy_Len $Start_Loc $Stop_Loc $Copy_Score/;

					}else{

						$Sbjct_Copies{$Sbjct} =~ s/$Q $L $B $E $S  //;	# if more than one copies have been beaten by the new copy, only replace the first one and remove the others

					}

				}

			}

		}



		if(($Same_Num == 0)&&(($Copy_Len/$Query_Len >= $Min_Pro - 0.2))) {

			$Sbjct_Copies{$Sbjct} .= $Query." ".$Copy_Len." ".$Start_Loc." ".$Stop_Loc." ".$Copy_Score."  ";

		}

	}else{

		if($Copy_Len/$Query_Len >= $Min_Pro - 0.2) {

			$Sbjct_Copies{$Sbjct} = $Query." ".$Copy_Len." ".$Start_Loc." ".$Stop_Loc." ".$Copy_Score."  ";

		}

	}

}



#---------------------- output alignments for each copy --------------------------------

open(GF, ">$Output.align_good")||die"$!\n";

open(BF, ">$Output.align_bad")||die"$!\n";



$Min_Pro_Hiden = $Realignment == 1 ? $Min_Pro - 0.2 : $Min_Pro; # -------------- This is because that a homolog's length may increase after Realignment



foreach(keys(%Sbjct_Copies)) {

	$Sbjct = $_;

	@Copies = split(/  /, $Sbjct_Copies{$_});

	foreach(@Copies) {

		($Query, $Copy_Len, $Start_Loc, $Stop_Loc) = split(/ /, $_);

		$Query_Len = $Query_Length{$Query};

		$Sbjct_Len = $Sbjct_Length{$Sbjct};

		@Matches = split(/  /, $Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Start_Loc});



		if($Copy_Len/$Query_Len < $Min_Pro_Hiden) {

			print (BF "Query: $Query Query length: $Query_Len Sbjct: $Sbjct Sbjct length: $Sbjct_Len\n");

			foreach(@Matches) {

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

				$Alignment = $Query_Sbjct_Aligns{$Query." ".$Sbjct." ".$S_Begin};

				if ($Alignment) {
                    print (BF "$Alignment\n");
                }

				$Query_Sbjct_Aligns{$Query." ".$Sbjct." ".$S_Begin} = ();

			}

		}elsif($Copy_Len/$Query_Len >= $Min_Pro) {

			print (GF "Query: $Query Query length: $Query_Len Sbjct: $Sbjct Sbjct length: $Sbjct_Len\n");

			foreach(@Matches) {

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

				$Alignment = $Query_Sbjct_Aligns{$Query." ".$Sbjct." ".$S_Begin};

				if ($Alignment) {
                    print (GF "$Alignment\n");
                }

				$Query_Sbjct_Aligns{$Query." ".$Sbjct." ".$S_Begin} = ();

			}

		}

	}

}



print (BF "\n---------------------------------\n");

print (BF "The followings are alignments that are filtered out in the second round BLAST\n");

print (BF "---------------------------------\n\n");



#--------------------------- output files -----------------

open(OF, ">$Output.aa")||die"$!\n" if($Protein == 1);

open(OF, ">$Output.dna")||die"$!\n" if($Protein == 0);

open(OF2, ">$Output.list")||die"$!\n";

open(OF3, ">$Output.dna")||die"$!\n" if($Protein == 1);

open(OF4, ">$Output.flank")||die"$!\n";



@OF_AA   = ();

@OF_DNA  = ();

@OF_List = ();

@OF_Len  = ();

$Copy_Num = 0;

$Sub_Flanking = $Max_Intron/2;

#----------------  for each sbjct ------------------

foreach(keys(%Sbjct_Copies)) {

	$Sbjct = $_;

	$Sbjct_Len = $Sbjct_Length{$Sbjct};

	#----------------- for each copy --------------------

	@Copies = split(/  /, $Sbjct_Copies{$Sbjct});

	foreach(@Copies) {

		$Curent_Copy = $_;

		($Query, $Copy_Len, $Copy_Begin, $Copy_End) = split(/ /, $_);

		$Query_Len = $Query_Length{$Query};

		$List = $Query." ".$Query_Len." ".$Database." ".$Sbjct."  ";



		#---------------------- filter out low quality copies --------------------------------

		if($Realignment == 0) {

			next if($Copy_Len/$Query_Len < $Min_Pro);

		}else{

			next if($Copy_Len/$Query_Len < $Min_Pro - 0.2);

		}



		@Matches = split(/  /, $Query_Sbjct_Copy_Infor{$Query." ".$Sbjct." ".$Copy_Begin});



		#---------------------- High resuluation realignment (to find small missing exons)------------

		if($Realignment == 1) {

			$Query_Seq = $Query_Name_Seq{$Query};	

			

			$Copy_First = 1; 

			foreach(@Matches) {

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

#				print "$Query: $Query_Seq: $_\n";

				if($Copy_First == 1) {

					$Dir = $S_Begin > $S_End ? "minus" : "plus";



					$Copy_Query_Begin = $Q_Begin;

					$Copy_Query_End   = $Q_End;

					if($Dir eq "plus") {

						$Copy_Sbjct_Begin = $S_Begin;

						$Copy_Sbjct_End   = $S_End;

					}else{

						$Copy_Sbjct_Begin = $S_End;

						$Copy_Sbjct_End   = $S_Begin;

					}



					$Copy_First = 0;

				}else{

					$Copy_Query_Begin = $Q_Begin < $Copy_Query_Begin ? $Q_Begin : $Copy_Query_Begin;

					$Copy_Query_End   = $Q_End > $Copy_Query_End ? $Q_End : $Copy_Query_End;

					if($Dir eq "plus") {

						$Copy_Sbjct_Begin = $S_Begin < $Copy_Sbjct_Begin ? $S_Begin : $Copy_Sbjct_Begin;

						$Copy_Sbjct_End   = $S_End > $Copy_Sbjct_End ? $S_End : $Copy_Sbjct_End;

					}else{

						$Copy_Sbjct_Begin = $S_End < $Copy_Sbjct_Begin ? $S_End : $Copy_Sbjct_Begin;

						$Copy_Sbjct_End   = $S_Begin > $Copy_Sbjct_End ? $S_Begin : $Copy_Sbjct_End;

					}

				}

			}

			

			$Begin_Range = 0;

			$End_Range   = 0;

#			print "$Copy_Query_Begin\t$Copy_Query_End\t$Copy_Sbjct_Begin\t$Copy_Sbjct_End\n";



			if(($Copy_Query_Begin >= $Min_Gap/2)&&($Protein == 1)) { # -------------- if it miss head ------------

				# -------------- check other copies' locations to void overlapping --------------

				$Head_Knock = 0;

				$Tail_Knock = 0;

				@Copies = split(/  /, $Sbjct_Copies{$Sbjct});

				foreach(@Copies) {

					next if($Curent_Copy eq $_);

					($Q, $CL, $Start, $Stop) = split(/ /, $_);

					if($Dir eq "plus") {

						if(($Copy_Sbjct_Begin > $Stop)&&($Copy_Sbjct_Begin - $Sub_Flanking < $Stop)) {

							$Begin_Range = $Copy_Sbjct_Begin - $Stop + 1;

							$Copy_Sbjct_Begin = $Stop + 1;

							$Head_Knock = 1;

						}

					}else{

						if(($Copy_Sbjct_End < $Start)&&($Copy_Sbjct_End + $Sub_Flanking > $Start)) {

							$End_Range = $Start - $Copy_Sbjct_End + 1;

							$Copy_Sbjct_End = $Start - 1;

							$Tail_Knock = 1;

						}

					}

				}

				if($Dir eq "plus") {

					if($Head_Knock == 0) {

						if($Copy_Sbjct_Begin - $Sub_Flanking <= 0) {

							$Begin_Range = $Copy_Sbjct_Begin;

							$Copy_Sbjct_Begin = 1;

						}else{

							$Copy_Sbjct_Begin = $Copy_Sbjct_Begin - $Sub_Flanking;

							$Begin_Range = $Sub_Flanking;

						}

					}

				}else{

					if($Tail_Knock == 0) {

						if($Copy_Sbjct_End + $Sub_Flanking > $Sbjct_Len) {

							$End_Range = $Sbjct_Len - $Copy_Sbjct_End;

							$Copy_Sbjct_End	= $Sbjct_Len;

						}else{

							$Copy_Sbjct_End	= $Copy_Sbjct_End + $Sub_Flanking;

							$End_Range = $Sub_Flanking;

						}

					}

				}

			}

#			print "$Copy_Sbjct_Begin\t$Copy_Sbjct_End\n";



			if(($Query_Len - $Copy_Query_End >= $Min_Gap/2)&&($Protein == 1)) { # -------------- if it miss some tail ------------

				# -------------- check other copies' locations to void overlapping --------------

				$Head_Knock = 0;

				$Tail_Knock = 0;

				@Copies = split(/  /, $Sbjct_Copies{$Sbjct});

				foreach(@Copies) {

					next if($Curent_Copy eq $_);

					($Q, $CL, $Start, $Stop) = split(/ /, $_);

					if($Dir eq "plus") {

						if(($Copy_Sbjct_End < $Start)&&($Copy_Sbjct_End + $Sub_Flanking > $Start)) {

							$End_Range  = $Start - $Copy_Sbjct_End + 1;

							$Copy_Sbjct_End = $Start - 1;

							$Tail_Knock = 1;

						}

					}else{

						if(($Copy_Sbjct_Begin > $Stop)&&($Copy_Sbjct_Begin - $Sub_Flanking < $Stop)) {

							$Begin_Range  = $Copy_Sbjct_Begin - $Stop + 1;

							$Copy_Sbjct_Begin = $Stop + 1;

							$Head_Knock = 1;

						}

					}

				}

				if($Dir eq "plus") {

					if($Tail_Knock == 0) {

						if($Copy_Sbjct_End + $Sub_Flanking > $Sbjct_Len) {

							$End_Range = $Sbjct_Len - $Copy_Sbjct_End;

							$Copy_Sbjct_End = $Sbjct_Len;

						}else{

							$Copy_Sbjct_End = $Copy_Sbjct_End + $Sub_Flanking;

							$End_Range = $Sub_Flanking;

						}

					}

				}else{

					if($Head_Knock == 0) {

						if($Copy_Sbjct_Begin - $Sub_Flanking <= 0) {

							$Begin_Range = $Copy_Sbjct_Begin;

							$Copy_Sbjct_Begin = 1;

						}else{

							$Copy_Sbjct_Begin = $Copy_Sbjct_Begin - $Sub_Flanking;

							$Begin_Range = $Sub_Flanking;

						}

					}

				}

			}

#			print "$Copy_Sbjct_Begin\t$Copy_Sbjct_End\t$Begin_Range\t$End_Range\n";



			$Sec_Sbjct_Seq = sbject_fasta_picker($Sbjct, $Copy_Sbjct_Begin, $Copy_Sbjct_End - $Copy_Sbjct_Begin + 1);	

#			print "$Sbjct, $Copy_Query_Begin, $Copy_Query_End, $Copy_Sbjct_Begin, $Copy_Sbjct_End\n";



			#----------------- realignment -----------------

			open(TQF, ">$Output.TQ")||die"$!\n";

			print (TQF ">TQ\n$Query_Seq\n");

			close(TQF);



			open(TSF, ">$Output.TS")||die"$!\n";

			print (TSF ">TS\n$Sec_Sbjct_Seq\n");

#			print "$Sec_Sbjct_Seq\n";

			close(TSF);

		

			system("formatdb -i $Output.TS -o F -p F\n");

			system("blastall -m 0 -i $Output.TQ -d $Output.TS -o $Output.TO -p tblastn -e $Max_Evalue -F F -C $Composition\n") if($Protein == 1);

#			system("formatdb -i $Output.TS -o F -p F\n");

#			system("blastall -i $Output.TQ -d $Output.TS -o $Output.TO -p tblastn -e $Max_Evalue -F F -C $Composition\n") if($Protein == 1);



			system("blastall -m 0 -i $Output.TQ -d $Output.TS -o $Output.TO -p blastn -e $Max_Evalue -F F\n") if($Protein == 0);



			$Query_Sbjct_Matches{"TQ TS"} = "none";   # clean the old sub matches

	

			tblastn_loader($Output.".TO");	# load sub tblastn result



			#----------------------- load sub matches that have the same direction of the raw matches

			@Sec_Match_Pairs = split(/   /, $Query_Sbjct_Matches{"TQ TS"});

			if($Query_Sbjct_Matches{"TQ TS"} eq "none") {

				print "Something wrong with the second blast\n";

				print "$Query $Sbjct\n";

				exit(0);

			}



			@Sec_Matches = ();

#			print "D: $Dir\t BR: $Begin_Range\t CQB:$Copy_Query_Begin\t CSE:$Copy_Sbjct_End\tER:$End_Range\n";

			foreach(@Sec_Match_Pairs) {

#				print "$_\n";

				($Sec_Query_Info, $Sec_Sbjct_Info)  = split(/  /, $_);             # $Match_Info = $Q_Begin." ".$Q_Match." ".$Q_End."  ".$S_Begin." ".$S_Match." ".$S_End."   ";

				($Sec_Q_Begin, $Sec_Q_Match, $Sec_Q_End) = split(/ /, $Sec_Query_Info);    # $Q_Begin." ".$Q_Match." ".$Q_End

				($Sec_S_Begin, $Sec_S_Match, $Sec_S_End, $Sec_Frame, $Sec_Score) = split(/ /, $Sec_Sbjct_Info);    # $Q_Begin." ".$Q_Match." ".$Q_End $Frame

				$Sec_Dir = $Sec_S_Begin > $Sec_S_End ? "minus" : "plus";

				#----------------------------- filter out alignments that have wrong directions

				next if($Dir ne $Sec_Dir);

				

				# --------------- filter out outrange matches ----------

				if($Sec_Dir eq "plus") {

					#-------- this outmatch is in the Head range ----

					if(($Sec_S_End < $Begin_Range)&&($Sec_Q_Begin > $Copy_Query_Begin)) {

#						print "1\n";

						next;

					}

					#-------- this outmatch is in the Tail range ----

					if(($Sec_S_Begin > $Copy_Sbjct_End - $Copy_Sbjct_Begin - $End_Range)&&($Sec_Q_End < $Copy_Query_End)) {

#						print "2\n";

						next;

					}

				}else{

					#-------- things outmatch is in the Head range ----

					if(($Sec_S_Begin < $Begin_Range)&&($Sec_Q_End < $Copy_Query_End)) {

#						print "3\n";

						next;

					}

					#-------- this outmatch is in the Tail range ----

					if(($Sec_S_End > $Copy_Sbjct_End - $Copy_Sbjct_Begin - $End_Range)&&($Sec_Q_Begin > $Copy_Query_Begin)) {

#						print "4\n";

						next;

					}

				}

				#----------------------------- filter out small alignments that are in big alignments

				$Small_In_Big = 0;

				foreach(@Sec_Match_Pairs) {

					($Sec_Query_Info_Sub, $Sec_Sbjct_Info_Sub)  = split(/  /, $_);             # $Match_Info = $Q_Begin." ".$Q_Match." ".$Q_End."  ".$S_Begin." ".$S_Match." ".$S_End."   ";

					($Sec_Q_Begin_Sub, $Sec_Q_Match_Sub, $Sec_Q_End_Sub) = split(/ /, $Sec_Query_Info_Sub);    # $Q_Begin." ".$Q_Match." ".$Q_End

					($Sec_S_Begin_Sub, $Sec_S_Match_Sub, $Sec_S_End_Sub, $Sec_Frame_Sub, $Sec_Score_Sub) = split(/ /, $Sec_Sbjct_Info_Sub);    # $Q_Begin." ".$Q_Match." ".$Q_End $Frame

					$Sec_Score_Sub = $Sec_Score_Sub;

					$Sec_Frame_Sub = $Sec_Frame_Sub;

					$Sec_Dir_Sub = $Sec_S_Begin_Sub > $Sec_S_End_Sub ? "minus" : "plus";

					#----------------------------- filter out alignments that have wrong directions

					next if($Dir ne $Sec_Dir_Sub);

					

					if(($Sec_Q_Begin + $Min_Gap > $Sec_Q_Begin_Sub)&&($Sec_Q_End < $Sec_Q_End_Sub + $Min_Gap)) {

						$Small_In_Big = 1 if(length($Sec_Q_Match_Sub) > length($Sec_Q_Match_Sub));

						last;

					}

					if(($Sec_S_Begin + $Min_Gap > $Sec_S_Begin_Sub)&&($Sec_S_End < $Sec_S_End_Sub + $Min_Gap)) {

						$Small_In_Big = 1 if(length($Sec_S_Match_Sub) > length($Sec_S_Match_Sub));;

						last;

					}

				}

				next if($Small_In_Big == 1);



				$Sec_S_Begin += $Copy_Sbjct_Begin - 1;

				$Sec_S_End += $Copy_Sbjct_Begin - 1;

				push(@Sec_Matches, $Sec_Q_Begin." ".$Sec_Q_Match." ".$Sec_Q_End." ".$Sec_S_Begin." ".$Sec_S_Match." ".$Sec_S_End." ".$Sec_Frame." ".$Sec_Score);

			}

			@Matches = @Sec_Matches;

		}



		#----------------------------- sorting exons. in this step, it will find the most reliable path along the query, neglect those false matches ------------

		$Empty = 1;

		foreach(@Matches) {

			$Empty = 0;

		}



		if($Empty == 1) {

			print "it is empty!\n";

			exit(0);

		}

		$Most_Possible_Path = path_finder();

		@Matches = ();

		@Sorted_Alignments = split(/ +/, $Most_Possible_Path);

		foreach(@Sorted_Alignments) {

			next if(($_ eq "Begin")||($_ eq "End"));

			$Match_Infor = $ID_Match{$_};

			push(@Matches, $Match_Infor);

		}



		#------------------  get rid of long "-"s (small intron) for each match --------------------

		if($Protein == 1) {

			@New_Matches = ();

			

			foreach(@Matches) {

				$Match = $_;

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

				$Dir = $S_Begin < $S_End ? "plus" : "minus";



				$Len   = length($Q_Match);

				@Q_AAs = split(//, $Q_Match);

				@S_AAs = split(//, $S_Match);

				@Intron_BPs = ();

				for($i = 0; $i < $Len; $i ++) {

					$Intron_BPs[$i] = 0;

				}



				$Small_Intron = 0;

				for($i = 0; $i < $Len; $i ++) {

					if($i + $Min_Intron < $Len) {

						$Sub_Query = substr($Q_Match, $i, $Min_Intron);

						if($Sub_Query =~ /^-+$/) {

							for($j = $i; $j < $i + $Min_Intron; $j ++) {

								$Intron_BPs[$j] = 1;

							}

							$Small_Intron = 1;

						}

					}

				}

			

				# ------------------------ if there is small introns ----------------------------

				if($Small_Intron == 1) {

					small_intron_killer();

				}else{

					push(@New_Matches, $Match);

				}

			}

		}

		

		#----------------------- dealing with small alignments that are in big ones -------------

		if($Protein == 1) {

			@Matches = redundancy_match_filter(@New_Matches);



		}else{

			@New_Matches = ();

			

			@Bad_Matches = ();

			foreach($i = 0; $i < @Matches; $i ++) {

#				print "$Matches[$i]\n";

				$Bad_Matches[$i] = 0;

			}

		

			foreach($i = 0; $i < @Matches; $i ++) {

				next if($Matches[$i] eq "none");

				($Q_B, $Q_M, $Q_E, $S_B, $S_M, $S_E, $F, $Score) = split(/ /, $Matches[$i]);

				foreach($j = 0; $j < @Matches; $j ++) {

					next if($Bad_Matches[$j] == 1);

					next if($j == $i);

					next if($Matches[$j] eq "none");

					($Q_B_N, $Q_M_N, $Q_E_N, $S_B_N, $S_M_N, $S_E_N, $F_N, $Score_N) = split(/ /, $Matches[$j]);

					if(($Q_B > $Q_B_N)&&($Q_E < $Q_E_N)) {

						$Bad_Matches[$i] = 1;

					}elsif(($Q_B == $Q_B_N)&&($Q_E < $Q_E_N)) {

						$Bad_Matches[$i] = 1;

					}elsif(($Q_B > $Q_B_N)&&($Q_E == $Q_E_N)){

						$Bad_Matches[$i] = 1;

					}elsif(($Q_B == $Q_B_N)&&($Q_E == $Q_E_N)){

						$Bad_Matches[$i] = 1;

					}

				}

			}

			foreach($i = 0; $i < @Matches; $i ++) {

#				print "\nBAD: $Bad_Matches[$i]: $Matches[$i]\n";

				push(@New_Matches, $Matches[$i]) if($Bad_Matches[$i] == 0);

			}



			@Matches = @New_Matches;

		}



		# ---------------------- sort the new matched segments in Query order & filter out matches flanked by two overlapping matches

#		print "@Matches\n";

		$Match_In_Overlap = 1;

		while($Match_In_Overlap == 1) {

			$Match_In_Overlap = 0;



			%B_E = ();

			@Begins = ();

			@Begins_Sorted = ();

			@New_Matches = ();

			$Match_Num   = 0;



			$Copy_Len = 0;

			foreach(@Matches) {

				next if($_ eq "none");

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Score) = split(/ /, $_);

				$B_E{$Q_Begin} = $Q_Begin." ".$Q_Match." ".$Q_End." ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame." ".$Score;

				push(@Begins, $Q_Begin);

				$Match_Num ++;

				$Copy_Len += length($Q_Match);

			}



			@Begins_Sorted = sort {$a <=> $b} @Begins;

			if($Match_Num >= 3) {

				for($i = 1; $i < $Match_Num - 1; $i ++) {

					($Q_B_B, $Q_M_B, $Q_E_B, $S_B_B, $S_M_B, $S_E_B, $F_B, $Score_B) = split(/ /, $B_E{$Begins_Sorted[$i - 1]});

					($Q_B, $Q_M, $Q_E, $S_B, $S_M, $S_E, $F, $Score) = split(/ /, $B_E{$Begins_Sorted[$i]});

					($Q_B_N, $Q_M_N, $Q_E_N, $S_B_N, $S_M_N, $S_E_N, $F_N, $Score_N) = split(/ /, $B_E{$Begins_Sorted[$i + 1]});

					push(@New_Matches, $Q_B_B." ".$Q_M_B." ".$Q_E_B." ".$S_B_B." ".$S_M_B." ".$S_E_B." ".$F_B." ".$Score_B) if($i == 1);

					if(($Q_B < $Q_E_B)&&($Q_E > $Q_B_N)) {

						if($Q_E_B > $Q_B_N) {

							$Match_In_Overlap = 1;

						}else{

							push(@New_Matches, $Q_B." ".$Q_M." ".$Q_E." ".$S_B." ".$S_M." ".$S_E." ".$F." ".$Score);

						}

					}else{

						push(@New_Matches, $Q_B." ".$Q_M." ".$Q_E." ".$S_B." ".$S_M." ".$S_E." ".$F." ".$Score);

					}

				}

				push(@New_Matches, $Q_B_N." ".$Q_M_N." ".$Q_E_N." ".$S_B_N." ".$S_M_N." ".$S_E_N." ".$F_N." ".$Score_N);

				@Matches = @New_Matches;

			}

		}

#		print "@Matches\n";



		# ------------------------- filter out short copies ---------------------

		$Copy_Len = 0;

		foreach(@Matches) {

			next if($_ eq "none");

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

			$Copy_Len += length($Q_Match);			

		}



		if($Copy_Len/$Query_Len < $Min_Pro) {

			foreach(@Matches) {

				next if($_ eq "none");

				($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $_);

				print (BF "Query:$Query Positions:($Q_Begin - $Q_End)     Sbject:$Sbjct Positions($S_Begin - $S_End)\n$Q_Match\n$S_Match\n\n");			

			}

			next;

		}



		# -------------------------------- get rid of the overlapings based on query ---------------------------------

		for($i = 0; $i < $Match_Num; $i ++) {

			$B = $Begins_Sorted[$i];

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $B_E{$B});

			$Dir = $S_Begin > $S_End ? "minus" : "plus";

			$Q_Match_Len = length($Q_Match);

			if($i + 1 < $Match_Num) {

				if($Q_End >= $Begins_Sorted[$i+1]) {

					($Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next) = split(/ /, $B_E{$Begins_Sorted[$i+1]});

#					print "$B_E{$B}\n$B_E{$Begins_Sorted[$i+1]}\n";

					$Q_Match_Len = length($Q_Match);

					$Overlap_Len = $Q_End - $Q_Begin_Next + 1;



					@Q_L_AAs = ();

					@S_L_AAs = ();

					$Q_L_Gap = 0;

					for($j = 0; $j < $Overlap_Len + $Q_L_Gap; $j ++) {

#						print "$Overlap_Len $Q_Match_Len - $j\n";

						$Q_AA = substr($Q_Match, $Q_Match_Len - $j - 1, 1);

						push(@Q_L_AAs, $Q_AA);

						$S_AA = substr($S_Match, $Q_Match_Len - $j - 1, 1);

						push(@S_L_AAs, $S_AA);

						$Q_L_Gap ++ if($Q_AA eq "-");

					}

					@Q_L_AAs = reverse(@Q_L_AAs);

					@S_L_AAs = reverse(@S_L_AAs);

					

					@Q_R_AAs = ();

					@S_R_AAs = ();

					$Q_R_Gap = 0;

					$OverlaP_R = 0;

					for($j = 0; $j < $Overlap_Len + $Q_R_Gap; $j ++) {

						$Q_AA = substr($Q_Match_Next, $j, 1);

						push(@Q_R_AAs, $Q_AA);

						$S_AA = substr($S_Match_Next, $j, 1);

						push(@S_R_AAs, $S_AA);

						$Q_R_Gap ++ if($Q_AA eq "-");

						$OverlaP_R ++;

					}

				

					$Max_Score = -100;



					for($j = 0; $j <= $Overlap_Len; $j ++) {

						$Score = 0;

						#----------- left match ---------------

						$Gap_Len = 0;

						for($k = 0; $k < $Overlap_Len - $j + $Gap_Len; $k ++) {

							if(!(defined($Q_L_AAs[$k]))) {

								#print "$Query\t$Sbjct\n@Matches\n";

								exit(0);

							}

							if(!(defined($S_L_AAs[$k]))) {

								#print "$Query\t$Sbjct\n@Matches\n";

								exit(0);

							}

							if(($Q_L_AAs[$k] eq "-")||($S_L_AAs[$k] eq "-")) {

								$Score += -10;

							}elsif(($Q_L_AAs[$k] eq "*")||($S_L_AAs[$k] eq "*")) {

								$Score += -10;

							}elsif(($Q_L_AAs[$k] eq "X")||($S_L_AAs[$k] eq "X")) {

								$Score += -10;

							}else{

								$Score += $AA_Score{$Q_L_AAs[$k]." ".$S_L_AAs[$k]};

							}

							$Gap_Len ++ if($Q_L_AAs[$k] eq "-");

						}	

						#----------- right match --------------

						$Gap_Len = 0;

						for($k = 0; $k < $j + $Gap_Len; $k ++) {

							if(!(defined($Q_R_AAs[$OverlaP_R - $k - 1]))) {

								#print "$Query\t$Sbjct\n@Matches\n";

								exit(0);

							}

							if(!(defined($S_R_AAs[$OverlaP_R - $k - 1]))) {

								#print "$Query\t$Sbjct\n@Matches\n";

								exit(0);

							}

							if(($Q_R_AAs[$OverlaP_R - $k - 1] eq "-")||($S_R_AAs[$OverlaP_R - $k - 1] eq "-")) {

								$Score += -10;

							}elsif(($Q_R_AAs[$OverlaP_R - $k - 1] eq "*")||($S_R_AAs[$OverlaP_R - $k - 1] eq "*")) {

								$Score += -10;

							}elsif(($Q_R_AAs[$OverlaP_R - $k - 1] eq "X")||($S_R_AAs[$OverlaP_R - $k - 1] eq "X")) {

								$Score += -10;

							}else{

								$Score += $AA_Score{$Q_R_AAs[$OverlaP_R - $k - 1]." ".$S_R_AAs[$OverlaP_R - $k - 1]};

							}

							$Gap_Len ++ if($Q_R_AAs[$OverlaP_R - $k - 1] eq "-");

						}	

						#------------ find the best split site ---------------

						if($Max_Score < $Score) {

							$Left	= $Overlap_Len - $j;

							$Right	= $j;

							$Max_Score = $Score;

						}

#						print "j: $j score $Score\n";

					}

#					print "$Left $Right\n\n\n";

						

					#----------------- trim off the additional AAs of cuurent match------------

					$Q_L_Gap = 0;

					$S_L_Trim = 0;

					for($j = 0; $j < ($Overlap_Len - $Left) + $Q_L_Gap; $j ++) {

						$Q_L_AA = substr($Q_Match, $Q_Match_Len - $j - 1, 1);

						$S_L_AA = substr($S_Match, $Q_Match_Len - $j - 1, 1);

						$Q_L_Gap ++ if($Q_L_AA eq "-");

						$S_L_Trim ++ if($S_L_AA ne "-");

					}

					$Q_End -= $Overlap_Len - $Left;

					if($Protein == 1) {

						$S_End -= $S_L_Trim * 3 if($Dir eq "plus");

						$S_End += $S_L_Trim * 3 if($Dir eq "minus");

					}else{

						$S_End -= $S_L_Trim if($Dir eq "plus");

						$S_End += $S_L_Trim if($Dir eq "minus");

					}

					$Q_Match = substr($Q_Match, 0, $Q_Match_Len - $j);

					$S_Match = substr($S_Match, 0, $Q_Match_Len - $j);



					#----------------- trim off the additional AAs of the next match------------

					$Q_R_Gap = 0;

					$S_R_Trim = 0;

					for($j = 0; $j < ($Overlap_Len - $Right) + $Q_R_Gap; $j ++) {

						$Q_R_AA = substr($Q_Match_Next, $j, 1);

						$S_R_AA = substr($S_Match_Next, $j, 1);

						$Q_R_Gap ++ if($Q_R_AA eq "-");

						$S_R_Trim ++ if($S_R_AA ne "-");

					}

					$Q_Begin_Next += $Overlap_Len - $Right;

					if($Protein == 1) {

						$S_Begin_Next += $S_R_Trim * 3 if($Dir eq "plus");

						$S_Begin_Next -= $S_R_Trim * 3 if($Dir eq "minus");

					}else{

						$S_Begin_Next += $S_R_Trim if($Dir eq "plus");

						$S_Begin_Next -= $S_R_Trim if($Dir eq "minus");

					}

					$Q_Match_Next = substr($Q_Match_Next, $j);

					$S_Match_Next = substr($S_Match_Next, $j);

					

					#----------------- replace the old match with the trimed one ---------

					$B_E{$B} = $Q_Begin." ".$Q_Match." ".$Q_End." ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame;

					#----------------- detect and trim the gaps that are on the end edge

					$Gap_Edge = 1;

					while($Gap_Edge == 1) {

						if($Q_Match =~ /\-$/) {

							$Q_Match = substr($Q_Match, 0, length($Q_Match) - 1);

							$S_Match = substr($S_Match, 0, length($S_Match) - 1);

							if($Protein == 1) {

								$S_End -= 3 if($Dir eq "plus");

								$S_End += 3 if($Dir eq "minus");

							}else{

								$S_End -- if($Dir eq "plus");

								$S_End ++ if($Dir eq "minus");

							}

						}else{

							$Gap_Edge = 0;

						}

					}

					$B_E{$Begins_Sorted[$i]} = $Q_Begin." ".$Q_Match." ".$Q_End." ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame;

					#----------------- detect and trim the gaps that are on the begin edge

					$Gap_Edge = 1;

					while($Gap_Edge == 1) {

						if($Q_Match_Next =~ /^\-/) {

							$Q_Match_Next = substr($Q_Match_Next, 1);

							$S_Match_Next = substr($S_Match_Next, 1);

							if($Protein == 1) {

								$S_Begin_Next += 3 if($Dir eq "plus");

								$S_Begin_Next -= 3 if($Dir eq "minus");

							}else{

								$S_Begin_Next ++ if($Dir eq "plus");

								$S_Begin_Next -- if($Dir eq "minus");

							}

						}else{

							$Gap_Edge = 0;

						}

					}

					$B_E{$Begins_Sorted[$i+1]} = $Q_Begin_Next." ".$Q_Match_Next." ".$Q_End_Next." ".$S_Begin_Next." ".$S_Match_Next." ".$S_End_Next." ".$Frame_Next;

				}

			}

		}





		# -------------------------------- get rid of the overlapings based on query ---------------------------------

		#-------------- sometimes between two exons, the query sequences are continuous but the sbjct sequneces are overlap instead of continous or gap ------------

		#-------------- in this case, the overlaped sequence in the second alignment will be cut off, so be the corrosponding query sequence part, result in a frameshift position

		#-------------- I suspect this overlap problem is so serious that the overlaping region may extend to more than 10 aas, but further improvement like consindering

		# which aligment to be cut off may be necessary.



		for($i = 0; $i < $Match_Num; $i ++) {

			$B = $Begins_Sorted[$i];

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $B_E{$B});

			$Dir = $S_Begin > $S_End ? "minus" : "plus";

			$Q_Match_Len = length($Q_Match);



			if($i + 1 < $Match_Num) {

				($Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next) = split(/ /, $B_E{$Begins_Sorted[$i+1]});



				if($Dir eq "plus") {

					if($S_End > $S_Begin_Next) {

#						print "$Query $Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next\n";

						$Sbjct_Overlap_Len = $S_End - $S_Begin_Next + 1;

						$Query_Overlap_Len = int($Sbjct_Overlap_Len/3);

						$S_Begin_Next += $Sbjct_Overlap_Len;

						$Q_Begin_Next += $Query_Overlap_Len;

						$Q_Match_Next = substr($Q_Match_Next, $Query_Overlap_Len);

						$S_Match_Next = substr($S_Match_Next, $Query_Overlap_Len);

						$B_E{$Begins_Sorted[$i+1]} = $Q_Begin_Next." ".$Q_Match_Next." ".$Q_End_Next." ".$S_Begin_Next." ".$S_Match_Next." ".$S_End_Next." ".$Frame_Next;

					}

				}else{

					if($S_End < $S_Begin_Next) {

#						print "$Query $Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next\n";

						$Sbjct_Overlap_Len = $S_Begin_Next - $S_End + 1;

						$Query_Overlap_Len = int($Sbjct_Overlap_Len/3);

						$S_Begin_Next -= $Sbjct_Overlap_Len;

						$Q_Begin_Next += $Query_Overlap_Len;

						$Q_Match_Next = substr($Q_Match_Next, $Query_Overlap_Len);

						$S_Match_Next = substr($S_Match_Next, $Query_Overlap_Len);

						$B_E{$Begins_Sorted[$i+1]} = $Q_Begin_Next." ".$Q_Match_Next." ".$Q_End_Next." ".$S_Begin_Next." ".$S_Match_Next." ".$S_End_Next." ".$Frame_Next;

					}

				}

			}

		}



		# -------------------------------- find corrosposing DNA sequences  ---------------------------------

		$Q_Protein = "";

		$S_Protein = "";

		$S_DNA     = "";

		for($i = 0; $i < $Match_Num; $i ++) {

			$B = $Begins_Sorted[$i];

			($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $B_E{$B});

			if($S_Begin > $S_End) {

				$Copy_Begin = $Copy_Begin > $S_End ? $S_End : $Copy_Begin;

				$Copy_End   = $Copy_End < $S_Begin ? $S_Begin : $Copy_End;

			}else{

				$Copy_Begin = $Copy_Begin > $S_Begin ? $S_Begin : $Copy_Begin;

				$Copy_End   = $Copy_End < $S_End ? $S_End : $Copy_End;

			}



			$Q_Protein .= $Q_Match;

			$S_Protein .= $S_Match;

			if($Protein == 1) {

				if($S_Begin < $S_End) {

					$S_Match_DNA = sbject_fasta_picker($Sbjct, $S_Begin, $S_End - $S_Begin + 1);

				}else{

					$S_Match_DNA = DNA_reverser(sbject_fasta_picker($Sbjct, $S_End, $S_Begin - $S_End + 1));

				}

				$S_DNA .= $S_Match_DNA;

			}

			$List .= "$Q_Begin $Q_End $S_Begin $S_End $Frame $Q_Match $S_Match  ";

		}

		#----------------------- filter out sequences that have "*"s (putative pseudogenes) ---------------

		if($Pseu_Filter == 1) {

			if($S_Protein =~ /\*/) {

				next;

			}

		}



		#--------------------------- find DNA copies with flanking regions ----------------------

		if(($Copy_Begin - $Flank < 0)||($Copy_End + $Flank > $Sbjct_Len)) {



		}else{

			if($S_Begin < $S_End) {

				$S_Flank_DNA = sbject_fasta_picker($Sbjct, int($Copy_Begin) - $Flank, $Copy_End - int($Copy_Begin) + 1 + $Flank * 2);

			}else{

				$S_Flank_DNA = DNA_reverser(sbject_fasta_picker($Sbjct, int($Copy_Begin) - $Flank, $Copy_End - int($Copy_Begin) + 1 + $Flank * 2));

			}

		}



		#----------------------- filter out premature stop codons and gaps, which may cause errors in tree generation step -----------------

		$S_Protein =~ s/-//g;

		$S_Protein =~ s/\*//g;

		$S_DNA =~ s/-//g;

		$S_DNA =~ s/\*//g;

		#------------------------------------------------------------------

		push(@OF_AA, $S_Protein);

		push(@OF_DNA, $S_DNA) if($Protein == 1);

		push(@OF_List, $List);

		push(@OF_Info, "Query:".$Query." Sbjct:".$Sbjct." Length:".$Copy_Len." Location:(".int($Copy_Begin)." - ".$Copy_End.") Direction:".$Dir);

		$Homolog_Len = abs(int($Copy_Begin) - $Copy_End);

		push(@OF_Len, $Homolog_Len);

		push(@OF_Flank, $S_Flank_DNA);

		$Copy_Num ++;

		last if($Copy_Num >= $Max_Output);

	}

	last if($Copy_Num >= $Max_Output);

}



#--------------------------------- sort the homologs from the longest to the shortest in the output file

$Homolog_Num = 0;

$Longest_AA  = 1;

while($Longest_AA != 0) {

	$Longest_AA = 0;

	$Longest_I  = 0;

	for($i = 0; $i < $Copy_Num; $i ++) {

		$Len = $OF_Len[$i];

		if($Len > $Longest_AA) {

			$Longest_AA = $Len;

			$Longest_I  = $i;

		}

	}

	$OF_Len[$Longest_I] = 0;



	if($Longest_AA != 0) {

		$Homolog_Num ++;

		if($GenBank == 1) {

			$OF_Info[$Longest_I] =~ /Sbjct\:(\S+)\s/;

			$Sbjct      = $1;

			$OF_Info[$Longest_I] =~ /Location\:\((\d+) \- (\d+)\)/;

			$Copy_Begin = $1;

			$Copy_End   = $2;

			foreach(keys(%GenBank_Info)) {

				$Gene_ID = $_;

				$Found   = 0;

				($Accession, $Gene_Start, $Gene_Stop) = split(/ /, $GenBank_Info{$_});

				if($Sbjct =~ /$Accession/) {

					if(($Copy_Begin >= $Gene_Start)&&($Copy_End <= $Gene_Stop)) {

						$Full_Name = $Gene_ID;

						$Found = 1;

					}

				}

				last if($Found == 1);

			}

			if($Found == 0) {

				#$Full_Name = $Database_ID."_".$ID_Tag."_".$Homolog_Num;

				#$Full_Name = $Homolog_Num."_".$Database_ID."_".$ID_Tag;

				$Full_Name = $Homolog_Num."_".$Database_ID;

			}

		}else{

			#$Full_Name = $Database_ID."_".$ID_Tag."_".$Homolog_Num;

			#$Full_Name = $Homolog_Num."_".$Database_ID."_".$ID_Tag;

			$Full_Name = $Homolog_Num."_".$Database_ID;

		}

		#print (OF ">$Full_Name $OF_Info[$Longest_I]\n$OF_AA[$Longest_I]\n");

		print (OF ">$Full_Name \n$OF_AA[$Longest_I]\n");

		print (OF2 "$Full_Name $OF_List[$Longest_I]\n");

		print (OF3 ">$Full_Name $OF_Info[$Longest_I]\n$OF_DNA[$Longest_I]\n") if($Protein == 1);

		print (OF4 ">$Full_Name $OF_Info[$Longest_I]\n$OF_Flank[$Longest_I]\n");

	}

}



close(OF);

close(OF2);

close(OF3);

close(OF4);

close(SF);



close(GF);

close(BF);

#-----------------------------------------------------

#-----------------------------------------------------

sub usuage {

    print "\nHi, need some help?\n";

    print STDERR <<"    _EOT_";



    Usage :tblastn_copy_finder.pl <options> <specification> <default>



	\$Input       = defined \$opt_i ? \$opt_i : "";

	\$Matrix      = defined \$opt_t ? \$opt_t : "/scratch/vehell/Tools/BLOSUM62.txt";

	\$Query_File  = defined \$opt_q ? \$opt_q : "";

	\$Protein     = defined \$opt_P ? \$opt_P : 1;

	\$Database    = defined \$opt_D ? \$opt_D : "";

	\$Max_Evalue  = defined \$opt_e ? \$opt_e : 0.01; 

	\$Min_Pro     = defined \$opt_M ? \$opt_M : 0.7;

	\$Max_Intron  = defined \$opt_d ? \$opt_d : 8000;

	\$Min_Intron  = defined \$opt_g ? \$opt_g : 60;

	#\$ID_Tag      = defined \$opt_I ? \$opt_I : "ZMmar";

	\$Max_Output  = defined \$opt_n ? \$opt_n : 100;

	\$Realignment = defined \$opt_R ? \$opt_R : 0;

	\$Composition = defined \$opt_c ? \$opt_c : T;

	\$Flank       = defined \$opt_f ? \$opt_f : 0;

	\$Pseu_Filter = defined \$opt_p ? \$opt_p : 0;

	\$GenBank     = defined \$opt_G ? \$opt_G : "none";

	\$Output      = defined \$opt_o ? \$opt_o : "protein_copies";

	\$Help        = defined \$opt_h ? \$opt_h : "";



    _EOT_

    exit(1)

}



#-----------------------------------------------------

sub tblastn_loader {

	my($Input) = @_;

	my($Head, $Q_Match, $S_Match, $Line, $Query, $No_Hits, $Query_Len, $Sbjct, $Sbject_Len, $Database);

	my($E_Value, $Alignment_Loader, $Blast_Align, $Frame, $Match_Begin, $Q_Begin, $S_Begin, $Q_End, $S_End);

	open(IF, $Input)||die"$!\n";

	$Head = 1;

	$Q_Match = "";

	$S_Match = "";

	while(<IF>) {

		$Line = $_;

		if($Line =~ /Query= (\S+)/) {

			$New_Query = $1;

			 if($Head == 0) {

				Query_Sbject_Info_Loader($E_Value, $No_Hits, $Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Blast_Align, $Query, $Sbjct, $Score);

				$Head = 1;

			 }

			$Query = $New_Query;

			$No_Hits = 0;

		}

		if($Line =~ /No hits found/) {

			$No_Hits = 1;

		}



		if($Line =~ /^>(\S+)/) {

			if($Head == 0) {

				Query_Sbject_Info_Loader($E_Value, $No_Hits, $Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Blast_Align, $Query, $Sbjct, $Score);

				$Head = 1;

			}

			$Sbjct = $1;

			$Head = 1;

		}

		if($Line =~ /Length = (\d+)/) {

			$Sbject_Len = $1;

			$Sbjct_Length{$Sbjct} = $Sbject_Len;

		}

		if($Line =~ /Database: (\S+)/) {

			$Database = $1;

		}

		if($Line =~ /Score =/) {

			if($Head == 0) {

				Query_Sbject_Info_Loader($E_Value, $No_Hits, $Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Blast_Align, $Query, $Sbjct, $Score);

			}else{

				$Head = 0;

			}

			

			if($Line =~ /Score = +(\S+) +/) {

				$Score = $1;

			}else{

				print "error while loading score\n$_\n";

				exit();

			}



			if($Line =~ /Expect/) {

				if($Line=~ /Expect\s+=\s+(\S+)/) {

					$E_Value = $1;

				}elsif($Line=~ /Expect\(\d+\)\s+=\s+(\S+)/){

					$E_Value = $1;

				}else{

					print "error while loading e value\n$_\n";

					exit();

				}

				$E_Value =~ s/\,//;

			}

		

			$Alignment_Loader = 1;



			$Match_Begin = 1;

			$Q_Match = "";

			$S_Match = "";

		}



		if(($Line =~ /Frame = (\-\d+)/)||($Line =~ /Frame = (\+\d+)/)) {

			$Frame = $1;

		}



		if($Protein == 0) {

			if($Line =~ /Strand = Plus/){

				if($Line =~ /Minus/) {

					$Frame = -1;

				}else{

					$Frame = 1;

				}

			}

		}

		if($Line =~ /Query: (\d+)\s+(\S+) (\d+)/) {

			if($Match_Begin == 1) {

				$Q_Begin     = $1;

			}

			$Q_Match .= $2;

			$Q_End   = $3;

		}

		if(($Line =~ /Sbjct: (\d+)\s+(\S+) (\d+)/)||($Line =~ /Sbjct: (\d+)(\S+) (\d+)/)) {

			if($Match_Begin == 1) {

				$S_Begin     = $1;

				$Match_Begin = 0;

			}

			$S_Match .= $2;

			$S_End   = $3;

		}

		

		if(($Line =~ /TBLASTN/)||($Line =~ /BLASTN/)) {

			$Alignment_Loader = 0;

		}

	

		if($Alignment_Loader == 1) {

			$Blast_Align = $_; 

			$Alignment_Loader = 2;

		}elsif($Alignment_Loader == 2){

			$Blast_Align .= $_;

		}

	}

	Query_Sbject_Info_Loader($E_Value, $No_Hits, $Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Blast_Align, $Query, $Sbjct, $Score);

	close(IF);

	return($Database);

}



#-----------------------------------------------------

sub Query_Sbject_Info_Loader {

	my($E_Value, $No_Hits, $Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame, $Blast_Align, $Query, $Sbjct, $Score) = @_;

	my($E_Qualify, $E_Head, $E_Tail);

	$E_Qualify = 0;

	



	$E_Qualify = E_value_comparison($Max_Evalue, $E_Value);



	if(($No_Hits == 0)&&(($E_Qualify eq ">")||($E_Qualify eq "="))) {

		if(($Query eq "TQ")&&($Sbjct eq "TS")) {

			$Query_Sbjct_Matches{"TQ TS"} = "" if($Query_Sbjct_Matches{"TQ TS"} eq "none");

		}

		# ---------------------------------------- query -> sbjcts

		if(defined($Query_Sbjcts{$Query})) {

			@Sbjcts = split(/ /, $Query_Sbjcts{$Query});

			$Have = 0;

			foreach(@Sbjcts) {

				if ($_ eq $Sbjct) {

					$Have = 1;

					last;

				}

			}

			if($Have == 0) {

				$Query_Sbjcts{$Query} .= $Sbjct." ";

			}

		}else{

			$Query_Sbjcts{$Query} = $Sbjct." ";

		}

		# ---------------------------------------- query_sbjct -> locs

		$Query_Sbjct_Matches{$Query." ".$Sbjct} .= $Q_Begin." ".$Q_Match." ".$Q_End."  ".$S_Begin." ".$S_Match." ".$S_End." ".$Frame." ".$Score."   ";	

		#----------------------------------------- query sbjct S_Start -> blast alignment

		$Query_Sbjct_Aligns{$Query." ".$Sbjct." ".$S_Begin} = $Blast_Align;	

	}

}



#-----------------------------------------------------

sub redundancy_match_filter {

	#--------- this function is designed to deal with the small alignments that are covered by bigger alignments

	#--------- future improvment: add overlap after break the big alignments to find real splice point in the next step

	my(@Matches) = @_;

	my($i, $j, $k, $Dir, $Overlap_Len, $No_Gap_Len, $Q_Match_Next_Len, $Q_AA, $S_AA, $Match_Len, $Q_Gap);

	my($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame);

	my($Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next);

	my(@Q_Overlap_AAs, @S_Overlap_AAs, $Overlap_Score, $Overlap_Score_Next);

	my($Q_Begin_Next_L, $Q_Match_Next_L, $Q_End_Next_L, $S_Begin_Next_L, $S_Match_Next_L, $S_End_Next_L);

	my($Q_Begin_Next_R, $Q_Match_Next_R, $Q_End_Next_R, $S_Begin_Next_R, $S_Match_Next_R, $S_End_Next_R);



	# ------------------  filter out redundancy overlap, which will cause problem for sorting the order of exons

	for($i = 0; $i < @Matches; $i ++) {

		next if($Matches[$i] eq "none");

		($Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame) = split(/ /, $Matches[$i]);	

		for($j = 0; $j < @Matches; $j ++) {

			next if($Matches[$j] eq "none");

			next if($i == $j);

			($Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next) = split(/ /, $Matches[$j]);

			

			if(($Q_Begin >= $Q_Begin_Next)&&($Q_End <= $Q_End_Next)) {  # curent be included in another alignment

				$Dir = $S_Begin > $S_End ? "minus" : "plus";

				$Overlap_Len = $Q_End - $Q_Begin + 1;

				$Overlap_Score_Next = 0;

				$Match_Len = 0;

					

				$No_Gap_Len = 0;

				$Q_Match_Next_Len = length($Q_Match_Next);

#				print "$Overlap_Len\n";

				#------------- caculate the score of the next match ---------------

				for($k = 0; $k < $Q_Match_Next_Len; $k ++) {

					$Q_AA = substr($Q_Match_Next, $k, 1);

					$S_AA = substr($S_Match_Next, $k, 1);

					$No_Gap_Len ++ if ($Q_AA ne "-");

					next if($No_Gap_Len -1 < $Q_Begin - $Q_Begin_Next);

					$Match_Len ++ if($Q_AA ne "-");

					if(($Q_AA eq "-")||($S_AA eq "-")) {

						$Overlap_Score_Next += -10;

					}elsif(($Q_AA eq "*")||($S_AA eq "*")) {

						$Overlap_Score_Next += -10;

					}elsif(($Q_AA eq "X")||($S_AA eq "X")) {

						$Overlap_Score_Next += -10;

					}else{

						$Overlap_Score_Next += $AA_Score{$Q_AA." ".$S_AA};

						if(!(defined($AA_Score{$Q_AA." ".$S_AA}))) {

							#print "$Q_AA. .$S_AA\n";

							#print "$Query, $Sbjct, @Matches, \n $Q_Begin, $Q_End, $S_Begin, $S_End, $Frame\n";

							exit(0);

						}

					}

#					print "$k\t$Q_AA $S_AA\t$Overlap_Score_Next\n";

					if($Match_Len == $Overlap_Len) {

						last;

					}

				}



				#------------- caculate the score of the next match ---------------

				@Q_Overlap_AAs = split(//, $Q_Match);

				@S_Overlap_AAs = split(//, $S_Match);

				$Overlap_Score = 0;



				for($k = 0; $k < @Q_Overlap_AAs; $k ++) {

					if(($Q_Overlap_AAs[$k] eq "-")||($S_Overlap_AAs[$k] eq "-")) {

						$Overlap_Score += -10;

					}elsif(($Q_Overlap_AAs[$k] eq "*")||($S_Overlap_AAs[$k] eq "*")) {

						$Overlap_Score += -10;

					}elsif(($Q_Overlap_AAs[$k] eq "X")||($S_Overlap_AAs[$k] eq "X")) {

						$Overlap_Score += -10;

					}else{

						$Overlap_Score += $AA_Score{$Q_Overlap_AAs[$k]." ".$S_Overlap_AAs[$k]};

					}

				}

#				print "$Q_Begin, $Q_Match, $Q_End, $S_Begin, $S_Match, $S_End, $Frame\n$Q_Begin_Next, $Q_Match_Next, $Q_End_Next, $S_Begin_Next, $S_Match_Next, $S_End_Next, $Frame_Next\n$Overlap_Score > $Overlap_Score_Next\n";



				#------------- make a judge which one to keep ---------------

				if($Overlap_Score <= $Overlap_Score_Next) {

					# ------------------- get rid of the small match ----------------

					$Matches[$i] = "none";

				}else{

					# ------------------- get the info of the sub next match left part -------------

					$Q_Begin_Next_L = $Q_Begin_Next;

					$S_Begin_Next_L = $S_Begin_Next;

					$Q_Match_Next_L = "";

					$S_Match_Next_L = "";

					$Q_End_Next_L = $Q_Begin - 1;

					$S_End_Next_L = $S_Begin_Next;

					$Q_Gap = 0;



					for($p = 0; $p < $Q_Begin - $Q_Begin_Next + $Q_Gap; $p ++) {

						$Q_AA = substr($Q_Match_Next, $p, 1);

						$S_AA = substr($S_Match_Next, $p, 1);

						$Q_Match_Next_L .= $Q_AA;

						$S_Match_Next_L .= $S_AA;

						if($Dir eq "plus") {

							$S_End_Next_L += 3 if($S_AA ne "-")

						}else{

							$S_End_Next_L -= 3 if($S_AA ne "-")

						}

						$Q_Gap ++ if($Q_AA eq "-");

					}

#					print "$Q_Begin_Next_L..$Q_End_Next_L. .$S_Begin_Next_L..$S_End_Next_L. .$Frame_Next\n";



					# ------------------- get the info of the sub next match right part -------------

					$Q_Begin_Next_R = $Q_Begin_Next;

					$S_Begin_Next_R = $S_Begin_Next;

					$Q_Match_Next_R = "";

					$S_Match_Next_R = "";

					$Q_End_Next_R = $Q_End_Next;

					$S_End_Next_R = $S_End_Next;

						

					$No_Gap_Len = 0;

					$Match_Len  = 0;

					for($p = 0; $p < $Q_Match_Next_Len; $p ++) {

						$Q_AA = substr($Q_Match_Next, $p, 1);

						$S_AA = substr($S_Match_Next, $p, 1);

						$No_Gap_Len ++ if ($Q_AA ne "-");

						$Q_Begin_Next_R ++ if($Q_AA ne "-");

						if($Dir eq "plus") {

							$S_Begin_Next_R += 3 if($S_AA ne "-")

						}else{

							$S_Begin_Next_R -= 3 if($S_AA ne "-")

						}

						next if($No_Gap_Len -1 < $Q_Begin - $Q_Begin_Next);

						$Match_Len ++ if($Q_AA ne "-");

						if($Match_Len == $Overlap_Len) {

							$Q_Match_Next_R = substr($Q_Match_Next, $p+1);

							$S_Match_Next_R = substr($S_Match_Next, $p+1);

							last;

						}

					}

#					print "$Q_Begin_Next_R. .$Q_End_Next_R. .$S_Begin_Next_R. .$S_End_Next_R\n";



					# ------------------- judge which part to keep  ------------------

					if($Dir eq "plus") {

						if($S_End < $S_Begin_Next_R) {

							$Matches[$j] = $Q_Begin_Next_R." ".$Q_Match_Next_R." ".$Q_End_Next_R." ".$S_Begin_Next_R." ".$S_Match_Next_R." ".$S_End_Next_R." ".$Frame_Next;

						}elsif($S_Begin > $S_End_Next_R) {

							$Matches[$j] = $Q_Begin_Next_L." ".$Q_Match_Next_L." ".$Q_End_Next_L." ".$S_Begin_Next_L." ".$S_Match_Next_L." ".$S_End_Next_L." ".$Frame_Next;

						}else{

							$Matches[$i] = "none";

						}

					}else{

						if($S_Begin < $S_End_Next_R) {

							$Matches[$j] = $Q_Begin_Next_L." ".$Q_Match_Next_L." ".$Q_End_Next_L." ".$S_Begin_Next_L." ".$S_Match_Next_L." ".$S_End_Next_L." ".$Frame_Next;

						}elsif($S_End > $S_Begin_Next_R) {

							$Matches[$j] = $Q_Begin_Next_R." ".$Q_Match_Next_R." ".$Q_End_Next_R." ".$S_Begin_Next_R." ".$S_Match_Next_R." ".$S_End_Next_R." ".$Frame_Next;

						}else{

							$Matches[$i] = "none";

						}

					}

				}

			}

		}

	}

	return(@Matches);

}



#-----------------------------------------------------

sub sbject_fasta_picker {

    my($Name, $Loc, $Len) = @_;

	my($Seq, $Zero_Loc);

    if(!(defined($Location_Info{$Name}))) {

		print "$Name\n";

       die "can not find sbject fasta seq\n";

    }

    $Zero_Loc = $Location_Info{$Name};

    sysseek(SF, $Loc + $Zero_Loc - 1, 0); 

    sysread(SF, $Seq, $Len);



    return($Seq);

}



#-----------------------------------------------------

sub DNA_reverser {

    my($Seq) = @_;

	$Seq = reverse $Seq;

	$Seq =~ tr/ACGTacgt/TGCAtgca/;

    return($Seq);

}



#-----------------------------------------------------

sub small_intron_killer {

	$Q_Gap = 0;

	$S_Gap = 0;

	$Sub_Q_Match = "";

	$Sub_S_Match = "";

	$Sub_Begin = 1;

	for($i = 0; $i < $Len; $i ++) {

		if($Intron_BPs[$i] == 1) {

			if($Sub_Begin == 0) {

				$Sub_Q_End = $Q_Begin + ($i - 1) - $Q_Gap;

				$Non_AA_Num = ($Sub_Q_Match =~ tr/\-|\*//);



				if($Dir eq "plus") {

					$Sub_S_End = $S_Begin + ($i - $S_Gap)*3 - 1;

				}else{

					$Sub_S_End = $S_Begin - ($i - $S_Gap)*3 + 1;

				}



				$Gap_Edge = 1;

				while($Gap_Edge == 1) {

					if($Sub_Q_Match =~ /^\-/) {

						$Sub_Q_Match = substr($Sub_Q_Match, 1);

						$Sub_S_Match = substr($Sub_S_Match, 1);

						$Sub_S_Begin += 3 if($Dir eq "plus");

						$Sub_S_Begin += 3 if($Dir eq "minus");

					}else{

						$Gap_Edge = 0;

					}

				}

				$Score = score_caculater($Sub_Q_Match, $Sub_S_Match);

				push(@New_Matches, $Sub_Q_Begin." ".$Sub_Q_Match." ".$Sub_Q_End." ".$Sub_S_Begin." ".$Sub_S_Match." ".$Sub_S_End." ".$Frame." ".$Score);



				$Sub_Q_Match = "";

				$Sub_S_Match = "";

				$Sub_Begin = 1;

			}

		}else{

			if($Sub_Begin == 1) {

				$Sub_Q_Begin = $Q_Begin + $i - $Q_Gap;

				if($Dir eq "plus") {

					$Sub_S_Begin = $S_Begin + ($i - $S_Gap)*3;

				}else{

					$Sub_S_Begin = $S_Begin - ($i - $S_Gap)*3;

				}

				$Sub_Begin   = 0;

			}

			$Sub_Q_Match .= $Q_AAs[$i];

			$Sub_S_Match .= $S_AAs[$i];;

		}

		$Q_Gap ++ if($Q_AAs[$i] eq "-");

		$S_Gap ++ if($S_AAs[$i] eq "-");

	}

	

	# ------------------- last submatch -----------------------

	$Sub_Q_End = $Q_Begin + ($i - 1) - $Q_Gap;

	$Non_AA_Num = ($Sub_Q_Match =~ tr/\-|\*//);



	if($Dir eq "plus") {

		$Sub_S_End = $S_Begin + ($i - $S_Gap)*3 - 1;

	}else{

		$Sub_S_End = $S_Begin - ($i - $S_Gap)*3 + 1;

	}

					

	$Gap_Edge = 1;

	while($Gap_Edge == 1) {

		if($Sub_Q_Match =~ /^\-/) {

			$Sub_Q_Match = substr($Sub_Q_Match, 1);

			$Sub_S_Match = substr($Sub_S_Match, 1);

			$Sub_S_Begin += 3 if($Dir eq "plus");

			$Sub_S_Begin -= 3 if($Dir eq "minus");

		}else{

			$Gap_Edge = 0;

		}

	}

	$Score = score_caculater($Sub_Q_Match, $Sub_S_Match);

	push(@New_Matches, $Sub_Q_Begin." ".$Sub_Q_Match." ".$Sub_Q_End." ".$Sub_S_Begin." ".$Sub_S_Match." ".$Sub_S_End." ".$Frame." ".$Score);



	$Sub_Q_Match = "";

	$Sub_S_Match = "";

	$Sub_Begin = 1;

}



#-------------------------------------------

sub E_value_comparison {

	my($E1, $E2) = @_;

	my($Compare_Result, $E1_Head, $E1_Tail, $E2_Head, $E2_Tail);

	$Compare_Result = ">";

	if($E1 =~ /e-/) {

		if($E1 =~ /(\d+)e-(\d+)/) {

			$E1_Head = $1;

			$E1_Tail = $2;

		}elsif($E1 =~ /e-(\d+)/){

			$E1_Head = 1;

			$E1_Tail = $1;

		}else{

			print "Strange E value\n";

			exit(0);

		}



		if($E2 =~ /e-/) {

			if($E2 =~ /(\d+)e-(\d+)/) {

				$E2_Head = $1;

				$E2_Tail = $2;

			}elsif($E2 =~ /e-(\d+)/){

				$E2_Head = 1;

				$E2_Tail = $1;

			}else{

				print "Strange E value\n";

				exit(0);

			}

	

			if($E1_Tail > $E2_Tail) {

				$Compare_Result = "<";

			}elsif($E1_Tail < $E2_Tail) {

				$Compare_Result = ">";

			}else{

				if($E1_Head > $E2_Head){

					$Compare_Result = ">";

				}elsif($E1_Head < $E2_Head){

					$Compare_Result = "<";

				}else{

					$Compare_Result = "=";

				}

			}

		}else{

			$Compare_Result = "<";

		}

	}else{

		if($E2 =~ /e-/) {

			$Compare_Result = ">";

		}else{

			if($E1 > $E2) {

				$Compare_Result = ">";

			}elsif($E1 < $E2){

				$Compare_Result = "<";

			}else{

				$Compare_Result = "=";

			}

		}

	}



	return($Compare_Result);

}



# --------------------------------------

sub path_finder {

	# ---------------- find the start alignment as the seed of paths--------

	$Max_Score = 0;

	%ID_Match = ();

	%ID_Score = ();

	@Paths = ();

	for($i = 0; $i < @Matches; $i ++) {

		$ID_Match{$i} = $Matches[$i];

#		print "Path $i : $Matches[$i]\n";

		($Q_B, $Q_M, $Q_E, $S_B, $S_M, $S_E, $Dir, $Score) = split(/ /, $Matches[$i]);

		$ID_Score{$i} = $Score;

		if($Score > $Max_Score) {

			$Max_Score = $Score;

			$Seed_ID   = $i;

		}

	}

#	print "Seed: $Seed_ID\n";

	push(@Paths, " ".$Seed_ID." ");



	# ----------- extend to the begining --------

	$Head = 0;

	while($Head == 0) {

		@New_Paths = ();

		$Head = 1;

		foreach(@Paths) {

			$Path = $_;

#			print "Head Paths: $_\n";

			$Path_End = 1;

			if($Path =~ /^Begin /) {

				push(@New_Paths, $Path);

				next;

			}

			$Path =~ /^ (\d+) /;

			$Current_ID = $1;



			($Q_B, $Q_M, $Q_E, $S_B, $S_M, $S_E, $Dir, $Score) = split(/ /, $ID_Match{$Current_ID});

			$Dir = $Dir > 0 ? "plus" : "minus";

			for($i = 0; $i < @Matches; $i ++) {

				next if($Path =~ / $i /);

				($Q_BN, $Q_MN, $Q_EN, $S_BN, $S_MN, $S_EN, $Dir_N, $Score_N) = split(/ /, $Matches[$i]);

				$Dir_N = $Dir_N > 0 ? "plus" : "minus";

				next if($Dir_N ne $Dir);



				if($Dir eq "plus") {

					if($S_BN >= $S_B) {

						next;

					}elsif(($S_BN < $S_B)&&($S_EN >= $S_E)) {  # may be too loose

						next;

					}

				}else{

					if($S_BN <= $S_B) {

						next;

					}elsif(($S_BN > $S_B)&&($S_EN <= $S_E)) {  # may be too loose

						next;

					}

				}



				if($Q_EN <= $Q_B) {

					push(@New_Paths, " ".$i." ".$Path);

					$Path_End = 0;

					$Head = 0;

				}elsif(($Q_BN < $Q_B)&&($Q_EN > $Q_B)) {

					push(@New_Paths, " ".$i." ".$Path);

					$Path_End = 0;

					$Head = 0;

				}elsif(($Q_BN <= $Q_B)&&($Q_EN >= $Q_E)) {

					next;

				}elsif(($Q_BN >= $Q_B)&&($Q_EN <= $Q_E)) {

					next;

				}elsif(($Q_BN >= $Q_B)&&($Q_EN > $Q_E)) {

					next;

				}elsif($Q_BN >= $Q_E) {

					next;

				}else{

					print "unknown relationship\n$ID_Match{$Current_ID}\n$Matches[$i]\n";

					exit(0);

				}

			}



			if($Path_End == 1) {

				push(@New_Paths, "Begin ".$Path);

			}

		}

		@Paths = @New_Paths;

#		print "\n";



	}



	# ----------- extend to the end -------------

	$Tail = 0;

	while($Tail == 0) {

		@New_Paths = ();

		$Tail = 1;

		foreach(@Paths) {

			$Path = $_;

#			print "Tail Paths: $_\n";

			$Path_End = 1;

			if($Path =~ /End$/) {

				push(@New_Paths, $Path);

				next;

			}

			$Path =~ / (\d+) $/;

			$Current_ID = $1;



			($Q_B, $Q_M, $Q_E, $S_B, $S_M, $S_E, $Dir, $Score) = split(/ /, $ID_Match{$Current_ID});

			$Dir = $Dir > 0 ? "plus" : "minus";



			for($i = 0; $i < @Matches; $i ++) {

				next if($Path =~ / $i /);

				($Q_BN, $Q_MN, $Q_EN, $S_BN, $S_MN, $S_EN, $Dir_N, $Score_N) = split(/ /, $Matches[$i]);

				$Dir_N = $Dir_N > 0 ? "plus" : "minus";

				next if($Dir_N ne $Dir);



				if($Dir eq "plus") {

					if($S_EN <= $S_E) {

						next;

					}elsif(($S_EN > $S_E)&&($S_BN <= $S_B)) {  # may be too loose

						next;

					}

				}else{

					if($S_EN >= $S_E) {

						next;

					}elsif(($S_EN < $S_E)&&($S_BN >= $S_B)) {  # may be too loose

						next;

					}

				}



				if($Q_EN <= $Q_B) {

					next;

				}elsif(($Q_BN < $Q_B)&&($Q_EN >= $Q_B)) {

					next;

				}elsif(($Q_BN <= $Q_B)&&($Q_EN >= $Q_E)) {

					next;

				}elsif(($Q_BN >= $Q_B)&&($Q_EN <= $Q_E)) {

					next;

				}elsif(($Q_BN >= $Q_B)&&($Q_EN > $Q_E)) {

					push(@New_Paths, $Path." ".$i." ");

					$Path_End = 0;

					$Tail = 0;

				}elsif($Q_BN >= $Q_E) {

					push(@New_Paths, $Path." ".$i." ");

					$Path_End = 0;

					$Tail = 0;

				}else{

					print "unknown relationship\n";

					exit(0);

				}

			}



			if($Path_End == 1) {

				push(@New_Paths, $Path." End");

			}

		}

		@Paths = @New_Paths;

#		print "\n";

	}



	# ----------- find the path that has the highest score -------------

	$Most_Reliable_Path = "";

	$Highest_Score = 0;

	foreach(@Paths) {

		$Path = $_;

		@Content = split(/ +/, $Path);

		$Total_Score = 0;

		foreach(@Content) {

			next if(($_ eq "Begin")||($_ eq "End"));

			$Score = $ID_Score{$_};

			$Total_Score = $Total_Score + $Score;

		}

#		print "Final Path: $Path \t $Total_Score\n";



		if($Total_Score > $Highest_Score) {

			$Most_Reliable_Path = $Path;

			$Highest_Score = $Total_Score;

		}

	}

#	print "\n most likely: $Most_Reliable_Path\n";



#	exit(0);



	return($Most_Reliable_Path);

}



#--------------------------------------------------------------

sub score_caculater {

	my($Q_Match, $S_Match) = @_;

	my(@Q_BPs, @S_BPs, $Q_BP, $S_BP, $i, $Score, $Len);

	$Score = 0;

	$Len = length($Q_Match);

	@Q_BPs = split(//, $Q_Match);

	@S_BPs = split(//, $S_Match);

	for($i = 0; $i < $Len; $i ++) {

		$Q_BP = $Q_BPs[$i];

		$S_BP = $S_BPs[$i];

		if(($Q_BP eq "-")||($S_BP eq "-")) {

			$Score += -10;

		}elsif(($Q_BP eq "*")||($S_BP eq "*")) {

			$Score += -10;

		}elsif(($Q_BP eq "X")||($S_BP eq "X")) {

			$Score += -10;

		}else{

			$Score += $AA_Score{$Q_BP." ".$S_BP};

		}

	}

	return($Score);

}

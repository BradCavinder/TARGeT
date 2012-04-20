#!/usr/local/bin/perl -w
# Yujun Han
# The difference from version 2 to 1 is 2 can detect frameshift. 
# Note the input file is also different for the additional frameshift information
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:w:m:M:s:P:n:o:h:");

$Input    = defined $opt_i ? $opt_i : "";
$User_W   = defined $opt_w ? $opt_w : 600;
$Matrix   = defined $opt_m ? $opt_m : "/Library/WebServer/Tools/BLOSUM62.txt";
$Min_Draw = defined $opt_M ? $opt_M : 15;
$Single   = defined $opt_s ? $opt_s : "none";
$Protein  = defined $opt_P ? $opt_P : 1;
$Max_Num  = defined $opt_n ? $opt_n : 100;
$Output   = defined $opt_o ? $opt_o : "result";
$Help     = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
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
}
# ----------------------------------------------------
$Max_Width = 0;
$Sub_Figure_Num = 0;
$Line_Num = 0;
open(IF, $Input)||die"$!\n";
while(<IF>) {
	chomp;
	$Line_Num ++;
	next if($Line_Num > $Max_Num);
	@All_Infor = split(/  /, $_);

	@Query_Begins = ();
	@Query_Ends = ();
	@Query_Matches = ();
	@Sbjct_Begins = ();
	@Sbjct_Ends = ();
	@Sbjct_Matches = ();

	for($i = 0; $i < @All_Infor; $i ++) {
		next if(!$All_Infor[$i]);
		if($i == 0) {
#			print "$All_Infor[$i]\n";
			($Name, $Query, $Query_Len, $Database, $Sbject) = split(/ /, $All_Infor[$i]);
			if($Single ne "none") {
				last if($Name ne $Single);   # ------------------- draw single copy, if required --------------------
			}
			$Sub_Figure_Num ++;

		}else{
			($Q_Begin, $Q_End, $S_Begin, $S_End, $Frame, $Q_Match, $S_Match) = split(/ /, $All_Infor[$i]);
#			print "$Q_Begin, $Q_End, $S_Begin, $S_End\n";

			push(@Query_Begins, $Q_Begin);
			push(@Query_Ends, $Q_End);
			push(@Query_Matches, $Q_Match);
			push(@Sbjct_Begins, $S_Begin);
			push(@Sbjct_Ends, $S_End);
			push(@Sbjct_Matches, $S_Match);
		}
	}
	if($Single ne "none") {
		next if($Name ne $Single);   # ------------------- draw single copy, if required --------------------
	}
	$Q_Start = $Query_Begins[0];
	$S_Start = $Sbjct_Begins[0];
		
	$Width = $Q_Start + abs(int(($S_End - $S_Start)/3)) + ($Query_Len - $Q_End) if($Protein == 1);
	$Width = $Q_Start + abs(int($S_End - $S_Start)) + ($Query_Len - $Q_End) if($Protein == 0);

	$Max_Width = $Width > $Max_Width ? $Width : $Max_Width;
}

#--------------------------------------------------- Draw SVG
$Full_Name = $Output.".svg";
#$X_Press = int($Max_Width/1800) + 1 if($Protein == 1);
$X_Press = int($Max_Width/$User_W) + 1;

$SVG_Width = $User_W + 280;
$Sub_Height = 50;
$SVG_Height = $Sub_Figure_Num * $Sub_Height + 45;

open(OF, ">$Full_Name")||die"$!\n";
print OF <<"    _SVG_";
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="$SVG_Width" height="$SVG_Height" xmlns="http://www.w3.org/2000/svg">
    _SVG_

#-------------------------- draw ruler ------------------
$Sub_Figure = 0;
$Start_X = 200;
$Start_Y = 20;
SVG_liner($Start_X, $Start_Y, $Start_X + $User_W, $Start_Y, 1, "stroke:black;");
for($i = 0; $i <= $User_W; $i ++) {
	SVG_liner($Start_X + $i, $Start_Y, $Start_X + $i, $Start_Y +5, 1, "stroke:black;")if($i%10 == 0);
	SVG_liner($Start_X + $i, $Start_Y, $Start_X + $i, $Start_Y +10, 1, "stroke:black;") if($i%50 == 0);
	if($Protein == 1) {
		SVG_texter($Start_X + $i, $Start_Y - 10, 10, "black", $i*$X_Press*3) if($i%100 == 0);
	}else{
		SVG_texter($Start_X + $i, $Start_Y - 10, 10, "black", $i*$X_Press) if($i%100 == 0);
	}
}
SVG_texter($Start_X + $i, $Start_Y, 12, "black", "(Base Pairs)");

open(IF, $Input)||die"$!\n";
@Copies = ();
while(<IF>) {
	push(@Copies, $_);
}


$Finish = 1;
while($Finish == 0) {
	$Finish = 1;
	$Max_Length = 0;
	for($i = 0; $i < @Copies; $i ++) {
		next if($Copies[$i] eq "none");
		@All_Infor = split(/  /, $Copies[$i]);
	
		@Copy_Loc = ();

		for($j = 1; $j < @All_Infor - 1; $j ++) {
			($Q_Begin, $Q_End, $S_Begin, $S_End, $Frame, $Q_Match, $S_Match) = split(/ /, $All_Infor[$j]);
			push(@Copy_Loc, $S_Begin);
			push(@Copy_Loc, $S_End);
		}
		$Copy_Len = abs($Copy_Loc[0] - $Copy_Loc[$#Copy_Loc]) + 1;
		if($Copy_Len > $Max_Length) {
			$Max_Length = $Copy_Len;
			$Max_ID     = $i;
		}
		$Finish = 0;
	}
	#push(@Copy_Big2Small, $Copies[$Max_ID]);
	$Copies[$Max_ID] = "none";
}

foreach(@Copies) {
	@All_Infor = split(/  /, $_);

	@Query_Begins = ();
	@Query_Ends = ();
	@Query_Matches = ();
	@Sbjct_Begins = ();
	@Sbjct_Ends = ();
	@Sbjct_Matches = ();
	@Frames = ();
	for($i = 0; $i < @All_Infor - 1; $i ++) {
		if($i == 0) {
			($Name, $Query, $Query_Len, $Database, $Sbject) = split(/ /, $All_Infor[$i]);
			if($Single ne "none") {
				last if($Name ne $Single);   # ------------------- draw single copy, if required --------------------
			}
		}else{
			($Q_Begin, $Q_End, $S_Begin, $S_End, $Frame, $Q_Match, $S_Match) = split(/ /, $All_Infor[$i]);
			push(@Query_Begins, $Q_Begin);
			push(@Query_Ends, $Q_End);
			push(@Query_Matches, $Q_Match);
			push(@Sbjct_Begins, $S_Begin);
			push(@Sbjct_Ends, $S_End);
			push(@Sbjct_Matches, $S_Match);
			push(@Frames, $Frame);
		}
	}
	if($Single ne "none") {
		next if($Name ne $Single);   # ------------------- draw single copy, if required --------------------
	}

	$Q_Start = $Query_Begins[0];
	$S_Start = $Sbjct_Begins[0];
		

	$Q_Head = $Query_Begins[0];
	$S_Head = $Sbjct_Begins[0];

	$Start_Y = 55 + $Sub_Figure * $Sub_Height;

	$Dir = $S_Begin > $S_End ? "-" : "+";
	SVG_texter(10, $Start_Y, 14, "red", $Name." (".$Dir.")");	
	SVG_texter(10, $Start_Y +20, 12, "black", $Sbject." (BPs)");	
	SVG_liner(10, $Start_Y +30, 200, $Start_Y + 30, 1, "stroke:black;");
	
	if($Protein == 1) {
		$Copy_End = $Q_Start + abs(int(($S_End - $S_Start)/3)) + ($Query_Len - $Q_End);
	}else{
		$Copy_End = $Q_Start + abs(int($S_End - $S_Start)) + ($Query_Len - $Q_End);
	}
	SVG_liner($Start_X, $Start_Y + 3, $Start_X + $Copy_End/$X_Press, $Start_Y + 3, 2, "stroke:grey;");

	$S_Gap = 0;
	for($i = 0; $i < @Query_Begins; $i ++) {
		$Q_B = $Query_Begins[$i];
		$Q_E = $Query_Ends[$i];
		$S_B = $Sbjct_Begins[$i];
		$S_E = $Sbjct_Ends[$i];

		$Q_M = $Query_Matches[$i];
		$S_M = $Sbjct_Matches[$i];

		@Q_AAs = split(//, $Q_M);
		@S_AAs = split(//, $S_M);

		$S_B_Real = $S_B;
		$S_E_Real = $S_E;
		
		$S_B = $Dir eq "-" ? $S_B - $S_Gap * 3 : $S_B + $S_Gap * 3 if($Protein == 1);
		$S_B = $Dir eq "-" ? $S_B - $S_Gap : $S_B + $S_Gap if($Protein == 0);

		for($j = 0; $j < @Q_AAs; $j ++) {
			$S_Gap ++ if($S_AAs[$j] eq "-");
		}
		$S_E = $Dir eq "-" ? $S_E - $S_Gap * 3 : $S_E + $S_Gap * 3 if($Protein == 1);
		$S_E = $Dir eq "-" ? $S_E - $S_Gap : $S_E + $S_Gap if($Protein == 0);

		# ----------------------- draw the miss matched begin part ----------------------

		if($i == 0) {
			SVG_liner($Start_X, $Start_Y + 3, $Start_X + $Q_Head/$X_Press, $Start_Y + 3, 2, "stroke:blue;");
#			SVG_texter($Start_X-10, $Start_Y + 7, 10, "red", 1);	
		}

		# -------------- draw the matche -------------------------------
		SVG_rectangle($Start_X + ($Q_Head + abs(int(($S_B - $S_Head)/3)))/$X_Press, $Start_Y, abs(int(($S_E - $S_B)/3))/$X_Press, 7, "grey", "grey")  if($Protein == 1);
		SVG_rectangle($Start_X + ($Q_Head + abs(int($S_B - $S_Head)))/$X_Press, $Start_Y, abs(int($S_E - $S_B))/$X_Press, 7, "grey", "grey")  if($Protein == 0);
		
		for($j = 0; $j < @Q_AAs; $j ++) {
			$X = $Q_Head + abs(int(($S_B - $S_Head)/3)) + $j  if($Protein == 1);
			$X = $Q_Head + abs(int($S_B - $S_Head)) + $j  if($Protein == 0);
			if($Protein == 1) {
				if(($S_AAs[$j] eq "X")||($S_AAs[$j] eq "*")) {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y-2, 1, "stroke:red;");
					SVG_circle($Start_X + $X/$X_Press, $Start_Y - 4, 2, "red", 1, "red");
				}elsif($S_AAs[$j] eq "!") {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y-2, 1, "stroke:blue;");
					SVG_circle($Start_X + $X/$X_Press, $Start_Y - 4, 2, "blue", 1, "blue");
				}elsif($S_AAs[$j] eq "-") {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y + 5, 1, "stroke:white;");
				}elsif($Q_AAs[$j] eq "-") {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 3, $Start_X + $X/$X_Press, $Start_Y, 1, "stroke:white;");
				}elsif($S_AAs[$j] eq $Q_AAs[$j]) {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y, 1, "stroke:black;");
				}else{
					$Red   = 170 - $AA_Score{$Q_AAs[$j]." ".$S_AAs[$j]} * 10;
					$Green = 170 - $AA_Score{$Q_AAs[$j]." ".$S_AAs[$j]} * 10;
					$Blue  = 170 - $AA_Score{$Q_AAs[$j]." ".$S_AAs[$j]} * 10;
					$Red = $Red < 0 ? 0 : $Red;
					$Red = $Red > 255 ? 255 : $Red;
					$Green = $Green < 0 ? 0 : $Green;
					$Green = $Green > 255 ? 255 : $Green;
					$Blue = $Blue < 0 ? 0 : $Blue;
					$Blue = $Blue > 255 ? 255 : $Blue;
					$Color = "stroke:rgb($Red, $Green, $Blue);";
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y, 1, $Color);				
				}
			}else{
				if($S_AAs[$j] eq $Q_AAs[$j]) {
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y, 1, "stroke:black;");
				}
			}
		}

		# ---------------- mark the start locations of query and sbject -------------------
		$X = $Q_Head + abs(int(($S_B - $S_Head)/3))  if($Protein == 1);
		$X = $Q_Head + abs(int($S_B - $S_Head))  if($Protein == 0);
#		$D = $Protein == 1 ? abs(($Sbjct_Ends[$i] - $Sbjct_Begins[$i])/3)/$X_Press : abs($Sbjct_Ends[$i] - $Sbjct_Begins[$i])/$X_Press;
		
		if($i == 0) {
			if($Q_B > 1) {
				SVG_texter($Start_X + $X/$X_Press, $Start_Y - 8, 9, "red", $Q_B);
			}else{
				SVG_texter($Start_X + $X/$X_Press, $Start_Y - 8, 9, "black", $Q_B);
			}
		}else{
			$Q_E_Pre = $Query_Ends[$i - 1];
			if($Q_B - $Q_E_Pre > 1) {
				SVG_texter($Start_X + $X/$X_Press, $Start_Y - 8, 9, "red", $Q_B) if($Q_E - $Q_B >= $Min_Draw);
			}else{
				SVG_texter($Start_X + $X/$X_Press, $Start_Y - 8, 9, "black", $Q_B) if($Q_E - $Q_B >= $Min_Draw);;
			}
		}
		SVG_liner($Start_X + $X/$X_Press, $Start_Y, $Start_X + $X/$X_Press, $Start_Y-6, 1, "stroke:grey;") if($Q_E - $Q_B >= $Min_Draw);;

		SVG_texter($Start_X - 50, $Start_Y+7, 9, "black", $S_B_Real) if($i == 0);	

		# ---------------- mark the stop locations of query and sbject --------------------
		$X = $Q_Head+abs(int(($S_B-$S_Head)/3))+abs(int(($S_E-$S_B)/3)) if($Protein == 1);
		$X = $Q_Head+abs(int($S_B-$S_Head))+abs(int($S_E-$S_B)) if($Protein == 0);
		if($i == @Query_Begins - 1) {      # ----------------------- last exon -------------------------
			if($Q_E + 1 < $Query_Len) {
				SVG_texter($Start_X+$X/$X_Press, $Start_Y+22, 9, "red", $Q_E);
			}else{
				SVG_texter($Start_X+$X/$X_Press, $Start_Y+22, 9, "black", $Q_E);
			}
			SVG_liner($Start_X+$X/$X_Press, $Start_Y, $Start_X+$X/$X_Press, $Start_Y+15, 1, "stroke:grey;");

			SVG_texter($Start_X+$Copy_End/$X_Press + 10, $Start_Y+7, 9, "black", $S_E_Real);
		}else{
			$Q_B_Next = $Query_Begins[$i + 1];
			
			if($Q_E + 1 < $Q_B_Next) {
				SVG_texter($Start_X+$X/$X_Press, $Start_Y+22, 9, "red", $Q_E) if($Q_E - $Q_B >= $Min_Draw);;
			}else{
				SVG_texter($Start_X+$X/$X_Press, $Start_Y+22, 9, "black", $Q_E) if($Q_E - $Q_B >= $Min_Draw);;
			}
			SVG_liner($Start_X+$X/$X_Press, $Start_Y, $Start_X+$X/$X_Press, $Start_Y+12, 1, "stroke:grey;") if($Q_E - $Q_B >= $Min_Draw);;
		}

		# ---------------- detect putative frameshift --------------------------
		if(($i + 1 < @Query_Begins)&&($Protein == 1)) { 
			$Q_E = $Query_Ends[$i];
			$S_E = $Sbjct_Ends[$i];
			$Frame = $Frames[$i];

			#$Q_B_N = $Query_Begins[$i+1];
			$S_B_N = $Sbjct_Begins[$i+1];
			$Frame_N = $Frames[$i+1];
			
			if(abs($S_B_N - $S_E) < 5) {
				if($Frame != $Frame_N) {
					$X = $Q_Head + abs(int(($S_B - $S_Head)/3)) + $j - 1;
					SVG_liner($Start_X + $X/$X_Press, $Start_Y + 7, $Start_X + $X/$X_Press, $Start_Y-2, 1, "stroke:blue;");
					SVG_circle($Start_X + $X/$X_Press, $Start_Y - 4, 2, "blue", 1, "blue");
				}
			}
		}
	}

	# ----------------------- draw the miss matched end part
	$X = $Q_Head + abs(int(($S_E - $S_Head)/3)) if($Protein == 1);
	$X = $Q_Head + abs(int($S_E - $S_Head)) if($Protein == 0);
	SVG_liner($Start_X + $X/$X_Press, $Start_Y + 3, $Start_X + ($X+($Query_Len-$Q_E))/$X_Press, $Start_Y + 3, 2, "stroke:blue;");
#	SVG_texter($Start_X + ($X+($Query_Len - $Q_E)+10)/$X_Press, $Start_Y + 10, 10, "red", $Query_Len);
	
	$Sub_Figure ++;
}
print(OF "</svg>\n");
close(IF);
close(OF);

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :tblastn_copy_finder.pl <options> <specification> <default>

	\$Input   = defined \$opt_i ? \$opt_i : "";
	\$User_W  = defined \$opt_w ? \$opt_w : 600;
	\$Matrix  = defined \$opt_m ? \$opt_m : "/Library/WebServer/Tools/BLOSUM62.txt";
	\$Single  = defined \$opt_s ? \$opt_s : "none";
	\$Protein = defined \$opt_P ? \$opt_P : 1;
	\$Output  = defined \$opt_o ? \$opt_o : "result";
	\$Help    = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}

#-----------------------------------------------------
sub SVG_liner {
	my($X1, $Y1, $X2, $Y2, $Width, $Color) = @_;
	print(OF '<line x1 = "');
	print(OF $X1);
	print(OF '" y1 = "');
	print(OF $Y1);
	print(OF '" x2 = "');
	print(OF $X2);
	print(OF '" y2 = "');
	print(OF $Y2);
	print(OF '" style = "');
	print(OF $Color);
	print(OF '" stroke-width = "');
	print(OF $Width);
	print(OF '"/>');
	print(OF "\n");
}

#-----------------------------------------------------
sub SVG_texter {
	my($X, $Y, $Font_Size, $Color, $Text) = @_;
	print(OF '<text x = "');
	print(OF $X);
	print(OF '" y = "');
	print(OF $Y);
	print(OF '" font-family = "');
	print(OF Arial);
	print(OF '" font-size = "');
	print(OF $Font_Size);
	print(OF '" fill = "');
	print(OF $Color);
	print(OF '">');
	print(OF $Text);
	print(OF '</text>');
	print(OF "\n");
}

#-----------------------------------------------------
sub SVG_rectangle {
	my($X, $Y, $W, $H, $Fill, $Stroke) = @_;
	print(OF '<rect x="');
	print(OF $X);
	print(OF '" y="');
	print(OF $Y);
	print(OF '" width="');
	print(OF $W);
	print(OF '" height="');
	print(OF $H);
	print(OF '" style = "');
	print(OF "fill:$Fill;");
	print(OF "stroke:$Stroke");
	print(OF '"/>');
	print(OF "\n");
}

#-----------------------------------------------------
sub SVG_circle {
	my($X, $Y, $R, $Stroke, $Stroke_W, $Fill) = @_;
	print(OF '<circle cx="');
	print(OF $X);
	print(OF '" cy="');
	print(OF $Y);
	print(OF '" r="');
	print(OF $R);
	print(OF '" stroke = "');
	print(OF $Stroke);
	print(OF '" stroke-width = "');
	print(OF $Stroke_W);
	print(OF '" fill = "');
	print(OF $Fill);
	print(OF '"/>');
	print(OF "\n");	
}

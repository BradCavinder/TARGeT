#!/usr/bin/env perl
use warnings;
# Yujun Han
#-----------------------------------------------------
use Getopt::Std;
#-----------------------------------------------------
getopts("i:W:H:o:h:");

$Input  = defined $opt_i ? $opt_i : "";
$Width  = defined $opt_W ? $opt_W : 560;
$Height = defined $opt_H ? $opt_H : 200;
$Output = defined $opt_o ? $opt_o : "result";
$Help   = defined $opt_h ? $opt_h : "";

usuage() if((!$Input)||($Help));
#-----------------------------------------------------
open(IF, $Input)||die"$!\n";
$Head = 1;
$Q_Match = "";
$S_Match = "";
$DNA = 0;

while(<IF>) {
	chomp;
	$Line = $_;
	if(($Line=~ /^BLASTN/)||($Line=~ /^TBLASTX/)) {
		$DNA = 1;
	}
	if($Line =~ /Query= (\S+)/) {
		$Query = $1;
		$No_Hits = 0;
	}
	if($Line =~ /No hits found/) {
		$No_Hits = 1;
	}
	if($Line =~ /\((\S+) letters\)/) {
		$Query_Len = $1;
		$Query_Len =~ s/\,//;
#		$Query_Length{$Query} = $Query_Len;
	}
	if($Line =~ /Database: (\S+)/) {
		$Database = $1;
	}
	if($Line =~ /^>(\S+)/) {
		if($Head == 0) {
			Sbject_Info_Loader();
			$Head = 1;
		}
		$Sbject = $1;
		$Head = 1;
	}
	if($Line =~ /Length = (\d+)/) {
#		$Sbject_Len = $1;
	}
	if($Line =~ /Score =/) {
		if($Head == 0) {
			Sbject_Info_Loader();
		}else{
			$Head = 0;
		}
		$Match_Begin = 1;
		$Q_Match = "";
		$S_Match = "";
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
}

Sbject_Info_Loader();
close(IF);

$SVG_ID = 1;
foreach(keys(%Query_Hits)) {
	$Query = $_;
	$Query_Len = $Query_Info{$_};
	@Query_BPs  = ();
	@Query_BPs_Same  = ();
	for($i = 1; $i < $Query_Len - 1; $i ++) {
		$Query_BPs[$i]  = 0;
		$Query_BPs_Same[$i]  = 0;
	}
	@Contents = split(/  /, $Query_Hits{$_});
	foreach(@Contents) {
		($Q_Begin, $Q_End, $Q_Match, $S_Match, $Direction) = split(/ /, $Blast_Hit{$_});
		for($i = $Q_Begin; $i < $Q_End; $i ++) {
			$Query_BPs[$i] ++;
		}
		@Q_Letters = split(//, $Q_Match);
		@S_Letters = split(//, $S_Match);
		$Gap_Additive = 0;
		for($i = $Q_Begin; $i - $Gap_Additive < $Q_End; $i ++) {
			if($Q_Letters[$i - $Q_Begin] eq "-") {
				$Gap_Additive ++;
			}else{
				$Query_BPs_Same[$i - $Gap_Additive]  ++ if($Q_Letters[$i - $Q_Begin] eq $S_Letters[$i - $Q_Begin]);
			}
		}
	}

	#___________________ Draw SVG
	$Full_Name = $Output."$SVG_ID".".svg";
	open(OF, ">$Full_Name")||die"$!\n";
	print OF <<"    _SVG_";
<?xml version="1.0" standalone="no"?>
<!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.0//EN" "http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd">
<svg width="$Width" height="$Height" xmlns="http://www.w3.org/2000/svg">
    _SVG_
	
	SVG_liner(10, 10, 550, 10, 1, "black");		#------------- draw the border
	SVG_liner(10, 10, 10, 190, 1, "black");
	SVG_liner(550, 10, 550, 190, 1, "black");
	SVG_liner(10, 190, 550, 190, 1, "black");

	$Red	= 200; 
	$Green	= 200;
	$Blue	= 200;

	$X_Times = int($Query_Len / 400) + 1;
	$Draw_len = int($Query_Len/$X_Times);

	@Draw_Data = ();
	@Same_Data = ();
	if($X_Times > 1) {						#------------------- caculate and press the data
		for($i = 1; $i < $Draw_len - 1; $i ++) {
			$Sum  = 0;
			$Same = 0;
			last if($i * $X_Times + $X_Times == $Query_Len - 1);
			for($j = $i * $X_Times; $j < $i * $X_Times + $X_Times; $j ++) {
				$Sum  += $Query_BPs[$j];
				$Same += $Query_BPs_Same[$j];
			}
			$Draw_Data[$i]  = int($Sum/$X_Times) + 1;
			$Same_Data[$i]  = int($Same/$X_Times) + 1;
		}
	}else{
		@Draw_Data = @Query_BPs;
		@Same_Data = @Query_BPs_Same;
	}

	$Max_Copy = 0;							# ------------------ find the highest copies
	for($i = 1; $i < $Draw_len - 1; $i ++) {
		$Max_Copy = $Draw_Data[$i] > $Max_Copy ? $Draw_Data[$i] : $Max_Copy;
	}
	$Y_Times = int($Max_Copy / 100) + 1;

	if(($Max_Copy < 50)&&($Max_Copy > 25)) {
		$Y_Times = 0.5;
	}elsif(($Max_Copy <= 25)&&($Max_Copy > 10)) {
		$Y_Times = 0.25;
	}elsif($Max_Copy <= 10) {
		$Y_Times = 0.1;
	}

	SVG_texter(140, 30, 15, "blue", "Query: ".$Query." (".$Query_Len." BPs)") if($DNA == 1);
	SVG_texter(140, 30, 15, "blue", "Query: ".$Query." (".$Query_Len." AAs)") if($DNA == 0);
	SVG_liner(60, 150, 60 + $Draw_len, 150, 1, "black");	#---------------- draw the legends
	SVG_liner(60, 50, 60, 150, 1, "black");
	SVG_texter(45, 150, 10, "black", 0);
	SVG_texter(240, 180, 14, "red", "BLAST Datbase: ".$Database);

	SVG_liner(60, 50, 68, 50, 1, "black");
	SVG_liner(60, 75, 68, 75, 1, "black");
	SVG_liner(60, 100, 68, 100, 1, "black");
	SVG_liner(60, 125, 68, 125, 1, "black");

	SVG_texter(30, 30, 12, "blue", "(Hits)");
	SVG_texter(30, 50, 12, "black", int(100 * $Y_Times));
	SVG_texter(30, 75, 12, "black", int(75 * $Y_Times));
	SVG_texter(30, 100, 12, "black", int(50 * $Y_Times));
	SVG_texter(30, 125, 12, "black", int(25 * $Y_Times));

	for($i = 1; $i < $Draw_len - 1; $i ++) {			#------------- Draw data
		next if($Draw_Data[$i] == 0);
		
		$X1 = 60 + $i;
		$Y1 = 149;
		$X2 = 60 + $i;
		$Y2 = 150 - $Draw_Data[$i]/$Y_Times;

		$Conserve_Ratio = $Same_Data[$i] / $Draw_Data[$i];
		$Red   = 200 - $Conserve_Ratio * 100 * 1.5;
		$Green = 200 - $Conserve_Ratio * 100 * 1.5;
		$Blue  = 200 - $Conserve_Ratio * 100 * 1.5;
		$Color = "stroke:rgb($Red, $Green, $Blue);";
		SVG_liner($X1, $Y1, $X2, $Y2, 1, $Color);
	}
	SVG_texter(70 + $i, 155, 12, "blue", "(Length, BPs)") if($DNA == 1);
	SVG_texter(70 + $i, 155, 12, "blue", "(Length, AAs)") if($DNA != 1);

	for($i = 1; $i < 10; $i ++) {			#------------------- draw X ruler
		$X_Mark = $i * 50 * $X_Times;
		next if($X_Mark > $Query_Len);
		SVG_liner(60 + $i * 50, 150, 60 + $i * 50, 145, 1, "black");
		SVG_texter(60 + $i * 50 - 10, 160, 12, "black", $X_Mark);
	}

	for($i = 0; $i < 100; $i ++) {			#------------------- draw conserve marker
		$Red   = 200 - $i * 1.5;
		$Green = 200 - $i * 1.5;
		$Blue  = 200 - $i * 1.5;
		SVG_liner(500, 150 - $i, 510, 150 - $i, 1, "stroke:rgb($Red, $Green, $Blue);");
	}
	SVG_texter(500, 40, 10, "black", "100%");
	SVG_texter(500, 160, 10, "grey", "0%");
	SVG_texter(480, 170, 12, "black", "Similarity");
	print(OF "</svg>\n");
	close(OF);
	$SVG_ID ++;
}

#-----------------------------------------------------
#-----------------------------------------------------
sub usuage {
    print "\nHi, need some help?\n";
    print STDERR <<"    _EOT_";

    Usage :tblastn_copy_finder.pl <options> <specification> <default>

	\$Input  = defined \$opt_i ? \$opt_i : "";
	\$Width  = defined \$opt_W ? \$opt_W : 620;
	\$Height = defined \$opt_H ? \$opt_H : 300;
	\$Output = defined \$opt_o ? \$opt_o : "result";
	\$Help   = defined \$opt_h ? \$opt_h : "";

    _EOT_
    exit(1)
}

#-------------------------------------------------------
sub Sbject_Info_Loader {
	if($No_Hits == 0) {
		$Direction  = $S_Begin > $S_End ? -1 : 1;
		$Query_Info{$Query}  = $Query_Len;
		$Query_Hits{$Query} .= $Sbject." ".$S_Begin."  ";
		$Blast_Hit{$Sbject." ".$S_Begin} = $Q_Begin." ".$Q_End." ".$Q_Match." ".$S_Match." ".$Direction;
	}
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


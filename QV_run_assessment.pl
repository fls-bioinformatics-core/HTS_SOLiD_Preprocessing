#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use Getopt::Long;

# Written by Ariella Sasson, Rutgers University, version 1, Nov 1, 2009 

# This script provides an analysis of the SOLiD QV file.  It samples 20% of the data available to speed up runtime and gives consistant results with a full analysis of the data.


$0 =~ s-.*/--g;
my ($inputF,$tLen,$negQV,$outputFile);
#input output option
GetOptions('i|input=s'        => \$inputF,
           't|trunc:i'        => \$tLen, 
           'm|misCall_sig:s'  => \$negQV,
           'o|output=s'       => \$outputFile,
           );
my @assign;
my $acnt;
@assign=&AssignInputs();
$acnt=@assign;

if ($acnt > 0){
    my $trash = join('',@assign);
    print "\n$trash";
    usage();
}

{
    open (iFP,$inputF)  or die "Can not open $inputF.\n";

    my $output=">".$outputFile;
    open (oFP,$output) or die "Can not open $output.\n";

    my $counter=0;
    my $Scounter=0;
    my @matrix;
    my $pic=0;
    my @rmNeg;
    my @avg20;
    my @p1_eoff;
    my @p3_eoff;
    my @p5_eoff;
    my @p1_e3;
    my @p1_e5;
    my @p5_e0;
    my @holdID;

#change these 2 loops to increase the matrix size

    for my $i ( 0 .. 55 ) {
	for my $j ( 0 .. 50 ) {
	    $matrix[$i][$j]=0;
	    $rmNeg[$i][$j]=0;
	    $avg20[$i][$j]=0;
	    $p1_eoff[$i][$j]=0;
	    $p3_eoff[$i][$j]=0;
	    $p5_eoff[$i][$j]=0;
	    $p1_e3[$i][$j]=0;
	    $p1_e5[$i][$j]=0;
	    $p5_e0[$i][$j]=0;
	}
    }

#reading line by line.  The bead numbers are ignored.
    while (<iFP>) {
	chomp;
	$b=">";
	if (/^($b)|^\#/) {
	    next;
	}
	my $qv = $_;
	$counter++;
	$pic++;

	if ($pic!=5){
	    next;
	}
	else{

	    $pic=0;
	    $Scounter++;

	    my @F3_arr=split (/\s/,$qv);
	    my $test=@F3_arr;

	    if ($tLen == -20){
		$tLen=$test;
	    }
	    elsif ($tLen > $test) {
		$tLen=$test;
            }


# Set up for p & e analysis
            for my $m ( 0 .. 2) {
                for my $n ( 0 .. $test-1 ) {
                    $holdID[$m][$n]=0;
                }
            }

# Full QV analysis No truncation 
	    for (0..$test-1){
		my $a = $F3_arr[$_];
		$a = $a + 1;
		$matrix[$_][$a]++;
		
		if ($a >= 26){
		    $holdID[0][$_]=1;
		}
		if ($a <= 11){
		    $holdID[1][$_]=1;
		}
		if ($a == 0) {
		    $holdID[2][$_]=1;
		}
	    }

# Average QV analysis min Score 20
	    my $sum = 0;
	    my $negInd=0;

	    for (0 .. $tLen-1){
		$sum=$sum + $F3_arr[$_];
		if (($F3_arr[$_] < 0) && ($negQV eq 'y')){
		    $negInd=1;
		}
	    }
	    my $avg_sc = $sum/$tLen;
	    
	    if ((($negQV eq 'y') && ($negInd == 0)) || ($negQV eq 'n')){
		if ($avg_sc >= 20) {
		    for (0..$tLen-1){
			my $a = $F3_arr[$_];
			$a = $a + 1;
			$avg20[$_][$a]++;
		    }
		}
	    }

# p & e analysis
	    my $psum = 0;
	    my $esum = 0;
	    my $negsum = 0;
	    
	    my $holdL = 9;
	    if ($tLen -1  < 9){
		$holdL = $tLen -1;
	    }

	    for (0 .. $holdL){
		$psum = $psum + $holdID[0][$_];
	    }  
	    for (0 .. $tLen-1){
		$esum = $esum + $holdID[1][$_];
		$negsum = $negsum + $holdID[2][$_];
	    }

# Removal of miscalls
	    if ($negsum ==0){
		for (0..$tLen-1){
                    my $a = $F3_arr[$_];
                    $a = $a + 1;
                    $rmNeg[$_][$a]++;
                }
	    }

# p1 e off
	    if ($negQV eq 'y'){
		if (($negsum ==0) && ($psum >=1)){
		    for (0..$tLen-1){
			my $a = $F3_arr[$_];
			$a = $a + 1;
			$p1_eoff[$_][$a]++;
		    }
		}
	    } else {
		if (($psum >=1)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p1_eoff[$_][$a]++;
                    }
                }
	    }

# p3 e off
            if ($negQV eq 'y'){
                if (($negsum ==0) && ($psum >=3)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p3_eoff[$_][$a]++;
                    }
                }
            } else {
                if (($psum >=3)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p3_eoff[$_][$a]++;
                    }
                }
            }

# p5 e off
            if ($negQV eq 'y'){
                if (($negsum ==0) && ($psum >=5)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p5_eoff[$_][$a]++;
                    }
                }
            } else {
                if (($psum >=5)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p5_eoff[$_][$a]++;
                    }
                }
            }

#p1 e5
            if ($negQV eq 'y'){
                if (($negsum ==0) && ($psum >=1) && ($esum <= 5)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p1_e5[$_][$a]++;
                    }
                }
            } else {
                if (($psum >=1) && ($esum <= 5)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p1_e5[$_][$a]++;
                    }
                }
            }

#p1 e3
            if ($negQV eq 'y'){
                if (($negsum ==0) && ($psum >=1) && ($esum <= 3)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p1_e3[$_][$a]++;
		    }
                }
            } else {
                if (($psum >=1) && ($esum <= 3)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p1_e3[$_][$a]++;
                    }
                }
            }


#p5 e0
            if ($negQV eq 'y'){
                if (($negsum ==0) && ($psum >=5) && ($esum == 0)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p5_e0[$_][$a]++;
                    }
                }
            } else {
                if (($psum >=5) && ($esum == 0)){
                    for (0..$tLen-1){
                        my $a = $F3_arr[$_];
                        $a = $a + 1;
                        $p5_e0[$_][$a]++;
                    }
                }
            }
	}
    }

    close(iFP);

    print oFP "Total Count of reads - $counter \n";
    print oFP "Total COunt of reads sampled - $Scounter \n";
    print oFP "\n\n";

    print oFP "Full QV Analysis\n";
    for my $i ( 0 .. $#matrix) {
	print oFP "@{$matrix[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Removal of miscalls\n";
    for my $i ( 0 .. $#rmNeg) {
        print oFP "@{$rmNeg[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where the mean >=20 \n";
    for my $i ( 0 .. $#avg20) {
        print oFP "@{$avg20[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=1 e=off \n";
    for my $i ( 0 .. $#p1_eoff) {
        print oFP "@{$p1_eoff[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=3 e=off \n";
    for my $i ( 0 .. $#p3_eoff) {
        print oFP "@{$p3_eoff[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=5 e=off \n";
    for my $i ( 0 .. $#p5_eoff) {
        print oFP "@{$p5_eoff[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=1 e=5 \n";
    for my $i ( 0 .. $#p1_e5) {
        print oFP "@{$p1_e5[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=1 e=3 \n";
    for my $i ( 0 .. $#p1_e3) {
        print oFP "@{$p1_e3[$i]}\n";
    }
    print oFP "\n\n";

    print oFP "Reads where p=5 e=0 \n";
    for my $i ( 0 .. $#p5_e0) {
        print oFP "@{$p5_e0[$i]}\n";
    }
    print oFP "\n\n";

    close(oFP);
}

sub AssignInputs{
    my @err;

    if ((!defined($inputF)) || ($inputF eq '')) {
        push(@err, "Error 1: Input File not defined\n");
    }
#------------------------------
    if ((!defined($negQV)) || ($negQV eq '') || ($negQV eq 'off') || ($negQV eq 'n')||($negQV eq 'no')){
        $negQV = 'n';
        print "Removal of Reads containing a miscall - $negQV\n";
    }
    elsif (($negQV eq 'on')||($negQV eq 'y')||($negQV eq 'yes')) {
	$negQV = 'y';
	print "Removal of reads containing a miscall - $negQV\n";
    }
    else {
        push(@err,"Error 2: Removal of reads with negative quality scores can either be 'on' or 'off'.  The default is off.\n");
    }
#------------------------------
    if ((!defined($tLen)) || ($tLen eq '')){
        $tLen = -20;
	print "No truncation desired\n";
    }
    elsif ($tLen > 0) {
	$tLen = $tLen;
	print "Truncation length - $tLen\n";
    }
    else {
        push(@err,"Error 3: Improper truncation length.  Truncation must be turned on and the length must be greater than 0.\n");
    }
#--------------------------------------------------------------
    if ((!defined($outputFile)) || ($outputFile eq '')) {
	push(@err,"Error 3: Output File not defined\n");
    }
    
    return @err;
}

sub usage {
    
    print "\nusage: $0 \n[-i SOLiD QV file you want to analyze\n -t the length of desired read \n -m remove from analysis any reads that contain a miscall (.) \n -o output file name]\n\n";
    print "outputs a series of matrices that represent QV analyses for various filtering parameters.  The rows represent the positions and the columns represent the quality scores.  ";
    print "The numbers represented in the matrix are the counts of that score at that positions.\n\n";
    print "example: $0 -i Frag.qual -t 35 -m y -o Frag_QV_analysis.txt\n\n";
    exit;

}

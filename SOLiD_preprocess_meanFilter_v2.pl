#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;

use Getopt::Long;

# Written by Ariella Sasson, Rutgers University, version 2.0, June 2010, and is distributed as Open Source software under the GPLv3.0.
# v2 is SAET friendly.  This will allow to filter traditional SOLiD reads right off the machine as well as post SAET.  
# In addition, this version does not maintain the comments within the output.

# this script takes the SOLiD csfasta and QV files as input.  It then returns quality filtered csfasta & QV files.  
# For Mate pair analysis this script finds the mates from the F3 and R3 files (SOLiD output) 
# and places them in MP files and the remainder sequences, the orphans,
# are put into orphan files. The F3 and R3 files are kept independant as in Solid output.  
# But all the files are quality filtered so that those
# that meet the filter are put in the files with a label of T and those that don't have a label of U.  
# The output field requires the begining of the file name and the script fills in the rest.
# For Fragment analysis the most you will output is 4 files csfasta (T_F3.csfasta & U_F3.csfasta) 
# QV (QV_T_F3.qual & QV_U_F3.qual)
#####################

$0 =~ s-.*/--g;

#input output option
my ($input_t, $inputF, $inputQVF,$inputR,$inputQVR,$avg_sig,$avg_thresh,$negQV,$trunc_sig,$trunc_len,$qv_an,$outputFile,$outputQVFlag);

GetOptions('i|input_type:s'    => \$input_t,
           'f|f3=s'            => \$inputF,
           'g|f3QV:s'          => \$inputQVF,
           'r|r3:s'            => \$inputR,
           's|r3QV:s'          => \$inputQVR,
	   'p|mean_analysis:s' => \$avg_sig,
	   'q|mean_qv:i'       => \$avg_thresh, 
	   'n|neg_qv:s'        => \$negQV,
	   't|trunc:s'         => \$trunc_sig, 
	   'u|tr_len:i'        => \$trunc_len,
	   'a|qv_analysis:s'   => \$qv_an,
	   'o|output=s'        => \$outputFile,
	   'v|output_qv:s'     => \$outputQVFlag,
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

if ($input_t eq 'MP'){

    my @holdIA =split /_/,$inputF;
    my $holdI =pop @holdIA;
    @holdIA =split /\./,$holdI;
    my $qf=shift @holdIA;

    @holdIA =split /_/,$inputR;
    $holdI=pop @holdIA;
    @holdIA =split /\./,$holdI; #/
    my $qr=shift @holdIA;

# opens the input files
    open (iFP,$inputF) or die "Can not open $inputF.\n";
    open (iFPJ,$inputR) or die "Can not open $inputR.\n";
    open (iFPQV,$inputQVF) or die "Can not open $inputQVF.\n";
    open (iFPJQV,$inputQVR)or die "Can not open $inputQVR.\n";

# creates all the output files
    my $outputMPF3=">".$outputFile."_T_mp_F3.csfasta";
    open (oFPMPF3,$outputMPF3) or die "Can not open $outputMPF3.\n";
    my $outputMPR3=">".$outputFile."_T_mp_R3.csfasta";
    open (oFPMPR3,$outputMPR3) or die "Can not open $outputMPR3.\n";   
    
    my $outputOF3=">".$outputFile."_T_orphans_F3.csfasta";
    open (oFPOF3,$outputOF3) or die "Can not open $outputOF3.\n";
    my $outputOR3=">".$outputFile."_T_orphans_R3.csfasta";
    open (oFPOR3,$outputOR3) or die "Can not open $outputOR3.\n";
    
    if ($outputQVFlag eq 'y'){
	my $outputMPQVF3=">".$outputFile."_T_mp_F3_QV.qual";
	open (oFPMPQVF3,$outputMPQVF3) or die "Can not open $outputMPQVF3.\n";
	my $outputMPQVR3=">".$outputFile."_T_mp_R3_QV.qual";
	open (oFPMPQVR3,$outputMPQVR3) or die "Can not open $outputMPQVR3.\n";
	
	my $outputOQVF3=">".$outputFile."_T_orphans_F3_QV.qual";
	open (oFPOQVF3,$outputOQVF3) or die "Can not open $outputOQVF3.\n";
	my $outputOQVR3=">".$outputFile."_T_orphans_R3_QV.qual"; 
	open (oFPOQVR3,$outputOQVR3) or die "Can not open $outputOQVR3.\n";
    }
    
    my $outputMPUF3=">".$outputFile."_U_mp_F3.csfasta";
    open (oFPMPUF3,$outputMPUF3) or die "Can not open $outputMPUF3.\n";
    my $outputMPUR3=">".$outputFile."_U_mp_R3.csfasta";
    open (oFPMPUR3,$outputMPUR3) or die "Can not open $outputMPUR3.\n";
    
    my $outputOUF3=">".$outputFile."_U_orphans_F3.csfasta";
    open (oFPOUF3,$outputOUF3) or die "Can not open $outputOUF3.\n";
    my $outputOUR3=">".$outputFile."_U_orphans_R3.csfasta"; 
    open (oFPOUR3,$outputOUR3) or die "Can not open $outputOUR3.\n";  
    
    if ($outputQVFlag eq 'y'){
	my $outputMPQVUF3=">".$outputFile."_U_mp_F3_QV.qual";
	open (oFPMPQVUF3,$outputMPQVUF3) or die "Can not open $outputMPQVUF3.\n";
	my $outputMPQVUR3=">".$outputFile."_U_mp_R3_QV.qual";  
	open (oFPMPQVUR3,$outputMPQVUR3) or die "Can not open $outputMPQVUR3.\n";  
	
	my $outputOQVUF3=">".$outputFile."_U_orphans_F3_QV.qual";
	open (oFPOQVUF3,$outputOQVUF3) or die "Can not open $outputOQVUF3.\n";
	my $outputOQVUR3=">".$outputFile."_U_orphans_R3_QV.qual";
	open (oFPOQVUR3,$outputOQVUR3) or die "Can not open $outputOQVUR3.\n";  
    }

    if ($qv_an eq 'y'){
        my $qv_an_output=">".$outputFile."_QV_analysis.txt";
        open (oFPQVAN,$qv_an_output) or die "Can not open $qv_an_output.\n";
    }

###### QV ANALYSIS #######

    my $mqv_cntr_F3=0;
    my $mqv_cntr_R3=0;
    my $mqv_cntr_pass_F3=0;
    my $mqv_cntr_pass_R3=0;
    
    my @matrix_F3;
    my @p_matrix_F3;
    my @matrix_R3;
    my @p_matrix_R3;
    
    for my $l ( 0 .. 55 ) {
	for my $m ( 0 .. 50 ) {
	    $matrix_F3[$l][$m]=0;
	    $matrix_R3[$l][$m]=0;
	    $p_matrix_F3[$l][$m]=0;
	    $p_matrix_R3[$l][$m]=0;
	}
    }
##############################
 
    my $counter=0;
    my $lF = "";
    my $lR = "";
    my $lFQV="";
    my $lRQV="";

    my $holdTitleA="";
    my $holdTitleB="";


# takes all the comments at the beginning of the csfasta file and then outputs it to the outputs
    while (!eof(iFP) && !eof(iFPJ) ){
	if ($counter == 0 ){
	    $counter =1;
	    $lF=<iFP>;
	    $lFQV=<iFPQV>;
	    $lR=<iFPJ>;
	    $lRQV=<iFPJQV>;
	    chomp $lF;
	    chomp $lFQV;
	    chomp $lR;
	    chomp $lRQV;
	}

	$_=$lF;
	while (/^\#/ && !eof(iFP)) {
	   $_=<iFP>;
	   chomp;
	   $lF=$_;
	}

	$_=$lFQV;
	while (/^\#/ && !eof(iFPQV)) {
	    $_=<iFPQV>;
	    chomp;
	    $lFQV=$_;
        }

        $_=$lR;
	while (/^\#/ && !eof(iFPJ)) {
	    $_=<iFPJ>;
	    chomp;
	    $lR=$_;
	}

        $_=$lRQV;
        while (/^\#/ && !eof(iFPJQV)) {
            $_=<iFPJQV>;
            chomp;
            $lRQV=$_;
        }

	last;
    }

    my $id="";
    my $seq="";
    my $qv="";
    my $idR="";
    my $seqR="";
    my $qvR="";
    my @values;
    my @valuesR;
    my $idQV;
    my $idRQV;
    $counter=1;
    my $QVIDF=0;
    my $QVIDR=0;
    my $q="";
    my $qR="";

    print "Starting Analysis\n";
# this is the part that does the work.  Both files are open and we use the fact that SOLiD keeps the list in bead order to our advantage
# we basically start going through the list if F3 > R3 then the R3 sequence gets written into the orphan file, if F3< R3 then the F3 sequence
# get written to the orphan file.  If they match they get written to the MP files.  The threshold gets calculated inside each check.
# the threshold is a separate subroutine at the bottom of the script.

    while (!eof(iFP) && !eof(iFPJ) ) {
	if ($lF=~/_/ && $lR=~/_/) {
	    if ($id !~/_/){
		$_=$lF;
		tr/^>//d;
		@values =split /_/,$_;
		$q=pop @values;
		$_=join "_",@values;
		$id = $_;

                $_=$lFQV;
		tr/^>//d;
		my @valuesQV =split /_/,$_;
		my $qQV=pop @valuesQV;
		$_=join "_",@valuesQV;
		$idQV = $_;

	    }

	    if ($seq!~/[ACGT][.0123]+/){
		$_ = <iFP>;
		chomp;
		$seq=$_;
		$_=<iFPQV>;
       		chomp;
       		$qv=$_;
	    }

	    if ($idR !~/_/){
		$_=$lR;
		tr/^>//d;
		@valuesR =split /_/,$_;
		$qR=pop @valuesR;
		$_=join "_",@valuesR;
		$idR = $_;

                $_=$lRQV;
		tr/^>//d;
		my @valuesQVR =split /_/,$_;
		my $qQVR=pop @valuesQVR;
		$_=join "_",@valuesQVR;
		$idRQV = $_;

	    }

	    if ($seqR!~/[ACGT][.0123]+/){
		$_ = <iFPJ>;
		chomp;
		$seqR=$_;
		$_=<iFPJQV>;
		chomp;
		$qvR=$_;
	    }

# if the the sequence has a match
	    if ($values[0]==$valuesR[0] && $values[1]==$valuesR[1] && $values[2]==$valuesR[2]){
		if ($qv_an eq 'y'){

		    $mqv_cntr_F3++;
		    $mqv_cntr_R3++;
		    
		    my $q_ah_F3=$qv;
		    my @hold_a_f3=split (/\s/,$q_ah_F3);
		    my $test=@hold_a_f3;
		    for (0..$test-1){
			my $a = shift @hold_a_f3;
#			print "$a\n";
			$a = $a + 1;
			$matrix_F3[$_][$a]++;
		    }
		    
		    my $q_ah_r3=$qvR;
		    my @hold_a_r3=split (/\s/,$q_ah_r3);
		    $test=@hold_a_r3;
		    for (0..$test-1){
			my $a = shift @hold_a_r3;
			$a = $a + 1;
			$matrix_R3[$_][$a]++;
		    }
		}

		if ($trunc_sig eq 'y'){
		    my @trunc_res=truncation($seq,$qv,$trunc_len);
		    $seq = shift @trunc_res;
		    $qv =shift @trunc_res;

                    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
                    $seqR = shift @trunc_resR;
                    $qvR =shift @trunc_resR;

		}

		my $strHoldF = ">".$id."_".$q."\n".$seq."\n";
		my $strHoldR = ">".$idR."_".$qR."\n".$seqR."\n";
		my $qualStrF = ">".$idQV."_".$q."\n".$qv."\n";
		my $qualStrR = ">".$idRQV."_".$qR."\n".$qvR."\n";

		$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
		$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#		print "mp -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR **************************************\n";

		if ($QVIDF == 0 && $QVIDR == 0){
		    print oFPMPF3 $strHoldF;
                    print oFPMPR3 $strHoldR;
		    if ($outputQVFlag eq 'y'){
			print oFPMPQVF3 $qualStrF;
			print oFPMPQVR3 $qualStrR;
		    }

		    if ($qv_an eq 'y'){
		
			$mqv_cntr_pass_F3++;
			$mqv_cntr_pass_R3++;
			
			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    
			my $q_h_r3=$qvR;
			my @hold_r3=split (/\s/,$q_h_r3);
			$test=@hold_r3;
			for (0..$test-1){
			    my $a = shift @hold_r3;
			    $a = $a + 1;
			    $p_matrix_R3[$_][$a]++;
			}
		    }
		}

		elsif ($QVIDF == 1 && $QVIDR == 0){
                    print oFPOUF3 $strHoldF;
                    print oFPOR3 $strHoldR;
                    if ($outputQVFlag eq 'y'){
			print oFPOQVUF3 $qualStrF;
			print oFPOQVR3 $qualStrR;
                    }

		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_R3++;
			
			my $q_h_r3=$qvR;
			my @hold_r3=split (/\s/,$q_h_r3);
			my $test=@hold_r3;
			for (0..$test-1){
			    my $a = shift @hold_r3;
			    $a = $a + 1;
			    $p_matrix_R3[$_][$a]++;
			}
		    }
                }

		elsif ($QVIDF == 0 && $QVIDR == 1){
                    print oFPOF3 $strHoldF;
                    print oFPOUR3 $strHoldR;
		    if ($outputQVFlag eq 'y'){
			print oFPOQVF3 $qualStrF;
			print oFPOQVUR3 $qualStrR;
                    }

		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_F3++;
			
			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    }
		}
                else{
                    print oFPMPUF3 $strHoldF;
                    print oFPMPUR3 $strHoldR;

                    if ($outputQVFlag eq 'y'){
                        print oFPMPQVUF3 $qualStrF;
                        print oFPMPQVUR3 $qualStrR;
                    }
                }

		$QVIDF=0;
		$QVIDR=0;

		$lF = $seq;
		$seq="";
		$id="";
		$lR = $seqR;
		$seqR="";
		$idR="";
		next;
	    }

	    elsif ($values[0]>$valuesR[0]){
		
		if ($qv_an eq 'y'){

		    $mqv_cntr_R3++;
		    
		    my $q_ah_r3=$qvR;
		    my @hold_a_r3=split (/\s/,$q_ah_r3);
		    my $test=@hold_a_r3;
		    for (0..$test-1){
			my $a = shift @hold_a_r3;
			$a = $a + 1;
			$matrix_R3[$_][$a]++;
		    }
		}

		if ($trunc_sig eq 'y'){
                    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
                    $seqR = shift @trunc_resR;
                    $qvR =shift @trunc_resR;
                }

		my $strHold = ">".$idR."_".$qR."\n".$seqR."\n";
		my $qualStr = ">".$idRQV."_".$qR."\n".$qvR."\n";
		
		$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#		print "R -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR ************************************\n";

		if ($QVIDR == 0){
		    print oFPOR3 $strHold;
		    if ($outputQVFlag eq 'y'){
			print oFPOQVR3 $qualStr;
		    }
		
		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_R3++;
			
			my $q_h_r3=$qvR;
			my @hold_r3=split (/\s/,$q_h_r3);
			my $test=@hold_r3;
			for (0..$test-1){
			    my $a = shift @hold_r3;
			    $a = $a + 1;
			    $p_matrix_R3[$_][$a]++;
			}
		    }
		} else {
		    print oFPOUR3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVUR3 $qualStr;
		     }
		}

		$QVIDR=0;

		$lF = $id;
		$lR = $seqR;
		$seqR="";
		$idR="";
		next;
	    }
            elsif ($values[0]<$valuesR[0]){
		
		if ($qv_an eq 'y'){

		    $mqv_cntr_F3++;
		    
		    my $q_ah_F3=$qv;
		    my @hold_a_f3=split (/\s/,$q_ah_F3);
		    my $test=@hold_a_f3;
		    for (0..$test-1){
			my $a = shift @hold_a_f3;
			$a = $a + 1;
			$matrix_F3[$_][$a]++;
		    }
		}
		
		if ($trunc_sig eq 'y'){
                    my @trunc_res=truncation($seq,$qv,$trunc_len);
                    $seq = shift @trunc_res;
                    $qv =shift @trunc_res;
                }

		my $strHold = ">".$id."_".$q."\n".$seq."\n";
		my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";

                $QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
#		print "F -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR ************************************ \n";

		if ($QVIDF == 0){
		
		    print oFPOF3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVF3 $qualStr;
		     }
		    
		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_F3++;
			
			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    }

		} else {
		    print oFPOUF3 $strHold;
		    if ($outputQVFlag eq 'y'){
			print oFPOQVUF3 $qualStr;
		    }

		} 

		$QVIDF=0;

		$lF = $seq;
		$seq="";
		$id="";
		$lR = $idR;
		next;
	    }
	    elsif ($values[0]==$valuesR[0] && $values[1]>$valuesR[1]){

		if ($qv_an eq 'y'){

		    $mqv_cntr_R3++;
		    
		    my $q_ah_r3=$qvR;
		    my @hold_a_r3=split (/\s/,$q_ah_r3);
		    my $test=@hold_a_r3;
		    for (0..$test-1){
			my $a = shift @hold_a_r3;
			$a = $a + 1;
			$matrix_R3[$_][$a]++;
		    }
		}
		
		if ($trunc_sig eq 'y'){
                    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
                    $seqR = shift @trunc_resR;
                    $qvR =shift @trunc_resR;
                }


                my $strHold = ">".$idR."_".$qR."\n".$seqR."\n";
		my $qualStr = ">".$idRQV."_".$qR."\n".$qvR."\n";

		$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#		print "R -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

		if ($QVIDR == 0){
		    print oFPOR3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVR3 $qualStr;
		     }
		
		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_R3++;
			
			my $q_h_r3=$qvR;
			my @hold_r3=split (/\s/,$q_h_r3);
			my $test=@hold_r3;
			for (0..$test-1){
			    my $a = shift @hold_r3;
			    $a = $a + 1;
			    $p_matrix_R3[$_][$a]++;
			}
		    }

		} else {
		    print oFPOUR3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVUR3 $qualStr;
		     }
		}

		$QVIDR=0;

		$lF = $id;
		$lR = $seqR;
		$seqR="";
		$idR="";
		next;
            }
            elsif ($values[0]==$valuesR[0] && $values[1]<$valuesR[1]){
		if ($qv_an eq 'y'){

		    $mqv_cntr_F3++;
		    
		    my $q_ah_F3=$qv;
		    my @hold_a_f3=split (/\s/,$q_ah_F3);
		    my $test=@hold_a_f3;
		    for (0..$test-1){
			my $a = shift @hold_a_f3;
			$a = $a + 1;
			$matrix_F3[$_][$a]++;
		    }
		}
		
		if ($trunc_sig eq 'y'){
                    my @trunc_res=truncation($seq,$qv,$trunc_len);
                    $seq = shift @trunc_res;
                    $qv =shift @trunc_res;
                }
 		
		my $strHold = ">".$id."_".$q."\n".$seq."\n";
		my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";

		$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
#		print "F -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

                if ($QVIDF == 0){
                    print oFPOF3 $strHold;
		    if ($outputQVFlag eq 'y'){
			print oFPOQVF3 $qualStr;
		    }

		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_F3++;
			
			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    }
                } else {
                    print oFPOUF3 $strHold;
                    if ($outputQVFlag eq 'y'){
                        print oFPOQVUF3 $qualStr;
                    }
                }

		$QVIDF=0;
		
		$lF = $seq;
		$seq="";
		$id="";
		$lR = $idR;
		next;
	    }
	    elsif ($values[0]==$valuesR[0] && $values[1]==$valuesR[1] && $values[2]>$valuesR[2]){

		if ($qv_an eq 'y'){

		    $mqv_cntr_R3++;
		    
		    my $q_ah_r3=$qvR;
		    my @hold_a_r3=split (/\s/,$q_ah_r3);
		    my $test=@hold_a_r3;
		    for (0..$test-1){
			my $a = shift @hold_a_r3;
			$a = $a + 1;
			$matrix_R3[$_][$a]++;
		    }
		}

		if ($trunc_sig eq 'y'){
                    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
                    $seqR = shift @trunc_resR;
                    $qvR =shift @trunc_resR;
                }

		my $strHold = ">".$idR."_".$qR."\n".$seqR."\n";
		my $qualStr = ">".$idRQV."_".$qR."\n".$qvR."\n";

		$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#		print "R -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

		if ($QVIDR == 0){
		    print oFPOR3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVR3 $qualStr;
		     }
		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_R3++;
			
			my $q_h_r3=$qvR;
			my @hold_r3=split (/\s/,$q_h_r3);
			my $test=@hold_r3;
			for (0..$test-1){
			    my $a = shift @hold_r3;
			    $a = $a + 1;
			    $p_matrix_R3[$_][$a]++;
			}
		    }
		} else {
		    print oFPOUR3 $strHold;
		     if ($outputQVFlag eq 'y'){
			 print oFPOQVUR3 $qualStr;
		     }
		}

		$QVIDR=0;

		$lF = $id;
		$lR = $seqR;
		$seqR="";
		$idR="";
		next;
	    }
	    elsif ($values[0]==$valuesR[0] && $values[1]==$valuesR[1] && $values[2]<$valuesR[2]){
		if ($qv_an eq 'y'){

		    $mqv_cntr_F3++;
		    
		    my $q_ah_F3=$qv;
		    my @hold_a_f3=split (/\s/,$q_ah_F3);
		    my $test=@hold_a_f3;
		    for (0..$test-1){
			my $a = shift @hold_a_f3;
			$a = $a + 1;
			$matrix_F3[$_][$a]++;
		    }
		}
		if ($trunc_sig eq 'y'){
                    my @trunc_res=truncation($seq,$qv,$trunc_len);
                    $seq = shift @trunc_res;
                    $qv =shift @trunc_res;
                }

      		my $strHold = ">".$id."_".$q."\n".$seq."\n";
		my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";
		
		$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
#		print "F -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

                if ($QVIDF == 0){
                    print oFPOF3 $strHold;
		    if ($outputQVFlag eq 'y'){
			print oFPOQVF3 $qualStr;
		    }
		    if ($qv_an eq 'y'){
			
			$mqv_cntr_pass_F3++;
			
			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    }
                } else {
                    print oFPOUF3 $strHold;
                    if ($outputQVFlag eq 'y'){
                        print oFPOQVUF3 $qualStr;
                    }

                } 

		$QVIDF=0;

		$lF = $seq;
		$lR = $idR;
		$seq="";
		$id="";
		next;
	    }

	} 
	elsif ($lF=~/_/ && $lR!~/_/ ){
	    $lR=<iFPJ>;
	    chomp $lR;
	    $lRQV=<iFPJQV>;
	    chomp $lRQV;
	    next;
	}
	elsif ($lF!~/_/ && $lR=~/_/){    
	    $lF=<iFP>;
	    chomp $lF;
	    $lFQV=<iFPQV>;
	    chomp $lFQV;
	    next;

	}
	else {
	    $lF=<iFP>;
	    $lR=<iFPJ>;
	    $lFQV=<iFPQV>;
	    $lRQV=<iFPJQV>;
	    chomp $lF;
	    chomp $lR;
	    chomp $lRQV;
	    chomp $lFQV;
	    next;

	}
    }
 
#closes the MP files because there can't be a match anymore
    close (oFPMPF3);
    close (oFPMPQVF3);
    close (oFPMPUF3);
    close (oFPMPQVUF3);
    close (oFPMPR3);
    close (oFPMPQVR3);
    close (oFPMPUR3);
    close (oFPMPQVUR3);

# if the end of file is in F3 places the last sequence in the appropriate files
    if (eof(iFP) && $id=~/_/ && $seq=~/[ACGT][.0123]+/){
	if ($qv_an eq 'y'){

	    $mqv_cntr_F3++;

	    my $q_ah_F3=$qv;
	    my @hold_a_f3=split (/\s/,$q_ah_F3);
	    my $test=@hold_a_f3;
	    for (0..$test-1){
		my $a = shift @hold_a_f3;
		$a = $a + 1;
		$matrix_F3[$_][$a]++;
	    }
	}

	if ($trunc_sig eq 'y'){
	    my @trunc_res=truncation($seq,$qv,$trunc_len);
	    $seq = shift @trunc_res;
	    $qv =shift @trunc_res;
	}

	my $strHold = ">".$id."_".$q."\n".$seq."\n";
	my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";

	$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
#	print "last in file F -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

	if ($QVIDF == 0){
	    print oFPOF3 $strHold;
	    if ($outputQVFlag eq 'y'){
		print oFPOQVF3 $qualStr;
	    }
	    if ($qv_an eq 'y'){

		$mqv_cntr_pass_F3++;
		
		my $q_h_F3=$qv;
		my @hold_f3=split (/\s/,$q_h_F3);
		my $test=@hold_f3;
		for (0..$test-1){
		    my $a = shift @hold_f3;
		    $a = $a + 1;
		    $p_matrix_F3[$_][$a]++;
		}
	    }
	} else {
	    print oFPOUF3 $strHold;
	    if ($outputQVFlag eq 'y'){
		print oFPOQVUF3 $qualStr;
	    }

	} 

	$QVIDF=0;

	$id="";
	$seq="";
    }

#if the eof is in R3 this places it in the appropriate file
    if (eof(iFPJ) && $idR=~/_/ && $seqR=~/[ACGT][.0123]+/){

	if ($qv_an eq 'y'){

	    $mqv_cntr_R3++;
	    
	    my $q_ah_r3=$qvR;
	    my @hold_a_r3=split (/\s/,$q_ah_r3);
	    my $test=@hold_a_r3;
	    for (0..$test-1){
		my $a = shift @hold_a_r3;
		$a = $a + 1;
		$matrix_R3[$_][$a]++;
	    }
	}

	if ($trunc_sig eq 'y'){
	    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
	    $seqR = shift @trunc_resR;
	    $qvR =shift @trunc_resR;
	}

	my $strHold = ">".$idR."_".$qR."\n".$seqR."\n";
	my $qualStr = ">".$idRQV."_".$qR."\n".$qvR."\n";
	
	$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#	print "last in file R -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";
	
	if ($QVIDR == 0){
	    print oFPOR3 $strHold;
	     if ($outputQVFlag eq 'y'){
		 print oFPOQVR3 $qualStr;
	     }
	    if ($qv_an eq 'y'){

		$mqv_cntr_pass_R3++;
		
		my $q_h_r3=$qvR;
		my @hold_r3=split (/\s/,$q_h_r3);
		my $test=@hold_r3;
		for (0..$test-1){
		    my $a = shift @hold_r3;
		    $a = $a + 1;
		    $p_matrix_R3[$_][$a]++;
		}
	    }
	} else {
	    print oFPOUR3 $strHold;
	     if ($outputQVFlag eq 'y'){
		 print oFPOQVUR3 $qualStr;
	     }
	}

	$QVIDR=0;
	
	$idR="";
	$seqR="";
    }

# this puts the remainder of the R3 sequences into the appropriate orphan files since the F3 file is done
    while (eof(iFP) && !eof(iFPJ) ) {
	if ($idR !~/_/){
	    $lR=<iFPJ>;
	    chomp $lR;
	    $_=$lR;
	    tr/^>//d;
	    @valuesR =split /_/,$_;
	    $qR=pop @valuesR;
	    $_=join "_",@valuesR;
	    $idR = $_;

	    $lRQV=<iFPJQV>;
	    chomp $lRQV;
	    $_=$lRQV;
	    tr/^>//d;
	    my @valuesRQV =split /_/,$_;
	    my $qRQV=pop @valuesRQV;
	    $_=join "_",@valuesRQV;
	    $idRQV = $_;
	    
	}
	if ($seqR!~/[ACGT][.0123]+/){
	    $_ = <iFPJ>;
	    chomp;
	    $seqR=$_;
	    $_=<iFPJQV>;
	    chomp;
	    $qvR=$_;
	}
	if ($qv_an eq 'y'){

	    $mqv_cntr_R3++;
	    
	    my $q_ah_r3=$qvR;
	    my @hold_a_r3=split (/\s/,$q_ah_r3);
	    my $test=@hold_a_r3;
	    for (0..$test-1){
		my $a = shift @hold_a_r3;
		$a = $a + 1;
		$matrix_R3[$_][$a]++;
	    }
	}
	if ($trunc_sig eq 'y'){
	    my @trunc_resR=truncation($seqR,$qvR,$trunc_len);
	    $seqR = shift @trunc_resR;
	    $qvR =shift @trunc_resR;
	}

	my $strHold = ">".$idR."_".$qR."\n".$seqR."\n";
	my $qualStr = ">".$idRQV."_".$qR."\n".$qvR."\n";

	$QVIDR=qualPass($qvR,$avg_sig,$avg_thresh,$negQV);
#	print "end frag file R -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

	if ($QVIDR == 0){
	    print oFPOR3 $strHold;
	     if ($outputQVFlag eq 'y'){
		 print oFPOQVR3 $qualStr;
	     }
	    if ($qv_an eq 'y'){

		$mqv_cntr_pass_R3++;

		my $q_h_r3=$qvR;
		my @hold_r3=split (/\s/,$q_h_r3);
		my $test=@hold_r3;
		for (0..$test-1){
		    my $a = shift @hold_r3;
		    $a = $a + 1;
		    $p_matrix_R3[$_][$a]++;
		}
	    }
	} else {
	    print oFPOUR3 $strHold;
	     if ($outputQVFlag eq 'y'){
		 print oFPOQVUR3 $qualStr;
	     }
	}

	$QVIDR=0;

	$idR="";
	$seqR="";

    }  
# this puts the remainder of the F3 sequences into the appropriate orphan files since the R3 file is done
    while (!eof(iFP) && eof(iFPJ) ) {
	if ($id !~/_/){
	    $lF=<iFP>;
	    chomp $lF;
	    $_=$lF;
	    tr/^>//d;
	    @values =split /_/,$_;
	    $q=pop @values;
	    $_=join "_",@values;
	    $id = $_;

	    $lFQV=<iFPQV>;
	    chomp $lFQV;
	    $_=$lFQV;
	    tr/^>//d;
	    my @valuesQV =split /_/,$_;
	    my $qQV=pop @valuesQV;
	    $_=join "_",@valuesQV;
	    $idQV = $_;
	}
	if ($seq!~/[ACGT][.0123]+/){
	    $_ = <iFP>;
	    chomp;
	    $seq=$_;
	    $_=<iFPQV>;
	    chomp;
	    $qv=$_;
	}
	if ($qv_an eq 'y'){

	    $mqv_cntr_F3++;
	    
	    my $q_ah_F3=$qv;
	    my @hold_a_f3=split (/\s/,$q_ah_F3);
	    my $test=@hold_a_f3;
	    for (0..$test-1){
		my $a = shift @hold_a_f3;
		$a = $a + 1;
		$matrix_F3[$_][$a]++;
	    }
	}
        if ($trunc_sig eq 'y'){
            my @trunc_res=truncation($seq,$qv,$trunc_len);
            $seq = shift @trunc_res;
            $qv =shift @trunc_res;
        }

	my $strHold = ">".$id."_".$q."\n".$seq."\n";
	my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";

	$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);
#	print "end frag file F -- $id -- $q -- $seq -- $qv -- $QVIDF [] $idR -- $qR -- $seqR -- $qvR -- $QVIDR  ************************************\n";

	if ($QVIDF == 0){
	    print oFPOF3 $strHold;
	    if ($outputQVFlag eq 'y'){
		print oFPOQVF3 $qualStr;
	    }
	    if ($qv_an eq 'y'){

		$mqv_cntr_pass_F3++;
		
		my $q_h_F3=$qv;
		my @hold_f3=split (/\s/,$q_h_F3);
		my $test=@hold_f3;
		for (0..$test-1){
		    my $a = shift @hold_f3;
		    $a = $a + 1;
		    $p_matrix_F3[$_][$a]++;
		}
	    }
	} else {
	    print oFPOUF3 $strHold;
	    if ($outputQVFlag eq 'y'){
		print oFPOQVUF3 $qualStr;
	    }

	} 

	$QVIDF=0;

	$id="";
	$seq="";
}

    close (iFP);
    close (iFPJ);
    close (iFPQV);
    close (iFPJQV);
    close (oFPOF3);
    close (oFPOQVF3);
    close (oFPOUF3);
    close (oFPOQVUF3);
    close (oFPOR3);
    close (oFPOQVR3);
    close (oFPOUR3);
    close (oFPOQVUR3);

    print "finished analyzing now starting to print QV analysis\n";

    if ($qv_an eq 'y'){

	print oFPQVAN "Count of reads in QV Analysis - F3 Total: $mqv_cntr_F3 \n";
	print oFPQVAN "Count of reads in QV Analysis - R3 Total: $mqv_cntr_R3 \n";
	print oFPQVAN "Count of reads in QV Analysis - F3 Passed: $mqv_cntr_pass_F3 \n";
	print oFPQVAN "Count of reads in QV Analysis - R3 Passed: $mqv_cntr_pass_R3 \n";
	print oFPQVAN "\n\n";
	
	print oFPQVAN "Full F3 QV Analysis\n";
	for my $i ( 0 .. $#matrix_F3) {
	    print oFPQVAN "@{$matrix_F3[$i]}\n";
	}
	print oFPQVAN "\n\n";
	
	print oFPQVAN "Full R3 QV Analysis\n";
	for my $i ( 0 .. $#matrix_R3) {
	    print oFPQVAN "@{$matrix_R3[$i]}\n";
	}
	print oFPQVAN "\n\n";
	
	print oFPQVAN "Passed F3 QV Analysis\n";
	for my $i ( 0 .. $#p_matrix_F3) {
	    print oFPQVAN "@{$p_matrix_F3[$i]}\n";
	}
	print oFPQVAN "\n\n";
	
	print oFPQVAN "Passed R3 QV Analysis\n";
	for my $i ( 0 .. $#p_matrix_R3) {
	    print oFPQVAN "@{$p_matrix_R3[$i]}\n";
	}
	
	close (oFPQVAN);
    }

} else {   ############ Fragment analysis###########################

    my @holdI=split /_/,$inputF;
    my $holdI=pop @holdI;
    @holdI=split /\./,$holdI;#/
    my $qf=shift @holdI;

# opens the input files
    open (iFP,$inputF) or die "Can not open $inputF.\n";
    open (iFPQV,$inputQVF) or die "Can not open $inputQVF.\n";

# creates all the output files
    my $outputMPF3=">".$outputFile."_T_F3.csfasta";
    open (oFPMPF3,$outputMPF3) or die "Can not open $outputMPF3.\n";
    
    if ($outputQVFlag eq 'y'){
        my $outputMPQVF3=">".$outputFile."_T_F3_QV.qual";
        open (oFPMPQVF3,$outputMPQVF3) or die "Can not open $outputMPQVF3.\n";
    }

    my $outputMPUF3=">".$outputFile."_U_F3.csfasta";
    open (oFPMPUF3,$outputMPUF3) or die "Can not open $outputMPUF3.\n";
    
    if ($outputQVFlag eq 'y'){
        my $outputMPQVUF3=">".$outputFile."_U_F3_QV.qual";
        open (oFPMPQVUF3,$outputMPQVUF3) or die "Can not open $outputMPQVUF3.\n";
    }

    if ($qv_an eq 'y'){
        my $qv_an_output=">".$outputFile."_QV_analysis.txt";
        open (oFPQVAN,$qv_an_output) or die "Can not open $qv_an_output.\n";
    }

###### QV ANALYSIS ####### 
    my $mqv_cntr_F3=0;
    my $mqv_cntr_pass_F3=0;

    my @matrix_F3;
    my @p_matrix_F3;

    for my $l ( 0 .. 55 ) {
        for my $m ( 0 .. 50 ) {
            $matrix_F3[$l][$m]=0;
            $p_matrix_F3[$l][$m]=0;
        }
    }
##############################

    my $counter=0;
    my $lF = "";
    my $lFQV="";
    
    my $holdTitleA="";
    my $holdTitleB="";
    
    while (!eof(iFP)){
        if ($counter == 0 ){
            $counter =1;
            $lF=<iFP>;
            $lFQV=<iFPQV>;
            chomp $lF;
            chomp $lFQV;
        }
	
        $_=$lF;
	while (/^\#/ && !eof(iFP)) {
            $_=<iFP>;
            chomp;
            $lF=$_;
	}

        $_=$lFQV;
        while (/^\#/ && !eof(iFPQV)) {
            $_=<iFPQV>;
            chomp;
            $lFQV=$_;
        }

        last;
    }

    my $id="";
    my $seq="";
    my $qv="";
    my @values;
    my $idQV;
    $counter=1;
    my $QVIDF=0;
    my $q;

    print "\n"."Starting Anaysis - Fragment\n";

# this is the part that does the work.  Both files are open and we use the fact that SOLiD keeps the list in bead order to our advantage
# we basically start going through the list if F3 > R3 then the R3 sequence gets written into the orphan file, if F3< R3 then the F3 sequence
# get written to the orphan file.  If they match they get written to the MP files.  The threshold gets calculated inside each check.
# the threshold is a separate subroutine at the bottom of the script.

    while (!eof(iFP)) {
        if ($lF=~/_/) {
            if ($id !~/_/){
                $_=$lF;
                tr/^>//d;
                @values =split /_/,$_;
                $q=pop @values;
                $_=join "_",@values;
                $id = $_;

                $_=$lFQV;
                tr/^>//d;
                my @valuesQV =split /_/,$_;
                my $qQV=pop @valuesQV;
                $_=join "_",@valuesQV;
                $idQV = $_;
		
##              die "lF and lFQV dont match $lF - $lFQV" unless $id == $idQV;                                    
            }

            if ($seq!~/[ACGT][.0123]+/){
                $_ = <iFP>;
                chomp;
                $seq=$_;
                $_=<iFPQV>;
                chomp;
                $qv=$_;
            }

	    if($seq=~/[ACGT][.0123]+/ && $lF=~/_/){

		if ($qv_an eq 'y'){

                    $mqv_cntr_F3++;
		    
                    my $q_ah_F3=$qv;
                    my @hold_a_f3=split (/\s/,$q_ah_F3);
                    my $test=@hold_a_f3;
                    for (0..$test-1){
                        my $a = shift @hold_a_f3;
#                        print "$a ";
			$a = $a + 1;
#			print "- $a \n";
			$matrix_F3[$_][$a]++;
                    }
                }
                if ($trunc_sig eq 'y'){
                    my @trunc_res=truncation($seq,$qv,$trunc_len);
                    $seq = shift @trunc_res;
                    $qv =shift @trunc_res;
                }
		
                my $strHold = ">".$id."_".$q."\n".$seq."\n";
                my $qualStr = ">".$idQV."_".$q."\n".$qv."\n";		
		
		$QVIDF=qualPass($qv,$avg_sig,$avg_thresh,$negQV);

                if ($QVIDF == 0){
                    print oFPMPF3 $strHold;
		    if ($outputQVFlag eq 'y'){
			print oFPMPQVF3 $qualStr;
		    }
		    if ($qv_an eq 'y'){

			$mqv_cntr_pass_F3++;

			my $q_h_F3=$qv;
			my @hold_f3=split (/\s/,$q_h_F3);
			my $test=@hold_f3;
			for (0..$test-1){
			    my $a = shift @hold_f3;
			    $a = $a + 1;
			    $p_matrix_F3[$_][$a]++;
			}
		    }
		} else {
		    print oFPMPUF3 $strHold;
                    if ($outputQVFlag eq 'y'){
                        print oFPMPQVUF3 $qualStr;
                    }
                }

                $QVIDF=0;

                $lF = $seq;
                $seq="";
                $id="";
                next;
            }
	}

        else {
            $lF=<iFP>;
            $lFQV=<iFPQV>;
            chomp $lF;
            chomp $lFQV;
            next;
	}
    }

    close (oFPMPF3);
    close (oFPMPQVF3);
    close (oFPMPUF3);
    close (oFPMPQVUF3);
    
    close (iFP);
    close (iFPQV);

    if ($qv_an eq 'y'){

        print oFPQVAN "Count of reads in QV Analysis - Fragment Total: $mqv_cntr_F3 \n";
        print oFPQVAN "Count of reads in QV Analysis - Fragment Passed: $mqv_cntr_pass_F3 \n";
        print oFPQVAN "\n\n";

        print oFPQVAN "Full QV Analysis\n";
        for my $i ( 0 .. $#matrix_F3) {
            print oFPQVAN "@{$matrix_F3[$i]}\n";
        }
        print oFPQVAN "\n\n";

        print oFPQVAN "Passed Threshold QV Analysis\n";
        for my $i ( 0 .. $#p_matrix_F3) {
            print oFPQVAN "@{$p_matrix_F3[$i]}\n";
        }
        print oFPQVAN "\n\n";

        close (oFPQVAN);
    }

}

}
############################################
########SUBROUTINES#########################
############################################

sub AssignInputs{
    my @err;

#------------------------------
    if ((!defined($input_t))||($input_t eq '')||($input_t eq 'F')||($input_t eq 'Frag')||($input_t eq 'Fragment')||($input_t eq 'f')||($input_t eq 'frag')||($input_t eq 'fragment')) {
	$input_t = 'F';
	print "Type of Input - $input_t\n";
    }
    elsif (($input_t eq 'mp')||($input_t eq 'MP')||($input_t eq 'mate-pair')||($input_t eq 'Mate-Pair')) {
	$input_t='MP';
	print "Type of Input - $input_t\n";
    }
    else {
	push(@err,"Error 1: Improper Input type please use F or Frag for fragment libraries or MP for mate pair libraries.\n Default is a Fragment library.\n");
    }
#------------------------------
    if ((!defined($inputF)) || ($inputF eq '')) {
	push(@err, "Error 2: Input File not defined\n");
	}
#------------------------------
 if (((!defined($inputQVF)) || ($inputQVF eq '')) && ((!defined($inputF)) || ($inputF eq ''))) {
     push(@err, "Error 3: Input File not defined\n");
 }
    elsif (((!defined($inputQVF)) || ($inputQVF eq ''))){
        my $string = $inputF;
        $string =~ s/.csfasta$/_QV.qual/;
        $inputQVF = $string;
    }
#------------------------------
    if (((!defined($inputR)) || ($inputR eq '')) && ($input_t eq 'MP')) {
        push(@err,"Error 4: Input File not defined\n");
    }
#------------------------------
    if (((!defined($inputQVR)) || ($inputQVR eq '')) && ($input_t eq 'MP') && ((!defined($inputR)) || ($inputR eq ''))) {
        push(@err,"Error 5: Input File not defined\n");
    }
    elsif (((!defined($inputQVR)) || ($inputQVR eq '')) && ($input_t eq 'MP')){
        my $string = $inputR;
        $string =~ s/.csfasta$/_QV.qual/;
        $inputQVR = $string;
    }
#------------------------------
    if ((!defined($avg_sig)) || ($avg_sig eq '') || ($avg_sig eq 'on') || ($avg_sig eq 'y')||($avg_sig eq 'yes')){
        $avg_sig = 'y';
        print "Average Analysis Filter - $avg_sig\n";
    }
    elsif (($avg_sig eq 'off')||($avg_sig eq 'n')||($avg_sig eq 'no')) {
        $avg_sig = 'n';
        print "Average Analysis Filter - $avg_sig\n";
    }
    else {
        push(@err,"Error 6: Average Analysis Filter can either be 'on' or 'off'.  The default is on.\n");
    }
#-----------------------------------------------------------------------------------------
    if ((!defined($avg_thresh)) || ($avg_thresh eq '')){
	$avg_thresh = '20';
	print "QV for Average Analysis - $avg_thresh\n";
    }
    elsif (($avg_thresh >= 0) &&  ($avg_thresh <= 34)) {
	$avg_thresh = $avg_thresh;
	print "QV for Average Analysis - $avg_thresh\n";
    }
    else {
	push(@err,"Error 7: Improper polyclonal QV value. The values must be with in the SOLiD QV range [0-34].\n");
    }
#------------------------------
    if ((!defined($negQV)) || ($negQV eq '') || ($negQV eq 'off') || ($negQV eq 'n')||($negQV eq 'no')){
	$negQV = 'n';
	print "Removal of reads with negative QV score - $negQV\n";
    }
     elsif (($negQV eq 'on')||($negQV eq 'y')||($negQV eq 'yes')) {
         $negQV = 'y';
         print "Removal of reads with negative QV score - $negQV\n";
     }
    else {
        push(@err,"Error 8: Removal of reads with negative quality scores can either be 'on' or 'off'.  The default is off.\n");
    }
#------------------------------
    if ((!defined($trunc_sig)) || ($trunc_sig eq '') || ($trunc_sig eq 'off') || ($trunc_sig eq 'n')||($trunc_sig eq 'no')){
	$trunc_sig = 'n';
	print "Truncation is off - $trunc_sig\n";
    }
    elsif (($trunc_sig eq 'on')||($trunc_sig eq 'y')||($trunc_sig eq 'yes')) {
        $trunc_sig = 'y';
        print "Truncation is on - $trunc_sig\n";
    }
    else {
        push(@err,"Error 9: Truncation can either be 'on' or 'off'.  The default is off.\n");
    }
#------------------------------
    if (((!defined($trunc_len)) || ($trunc_len eq '')) && ($trunc_sig eq 'n') ){
    }
    elsif (((!defined($trunc_len)) || ($trunc_len eq '')) && ($trunc_sig eq 'y') ){
	push(@err,"Error 10: Improper Truncation length.  Truncation is turned on please pass a desired length for the reads.  Any bases after that length will be eliminated. \n");
 
    }
    elsif (($trunc_len >= 0) && ($trunc_sig eq 'y')) {
        $trunc_len = $trunc_len;
        print "Truncation length - $trunc_len\n";
    }
    else {
        push(@err,"Error 10: Improper truncation length.  Truncation must be turned on and the length must be greater than or equal to 0.\n");
    }
#------------------------------
    if ((!defined($qv_an)) || ($qv_an eq '') || ($qv_an eq 'off') || ($qv_an eq 'n')||($qv_an eq 'no')){
        $qv_an = 'n';
        print "Quality Value analysis is off - $qv_an\n";
    }
    elsif (($qv_an eq 'on')||($qv_an eq 'y')||($qv_an eq 'yes')) {
        $qv_an = 'y';
        print "Quality Value analysis is on - $qv_an\n";
    }
    else {
        push(@err,"Error 11: Quality Value analysis can either be 'on' or 'off'.  The default is off.\n");
    }
#------------------------------
    if ((!defined($outputFile)) || ($outputFile eq '')) {
        push(@err,"Error 12: Output File not defined\n");
    }
#------------------------------
    if ((!defined($outputQVFlag)) || ($outputQVFlag eq '') || ($outputQVFlag eq 'on')){
        $outputQVFlag = 'y';
        print "Output QV files - $outputQVFlag\n";
    }
    elsif (($outputQVFlag eq 'off')){
        $outputQVFlag = 'n';
        print "Output QV files - $outputQVFlag\n";
    }
    else {
        push(@err,"Error 13: Improper QV Flag value.  Please use on or off to turn on or off QV file printing.\n");
    }
#------------------------------
    return @err;
}

sub qualPass{
    #checks if the quality score passes the poluclonal and the error tests

    my $qScores = shift;
    my $avg_s    = shift;
    my $avg_qv   = shift;
    my $neg_s   = shift;

    my @holdI =split / /,$qScores;
    my $holdVal=0;
    my $sum = 0;
    my $cntr=0;
    my $neg_sc_i=0;

    foreach my $t (@holdI){
	$cntr++;
	$sum = $sum + $t;
	if (($t < 0) && ($neg_s eq 'y')){
	    $neg_sc_i=1;
	}
    }

    my $avg = $sum/$cntr;
#    print "QualPass - $sum - $cntr - $avg ---------------------------------------\n";
    
    if ((($avg_s eq 'y') && ($avg < $avg_qv)) || ($neg_sc_i == 1)){
	$holdVal = 1;
    }

    return $holdVal;
}

sub truncation{
    
    my $seq1 = shift;
    my $qv1 =shift;
    my $t_len = shift;
    my @res;

    my $seq_len_f = length ($seq1);

    if ($t_len > $seq_len_f){
	$t_len = $seq_len_f;
	$res[0]=$seq1;
	$res[1]=$qv1;
	return @res;
    }
    
    my @hold_qv_f =split / /,$qv1;
    my @hold_seq_f = split //,$seq1;
    
    my @t_seq_f = @hold_seq_f[0..$t_len];
    my $t_qv_len = $t_len -1;
    my @t_qv_f = @hold_qv_f[0..$t_qv_len];

    $res[0]=join('',@t_seq_f);
    $res[1]=join(' ',@t_qv_f);

    return @res;
}

sub usage {

    print "\nusage: $0 \n[-f csfasta file for the forward strand \n -g quality file for the forward strand \n\t [-i input type \n\t -r csfasta file for the reverse strand \n\t -s quality file for the reverse strand \n\t -p Mean anlysis on/off \n\t -q Mean analysis min QV score \n\t -n Removal of all reads containing a missing colorcall on/off \n\t -t truncation on/off \n\t -u truncation length \n\t -a quality analysis on/off] \n -o output file name \n\t [-v Output matching QV files on/off]\n]\n\n";

    print "This is a filtering script that filters out the SOLiD csfasta and QV files \n"."based on the mean quality values of the read.  For Fragments it only separates based \n"."on passing or non passing scores. For Mate-pair data this script prints a \n"."mimic of the original Solid output but separating mate pairs from orphans \n"."with the additional criteris of passing the two Quality tests.  But in order \n"."for the mates to remain together they must have the same passing score otherwise \n"."they will fall separately into orphan files. Here is usage information\n\n";

    print "-f|f3 - csfasta file for the forward strand - required \n\t"."- this the name of the file coming from the SOLiD primary analysis \n"; 
    print "-g|f3QV - quality file for the forward strand - optional \n\t"."- this is the name of the file coming from the SOLiD primary analysis.  \n\t  If this is left blank, the name of the corresponding QV file must match \n\t  and be in the same location as the csfasta file except the ending will replace \n\t  the .csfasta with _QV.qual.  \n\n";

    print "-i|input_type - the input files type - optional (default fragment analysis) \n\t"."- identifies mate pair from fragment analysis. if the -i option is \n\t  left blank or excluded a Fragment analysis will be preformed and if there \n\t  is input on the -r and -s options it will be ignored.  For mate-pair analysis, \n\t  MP option must be selected. \n\t  [F,Frag,Fragment,f,frag,fragment,mp,MP,mate-pair,Mate-Pair] \n\n";

    print "-r|r3 - csfasta file for the reverse strand - optional \n\t"."- this is the name of the mate csfasta file to the f3 option. \n";
    print "-s|r3QV - quality file for the reverse strand - optional \n\t"."- this is the name of the mate's QV file. If this is left blank, the name of the \n\t  corresponding QV file must match and be in the same location as the csfasta \n\t  file except the ending will replace \n\t  the .csfasta with _QV.qual.  \n\n";

    print "-p|mean_analysis - Average anlysis on/off - optional (default on) \n\t"."- Average analysis filter looks at the average of the all the color calls. \n\t  It asks that the average of all the calls must exceed the qv score of \n\t  interest (p|mean_qv). \n\t  [on,yes,y,off,no,n]\n";  
    print "-q|mean_qv - Average analysis min QV score - optional (default 20) \n\t"."- This is the minimum score for the average  analysis.  \n\t  This is currently limited to a number between 0-34.\n\n";

    print "-n|neg_qv - Removal of all reads containing a negative quality score on/off \n\t"."- optional (default off) \n\t"."- Removes all reads with missing colorcalls (i.e. negative quality scores). \n\t  [on,yes,y,off,no,n]\n\n";

    print "-t|trunc - truncation on/off - optional (default off) \n\t"."- Turns truncation of the 3' end of the read on or off. \n\t  [on,yes,y,off,no,n]\n";
    print "-u|tr_len - desired length of read after truncation - optional \n\t"."- This is the length of the sequence desired, any color calls \n\t  after this length are removed.  This option must be filled in if \n\t  truncation is turned on and be an integer greater than 0. \n\n";

    print "-a|qv_analysis - turn on/off analysis of quality values - optional (default off) \n\t"."- Analysis of the quality values for all of the inputted reads and \n\t  the passing reads.  Analysis returns a file with a matrix of a count \n\t  of scores by position. \n\t  [on,yes,y,off,no,n]\n\n";

    print "-o|output - output file name - required \n\t"."- this is the begining of the name for the output information. \n\t  The endings are filled in as required.\n"; 
    print "-v|output_qv - Output matching QV files on/off - optional (default on) \n\t"."- this will print the matching QV files to the csfasta files. \n\t  [on,off]\n\n";

    print "example: $0 -i mp -f a_F3.csfasta -g a_F3_QV.qual \n -r b_R3.csfasta -s b_R3_QV.qual -p 3 -q 22 -e 10 -d 9 -v off -o test_test \n";

    exit;

}

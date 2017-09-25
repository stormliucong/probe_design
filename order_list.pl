#!usr/bin/perl -w
use strict;
use warnings;


my %count_gene;
open(IN,"$ARGV[0]") or die "$!";
while (my $line = <IN>) {
    if($line=~/>(.+?)\|/){
        my $gene = $1;        
        
        $line = <IN>;
        $line=~/cpg\#=(\d+)\|Tm=([\d,\.]+)\|length=(\d+)bp\|CG\%=([\d,\.]+)\%/;
        my $left_cpg = $1;
        my $left_tm = $2;
        my $left_length = $3;
        my $left_cg = $4;

        $line = <IN>;
        chomp $line;
        my $left_seq = $line;

        $line = <IN>;
        $line=~/cpg\#=(\d+)\|Tm=([\d,\.]+)\|length=(\d+)bp\|CG\%=([\d,\.]+)\%/;
        my $right_cpg = $1;
        my $right_tm = $2;
        my $right_length = $3;
        my $right_cg = $4;

        $line = <IN>;
        chomp $line;
        my $right_seq = $line;

        my $bs_seq = &bsConvert($left_seq.$right_seq);
        my $bs_left_seq = substr($bs_seq,0,length($left_seq));
        my $bs_right_seq = substr($bs_seq,length($left_seq),length($right_seq));

        if($left_cg > 30 or $right_cg > 30){
            if(not exists $count_gene{$gene}){
                $count_gene{$gene} = 1;
            }else{
                $count_gene{$gene} +=1;
            }
            my $gene_count = $count_gene{$gene};
            my $prime_name_mc = $gene."_".$gene_count."_mc";
            my $prime_name_F = $gene."_".$gene_count."-F-I";
            my $prime_name_R = $gene."_".$gene_count."-R-I";

            my $mc_seq = $left_seq.$right_seq;
            my $bp_count = length($mc_seq);
            print "1","\t",$prime_name_mc,"\t",$mc_seq,"\t",$bp_count,"\n";

            my $F_seq_1 = &reverseComplementDNA($bs_left_seq);
            my $F_seq_2 = "NNNNNNNN";
            my $F_seq_3 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC";
            my $F_seq = $F_seq_1.$F_seq_2.$F_seq_3;
            print "2","\t",$prime_name_F,"\t",$F_seq,"\t",length($F_seq),"\t",join("\t",$left_cpg,$left_tm,$left_length,$left_cg),"\n";

            my $R_seq_3 = &reverseComplementDNA($bs_right_seq);
            my $R_seq_2 = "NNNNNNNN";
            my $R_seq_1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
            my $R_seq = $R_seq_1.$R_seq_2.$R_seq_3;
            print "3","\t",$prime_name_R,"\t",$R_seq,"\t",length($R_seq),"\t",join("\t",$right_cpg,$right_tm,$right_length,$right_cg),"\n";
        }else{
            if(not exists $count_gene{$gene}){
                $count_gene{$gene} = 1;
            }else{
                if($count_gene{$gene} > 2){
                    next;
                }
                $count_gene{$gene} +=1;
            }
            
            my $gene_count = $count_gene{$gene};
            my $prime_name_mc = $gene."_".$gene_count."_mc";
            my $prime_name_F = $gene."_".$gene_count."-F-I";
            my $prime_name_R = $gene."_".$gene_count."-R-I";

            my $mc_seq = $left_seq.$right_seq;
            my $bp_count = length($mc_seq);
            print "1","\t",$prime_name_mc,"\t",$mc_seq,"\t",$bp_count,"\n";

            my $F_seq_1 = &reverseComplementDNA($bs_left_seq);
            my $F_seq_2 = "NNNNNNNN";
            my $F_seq_3 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTC";
            my $F_seq = $F_seq_1.$F_seq_2.$F_seq_3;
            print "2","\t",$prime_name_F,"\t",$F_seq,"\t",length($F_seq),"\t",join("\t",$left_cpg,$left_tm,$left_length,$left_cg),"\n";

            my $R_seq_3 = &reverseComplementDNA($bs_right_seq);
            my $R_seq_2 = "NNNNNNNN";
            my $R_seq_1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT";
            my $R_seq = $R_seq_1.$R_seq_2.$R_seq_3;
            print "3","\t",$prime_name_R,"\t",$R_seq,"\t",length($R_seq),"\t",join("\t",$right_cpg,$right_tm,$right_length,$right_cg),"\n";

        }

        
        
    }
    # body...
}

sub reverseComplementDNA (){
    my $dna = shift;
    my $revcom =reverse($dna);  
    $revcom=~tr/ACGTacgt/TGCAtgca/;  
    return $revcom;  
}

sub bsConvert (){
    my $seq = shift;
    $seq=~s/CG/MM/g;
    $seq=~s/C/T/g;
    $seq=~s/MM/CG/g;
    return($seq);
}

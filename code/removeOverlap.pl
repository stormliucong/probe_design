#!usr/bin/perl -w
use strict;

open(IN,$ARGV[0]) or die "$!";
my $previous_strand = "NULL";
my $previous_end = 0;
my $previous_gene = "NULL";
my $i = 0;
my $j = 0;
my %hash_output;
while (my $line1 = <IN>) {
    if($line1 =~/^>(.+?)\|(.+?)\|(.+?)\|(\d+)\|(\d+)/){
        my $gene = $1;
        my $chr = $2;
        my $strand = $3;
        my $c_pos = $4;
        my $g_pos = $5;
        my $line2 = <IN>;
        $line2=~/length=(\d+)bp/;
        my $l_length = $1;
        my $line3 = <IN>;
        my $line4 = <IN>;
        $line4=~/length=(\d+)bp/;
        my $line5 = <IN>;
        my $r_length = $1;
        my $start = $c_pos - $l_length;
        my $end = $g_pos + $r_length;
        if($previous_strand eq $strand and $previous_gene eq $gene){
            if($previous_end < $start){
                $i++;
                $j = 1;
                push @{$hash_output{$gene}{$strand}{$i}{$j}},$line1,$line2,$line3,$line4,$line5;
                #print $line1;
                #print $line2;
                #print $line3;
                #print $line4;
                #print $line5;
            }else{
                $j++;
                push @{$hash_output{$gene}{$strand}{$i}{$j}},$line1,$line2,$line3,$line4,$line5;

            }
        }else{
            $i = 1;
            $j = 1;
            push @{$hash_output{$gene}{$strand}{$i}{$j}},$line1,$line2,$line3,$line4,$line5;
            #print $line1;
            #print $line2;
            #print $line3;
            #print $line4;
            #print $line5;
        }
        $previous_end = $end;
        $previous_strand = $strand;
        $previous_gene = $gene;

    }
}
close IN or die "$!";

foreach my $gene (keys %hash_output){
    foreach my $strand (keys %{$hash_output{$gene}}){
        foreach my $i (keys %{$hash_output{$gene}{$strand}}){
            my @j = keys %{$hash_output{$gene}{$strand}{$i}};
            my $lower_limit = 0;
            my $upper_limit = $#j;
            my $random_number = int(rand($upper_limit-$lower_limit)) + $lower_limit;
            my $j = $j[$random_number];
            print ${$hash_output{$gene}{$strand}{$i}{$j}}[0];
            print ${$hash_output{$gene}{$strand}{$i}{$j}}[1];
            print ${$hash_output{$gene}{$strand}{$i}{$j}}[2];
            print ${$hash_output{$gene}{$strand}{$i}{$j}}[3];
            print ${$hash_output{$gene}{$strand}{$i}{$j}}[4];


        }
    }
}

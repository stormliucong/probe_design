#!usr/bin/perl -w
use strict;

# input: mctest_falist.txt is a list contains all sam file names (with path).
# output: a list of SEQXXXX.sam.flagStat.
# For each probe (each sequence in fa, summarize the number of three situations.
# (1) both ends are mapped by proper paired reads;
# (2) One ends are mapped by one end of reads, and the other end of reads are mapped to another probe;
# (3) One ends are mapped by one end of reads, and the other end of reads are not mapped to any other probe;

my @file_list;
open(FA,"mctest_falist.txt") or die "$!";
while(my $file = <FA>){
    chomp $file;
    push @file_list, $file;
}
close FA or die "$!";

foreach my $file (@file_list){
    print $file,"\n";
    $file=~/\/(SEQ\d+).sam/;
    my $sample_id = $1;
    my @ele;
    my (%stat,%pair_proper_map,%pair_inproper_map,%single_map);
    my $output = $sample_id."\.flagStat";
    open(IN,$file) or die "$!";
    print "reading ".$file,"\n";
    while(my $line = <IN>){
        if($line=~/^\@/){
            next;
        }else{
            chomp $line;
            @ele = split("\t",$line);
            my $flag = $ele[1];
            my $probe = $ele[2];
            #print $flag,"\t",$probe,"\n";
            if(not exists $stat{$probe}{$flag}){
                $stat{$probe}{$flag} = 0;
            }
            $stat{$probe}{$flag} = $stat{$probe}{$flag} + 1;
        }
    }
    close IN or die "$!";

    print "printing to ".$output,"\n";
    foreach my $probe (keys %stat){
        foreach my $flag (keys %{$stat{$probe}}){
            # read paired and read mapped in proper pair.
            if($flag eq "67" or $flag eq "83" or $flag eq "99" or $flag eq "115" or $flag eq "131" or $flag eq "147" or $flag eq "163" or $flag eq "179"){
                if(not exists $pair_proper_map{$probe}){
                    $pair_proper_map{$probe} = 0;
                }
                $pair_proper_map{$probe} = $stat{$probe}{$flag} + 1;
            }
            # read paired and read mapped in inproper pair.
            if($flag eq "65" or $flag eq "81" or $flag eq "97" or $flag eq "113" or $flag eq "129" or $flag eq "145" or $flag eq "161" or $flag eq "177"){
                if(not exists $pair_inproper_map{$probe}){
                    $pair_inproper_map{$probe} = 0;
                }
                $pair_inproper_map{$probe} = $stat{$probe}{$flag} + 1;
            }

            # read paired and read single mapped.
            if($flag eq "73" or $flag eq "89" or $flag eq "105" or $flag eq "121" or $flag eq "137" or $flag eq "153" or $flag eq "169" or $flag eq "185"){
                if(not exists $single_map{$probe}){
                    $single_map{$probe} = 0;
                }
                $single_map{$probe} = $stat{$probe}{$flag} + 1;
            }
        }
    }

    open(OUT,">".$output) or die "$!";
    foreach my $probe (keys %pair_proper_map){
        print OUT ($probe,"\t",$pair_proper_map{$probe},"\t",$pair_inproper_map{$probe},"\t",$single_map{$probe},"\n");
    }
    close OUT or die "$!";
}


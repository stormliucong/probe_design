#usr/bin/perl -w
# help to correct the bug caused by betaMatrixGenerator.pl.
use strict;
use List::MoreUtils qw(uniq);

my %bigMatHash;
my %fileHash;

my $line_count = 485577; # there are 485577 cpg in each files.
my $file_count = int(485577/5000) + 1;

opendir(DIR,$ARGV[0]) or die "$!";
my @d_name = readdir(DIR);

for my $file_seq (1..$file_count){

    print "processing file seq:".$file_seq."\n";
    foreach my $d_name (@d_name){
        if(-d $d_name & $d_name=~/(\w+)/){
            print "processing sample :".$d_name."\n";
        
            # whether it is a fold.
            opendir(SDIR,"./".$d_name) or die "$!";
            my @f_name = readdir(SDIR);
            closedir SDIR or die "$!";
            foreach my $f_name (@f_name){
                
                if($f_name=~/jhu-usc\.edu_(.+?)\.HumanMethylation450/){

                    # record file ids.
                    my $id = $d_name;
                    if(not exists $fileHash{$id}){
                            $fileHash{$id} = 1;
                    }
                }else{
                    print "./".$d_name."/".$f_name." DOES NOT fit the file name format\n";
                    next;
                }   
                # only read splitted files.
                if($f_name=~/_$file_seq\.tmp/){
                    my $id = $d_name;
                    
                    open(IN,"./".$d_name."/".$f_name) or die "$!";
                    my $line = <IN>; # that is a bug cause the number of cpg in each file is 4999 instead of 5000.
                    
                    chomp $line;
                    my @ele = split("\t",$line);
                    my $cpg = $ele[0];
                    my $beta = $ele[1];
                        #my $chr = $ele[2];
                        
                    $bigMatHash{$cpg}{$id} = $beta;
                    close IN or die "$!";

                }else{
                    next;
                }
            }
        }else{
            print $d_name." is not a sample.\n";
        }
    }
}
# print big matrix in splitted files.
# column seq is the same as row seq in meta data.
my @file_name;
open(IN,"betaMat_sample_group.txt") or die "$!";
while(my $line = <IN>){
    my @ele = split("\t",$line);
    push @file_name, $ele[0];
}
close IN or die "$!";

open(OUT,">betaMat"."_99"."\.txt") or die "$!";
foreach my $cpg (keys %bigMatHash){
print OUT ($cpg);
    foreach my $id (@file_name){
        if(not exists $bigMatHash{$cpg}{$id}){
            #print "NA-".$cpg."-".$id."\n";
            print OUT ("\t"."NA");
        }else{
            print OUT ("\t",$bigMatHash{$cpg}{$id});
        }
    }
    print OUT "\n";

}
close OUT or die "$!";









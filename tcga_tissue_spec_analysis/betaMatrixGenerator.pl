#usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

# define HASH;
my %tissueHash;
my %fileHash;
my @cpg_annotation;

# open the dir containing sample folders.
opendir(DIR,$ARGV[0]) or die "$!";


my @d_name = readdir(DIR);
closedir DIR or die "$!";

# first split file to save RAM.

my $file_1 = 0; # store cpg annotation only once.
open(CPG,">betaMat_cpg_annotation.txt") or die "$!";
foreach my $d_name (@d_name){
    if(-d $d_name & $d_name=~/(\w+)/){
        print "processing sample :".$d_name."\n";
        
        # whether it is a fold.
        opendir(SDIR,"./".$d_name) or die "$!";
        my @f_name = readdir(SDIR);
        closedir SDIR or die "$!";
        foreach my $f_name (@f_name){
            if(-f "./".$d_name."/".$f_name & $f_name=~/(\w+)/){
                if($f_name=~/jhu-usc\.edu_(.+?)\.HumanMethylation450/){

                    # record file ids.
                    my $id = $d_name;
                    if(not exists $fileHash{$id}){
                            $fileHash{$id} = 1;
                    }

                    # record tissue information.
                    my $tissue = $1;
                    if(not exists $tissueHash{$id}){
                        $tissueHash{$id} = $tissue;
                    }else{
                        if($tissueHash{$id} ne $tissue){
                            die "the information contained in splitted files not equal.\n";
                        }
                    }

                    # split files.
                    open(IN,"./".$d_name."/".$f_name) or die "$!";
                    my $line = <IN>;
                    my $line_count = 0;
                    my $file_count = 1;
                    $file_1++;

                    open(OUT,">"."./".$d_name."/".$f_name."_".$file_count.".tmp") or die "$!";    
                    while(my $line = <IN>){
                        print OUT $line;
                        $line_count++;


                        if($file_1 == 1){
                            # restore the cpg annotation only once.
                            my @ele = split("\t",$line);
                            $ele[1] = $file_count;
                            print CPG (join("\t",@ele));
                        }
                        

                        # start write to a new file for 5000 cpg.
                        if(($line_count % 5000) == 0){
                            $file_count++;
                            close OUT or die "$!";
                            open(OUT,">"."./".$d_name."/".$f_name."_".$file_count.".tmp") or die "$!";
                        }

                        
                        
                    }
                    close OUT or die "$!";
                }else{
                    print "./".$d_name."/".$f_name." DOES NOT fit the file name format\n";
                }
            }else{
                # it is a log folder;
                next;
            }
        }
    }else{
        print $d_name." is not a sample\n";
    }
}
close CPG or die "$!";


# then read and write for each splitted file. 

my $line_count = 485577; # there are 485577 cpg in each files.
my $file_count = int(485577/5000) + 1;

for my $file_seq (1..$file_count){

    # define hash here to reuse the memory.
    my %bigMatHash;
    #my %cpgLocalHash;

    print "processing file seq:".$file_seq."\n";
    foreach my $d_name (@d_name){
        if(-d $d_name & $d_name=~/(\w+)/){
            print "processing sample :".$d_name."\n";
        
            # whether it is a fold.
            opendir(SDIR,"./".$d_name) or die "$!";
            my @f_name = readdir(SDIR);
            closedir SDIR or die "$!";
            foreach my $f_name (@f_name){
                
                    
                # only read splitted files.
                if($f_name=~/_$file_seq\.tmp/){
                    my $id = $d_name;
                    
                    open(IN,"./".$d_name."/".$f_name) or die "$!";
                    my $line = <IN>; # that is a bug cause the number of cpg in each file is 4999 instead of 5000.
                    while($line = <IN>){
                        chomp $line;
                        my @ele = split("\t",$line);
                        my $cpg = $ele[0];
                        my $beta = $ele[1];
                        #my $chr = $ele[2];
                        
                        $bigMatHash{$cpg}{$id} = $beta;
                    
                    }
                    close IN or die "$!";

                }else{
                    next;
                }
            }
        }else{
            print $d_name." is not a sample.\n";
        }
    }
    # print big matrix in splitted files.
    # column seq is the same as row seq in meta data.
    open(OUT,">betaMat"."_".$file_seq."\.txt") or die "$!";
    foreach my $cpg (keys %bigMatHash){
        print OUT ($cpg);
        foreach my $id (keys %fileHash){
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
}


# print meta data.
open(OUT,">betaMat_sample_group.txt") or die "$!";
foreach my $id (keys %fileHash){
    print OUT ($id,"\t",$tissueHash{$id},"\n");
}
close OUT or die "$!";






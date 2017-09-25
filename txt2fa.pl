#!usr/bin/perl -w
use strict;

open(IN,$ARGV[0]) or die "$!";
# my $line = <IN>;
# print $line;die;
my ($ref_id,$seq_2,$seq_1);
my @ele;
while(my $line = <IN>){
    chomp $line;
    @ele = split(" ",$line);
    if($ele[0] eq "1"){
        $ref_id = $ele[1];
        print ">".$ref_id,"\n";
        $line = <IN>;
        chomp $line;
        @ele = split(" ",$line);
        if($ele[0] eq "2"){
            $ele[2]=~/^(\w+)NNNNNNNN/;
            $seq_2 = $1;
        }else{
            die "$!";
        }
        $line = <IN>;
        chomp $line;
        @ele = split(" ",$line);
        #print join("\t",@ele),"\n";
        #print $ele[0],"\n";
        if($ele[0] eq "3"){
            $ele[2]=~/NNNNNNNN(\w+)$/;
            $seq_1 = $1;

        }else{
            die "$!";
        }
        print $seq_1.$seq_2,"\n";
    }else{
        die "$!";
    }
}    

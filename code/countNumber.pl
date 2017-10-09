#!usr/bin/perl -w
use strict;
use warnings;

my %count;
my %length;
open(IN,$ARGV[0]) or die "$!";
while(my $line = <IN>){
	if($line =~/>(.+?)\|/){
		my $gene = $1;
		$count{$gene}++;
		$line =~/chr(.+?):(\d+)-(\d+)/;
		$length{$gene} = $3-$2;
	}
}
close IN or die "$!";

foreach my $gene (keys %count){
	print $count{$gene},"\t",$length{$gene},"\t",$gene,"\n";
}


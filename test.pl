#!usr/bin/perl -w
use strict;
use warnings;

my $seq = "GGGTCCGACACCGACGTCGGTCGTTCC";
substr($seq,-1) = "T";
print $seq,"\n";


# my @c_site;
# while($seq=~/CG/g){

#     push @c_site,pos($seq)-2;
# }
# print join(";",@c_site),"\n";
# foreach my $c_site (@c_site){
#     print substr($seq,$c_site,2),"\n";
# }
# my $rev = reverse($seq);
# print $rev,"\n";
# # sub reverseComplementDNA (){
# #     my $dna = shift;
# #     $dna=~tr/ACGTacgt/TGCAtgca/;  
# #     $dna=reverse($dna);  
# #     return $dna;  
# # }

# sub reverseComplementDNA (){
#     my $dna = shift;
#     print $dna,"\n";
#     my $revcom=reverse($dna);  
#     print $revcom,"\n";
#     $revcom=~tr/ACGTacgt/TGCAtgca/; 
#     print $revcom,"\n";
#     return $revcom;  
# }
# &reverseComplementDNA($seq);
#print &reverseComplementDNA($seq),"\n";
# # my @c_site;
# # my $seq= "CGGCACGTCGTGACAC";
# # while($seq=~/CG/g){
# #             push @c_site,pos($seq)-2;
# #         }
# # print join(";",@c_site),"\n";


# my %hash_fa = &readFA($ARGV[0]);
# my $seq = &fetchSeq("chr7",5568179,5568219,\%hash_fa);
# $seq = uc($seq); 
# print $seq,"\n";
# sub readFA (){
#     my $fa_file = shift;
#     my %hash_fa;
#     my $head;
#     open(IN,$fa_file) or die "$!";
#     while (my $line = <IN>) {
#         chomp $line;
#         if($line=~/^>(.+?)$/){
#             $head = $1;
#             if(not exists $hash_fa{$head}){
#                 $hash_fa{$head} = "";
#             }else{
#                 warn "Not proper file!\n";
#                 die "Same Head\n";
#             }
#         }else{
#             $hash_fa{$head} .= $line;
#         }
#     }
#     close IN or die "$!";
#     return %hash_fa;
# }
# sub fetchSeq (){
#     my $chrm = shift;
#     my $start = shift;
#     my $end = shift;
#     my $hash_fa_ref = shift;
#     my %hash_fa = %{$hash_fa_ref};
#     my $chr_seq = $hash_fa{$chrm};
#     print length($chr_seq),"\n";
#     print $start,"\n";
#     my $length = $end - $start + 1;
#     print substr($chr_seq,0,100),"\n";
#     print substr($chr_seq,-100),"\n";
#     my $seq = substr($chr_seq,$start,$length);
#     return($seq);
# }
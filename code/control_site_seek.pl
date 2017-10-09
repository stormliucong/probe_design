#!usr/bin/perl -w
use strict;
use warnings;


my $gtf_file = $ARGV[0];
my $fa_file = $ARGV[1];
my $cancer_marker_file = $ARGV[2];
my $output_file = $ARGV[3];


&main($gtf_file,$fa_file,$cancer_marker_file,$output_file);

# main function.
sub main (){
    my $gtf_file = shift;
    my $fa_file = shift;
    my $cancer_marker_file = shift;
    my $output_file = shift;
    my %hash_gtf;
    # read gtf and store the results in an hash.
    %hash_gtf = &readGTF($gtf_file);
    my %hash_fa;
    # read hg19 and store the results in a hash.
    %hash_fa = &readFA($fa_file);
    # read gene list and store the results in a hash.
    my @cancer_marker = &readList($cancer_marker_file);

    open(OUT,">".$output_file) or die "$!";
    # search candidate probe for each gene.
    foreach my $gene (@cancer_marker) {
        if(not exists $hash_gtf{$gene}){
            warn "$gene not exists in gtf files.\n";
            next;
        }
        # get the gene coordinate information from stored gtf hash.
        my $chrm = ${$hash_gtf{$gene}}[0];
        my $start = ${$hash_gtf{$gene}}[1];
        my $end = ${$hash_gtf{$gene}}[2];
        my $cor = $chrm.":".$start."-".$end;

        # skip very short or long genes.
        if($end - $start < 100){
            warn "$gene length is shorter than 100 bp.\n";
            next;
        }
        if($end - $start > 10000000){
            warn "$gene length is larger than 10M bp.\n";
            next;
        }

        # get the sequence of the gene + up and down stream 5kb seq.
        my $seq = &fetchSeq($chrm,$start-1500,$end+1500,\%hash_fa);
        $seq = uc($seq); # convert all seq letters to upper case (capitalized).

        # start from "+" strand and search candidate probe from flanking sequence of each CpG site.
        my $strand = "+";
        my $msg = join("|",$gene,$cor,$strand)."\n"; 
        warn "$msg";

        # select 50 candidate probes for each control gene.
        my $n = 0;
        my %info;
        while($n < 20){
            # select seq randomly.
            my $lower_limit = 0 + 20;
            my $upper_limit = length($seq) - 20;
            my $c_site = int(rand($upper_limit-$lower_limit)) + $lower_limit;
            my $g_site = $c_site + 1;

            # avoid randomly pick CpG site.
            if(substr($seq,$c_site,2) eq "CG"){
                next;
            }
            my ($hash_left_vote_ref,$hash_right_vote_ref,$left_vote,$right_vote) = &printCandidateProbe($c_site,$seq);

            if($hash_left_vote_ref){
                my %hash_left_vote = %{$hash_left_vote_ref};
                my %hash_right_vote = %{$hash_right_vote_ref};
                foreach my $left_key (keys %hash_left_vote){
                    if($left_key ne $left_vote){
                        next;
                    }
                    foreach my $right_key (keys %hash_right_vote){
                        if($right_key ne $right_vote){
                            next;
                        }
                        $info{$c_site} = ">".join("|",$gene,$cor,$strand,$c_site,$g_site)."\n"."L@".join("|",${$hash_left_vote{$left_key}}[1],${$hash_left_vote{$left_key}}[2],${$hash_left_vote{$left_key}}[3],${$hash_left_vote{$left_key}}[4])."\n".
${$hash_left_vote{$left_key}}[0]."\n"."R@".join("|",${$hash_right_vote{$right_key}}[1],${$hash_right_vote{$right_key}}[2],${$hash_right_vote{$right_key}}[3],${$hash_right_vote{$right_key}}[4])."\n".${$hash_right_vote{$right_key}}[0]."\n";
                        $n++;
                    }
                }
            }
        }
        # sort by c_site position.
        foreach my $k (sort {$a<=>$b} keys %info){
            print OUT $info{$k};
        }
    }
    close OUT or die "$!";
}

# get reverse complement DNA seq.
sub reverseComplementDNA (){
    my $dna = shift;
    my $revcom=reverse($dna);  
    $revcom=~tr/ACGTacgt/TGCAtgca/;  
    return $revcom;  
}

# collect candidate probe information.
sub printCandidateProbe (){
    my $c_site = shift;
    my $g_site = $c_site + 1;

    my $seq = shift;

    my $left_vote = 0;
    my $right_vote = 0;   

    my %hash_left_vote;
    my %hash_right_vote; 

    my $bs_seq = &bsConvert($seq);

    my $left_length = int(rand(20)) + 20;
    my $right_length = int(rand(20)) + 20;
    my @vote;

    my $left_start = $c_site - $left_length;
    my $right_end = $g_site + $right_length;

    my $left_bs_seq = substr($bs_seq,$left_start,$left_length);
    my $left_seq = substr($seq,$left_start,$left_length);
    @vote = &criteriaVote($left_bs_seq);
    if($vote[0]){
        $left_vote++;
        ${$hash_left_vote{$left_vote}}[0] = $left_seq."|".$left_bs_seq;
        @{$hash_left_vote{$left_vote}}[1..4] = @vote;            
    }

    my $right_bs_seq = substr($bs_seq,$g_site,$right_length);
    my $right_seq = substr($seq,$g_site,$right_length);
    @vote = &criteriaVote($right_bs_seq);
    if($vote[0]){
        $right_vote++;
        ${$hash_right_vote{$right_vote}}[0] = $right_seq."|".$right_bs_seq;
        @{$hash_right_vote{$right_vote}}[1..4] = @vote; 
    }
    #print $right_vote,"\t",$left_vote,"\n";
    if(exists $hash_left_vote{$left_vote} and exists $hash_right_vote{$right_vote}){
        return (\%hash_left_vote,\%hash_right_vote,$left_vote,$right_vote);
    }else{
        return(0);
    }
}

# get DNA seq based on chrm, start and end.
sub fetchSeq (){
    my $chrm = shift;
    my $start = shift;
    my $end = shift;
    my $hash_fa_ref = shift;
    my %hash_fa = %{$hash_fa_ref};
    my $chr_seq = $hash_fa{$chrm};
    my $length = $end - $start + 1;

    my $seq = substr($chr_seq,$start,$length);
    return($seq);
}

# read gene list.
sub readList (){
    my $cancer_marker_file = shift;
    open(IN,$cancer_marker_file) or die "$!";
    my @cancer_marker;
    while(my $line = <IN>){
        chomp $line;
        $line =~s/\r//g;
        push @cancer_marker,$line;
    }
    return @cancer_marker;
}

# read human genome.
sub readFA (){
    my $fa_file = shift;
    my %hash_fa;
    my $head;
    open(IN,$fa_file) or die "$!";
    while (my $line = <IN>) {
        chomp $line;
        if($line=~/^>(.+?)$/){
            $head = $1;
            if(not exists $hash_fa{$head}){
                $hash_fa{$head} = "";
            }else{
                warn "Not proper file!\n";
                die "Same Head\n";
            }
        }else{
            $hash_fa{$head} .= $line; # bug revised in 05-19-2016.
        }
    }
    close IN or die "$!";
    return %hash_fa;
}

# read gtf human genome annotation.
sub readGTF (){
    my $gtf_file = shift;
    my %hash_gtf;
    open(IN,$gtf_file) or die "$!";
    while (my $line = <IN>) {
        chomp $line;
        if($line=~/gene_name \"(.+?)\"/){
            my $gene = $1;
            my @ele = split("\t",$line);
            my $chrm = $ele[0];
            my $start = $ele[3];
            my $end = $ele[4];
            if(not exists $hash_gtf{$gene}){
                push @{$hash_gtf{$gene}}, $chrm,$start,$end;
            }else{
                if(${$hash_gtf{$gene}}[0] ne $chrm){
                    next;
                }
                if(${$hash_gtf{$gene}}[1] > $start){
                    ${$hash_gtf{$gene}}[1] = $start;
                }
                if(${$hash_gtf{$gene}}[2] < $end){
                    ${$hash_gtf{$gene}}[2] = $end;
                }
            }
        }
    }
    close IN or die "$!";
    return %hash_gtf;
}

# count c+g number.
sub cgCount (){
    my $seq = shift;
    my $c = 0;
    my $g = 0;
    while($seq=~/C/g){
        $c++;
    }
    while($seq=~/G/g){
        $g++;
    }
    return($c+$g);
}

# count cpg number.
sub cpgCount (){
    my $seq = shift;
    my $cpg = 0;
    while($seq=~/CG/g){
        $cpg++;
    }
    return($cpg);
}

# justify whether the candidate seq satisfy the criteria.
sub criteriaVote (){
    my $seq = shift;
    my $cpg = &cpgCount($seq);
    my $wak_tm = &WAKTm($seq);
    my $seq_length = length($seq);
    #print $bs_seq,"\n";
    #print &cgCount($seq),"\n";
    my $cg_content = &cgCount($seq)/$seq_length;


    if($cpg > 0){
        #warn "cpg=",$cpg,"\n";
        return(0);
    }
    if($wak_tm < 60 or $wak_tm > 70){
        #warn "tmp=",$wak_tm,"\n";
        return(0);
    }
    if($seq_length < 20 or $seq_length > 40){
        #warn "length=",$seq_length,"\n";
        return(0);
    }
    if($cg_content < 0.2 or $cg_content > 0.8){
        #warn "cg_content=",$cg_content,"\n";
        return(0);
    }
    $cpg = "cpg#=".$cpg;
    $wak_tm = "Tm=".sprintf("%.1f", $wak_tm);
    $cg_content = "CG%=".sprintf("%.1f", $cg_content*100)."%";
    $seq_length = "length=".$seq_length."bp";

    return($cpg,$wak_tm,$seq_length,$cg_content);
}

# bulsifite conversion of the seq.
sub bsConvert (){
    my $seq = shift;
    $seq=~s/CG/MM/g;
    $seq=~s/C/T/g;
    $seq=~s/MM/CG/g;
    return($seq);
}

# calculate the Tm for a given seq.
sub WAKTm_old (){
    my $seq = shift;
    my $cg_count = &cgCount($seq);
    my $na = 50;
    if (length($seq) > 0) {
        if (length($seq) < 14) {
            return (2*(length($seq)-$cg_count)+4*$cg_count+21.6+(7.21*log($na/1000)));
        }
        else {
            return (100.5+(0.41*$cg_count/length($seq)*100)-(820/length($seq))+(7.21*log($na/1000)));

        }
    }else {
            die("too short");
    }
}

# this calculation is based on operon.
sub WAKTm (){
    my $seq = shift;
    my $cg_count = &cgCount($seq);
    my $na = 0.33;
    if (length($seq) > 14){
        return (81.5+16.6*log($na)+41*$cg_count/length($seq)-500/length($seq));
    }else{
        warn("too short");
        return (0);
    }
}
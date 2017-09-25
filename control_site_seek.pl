#!usr/bin/perl -w
use strict;
use warnings;


my $gtf_file = $ARGV[0];
my $fa_file = $ARGV[1];
my $cancer_marker_file = $ARGV[2];


&main($gtf_file,$fa_file,$cancer_marker_file);

# main function.
sub main (){
    my $gtf_file = shift;
    my $fa_file = shift;
    my $cancer_marker_file = shift;
    my %hash_gtf;
    # read gtf and store the results in an hash.
    %hash_gtf = &readGTF($gtf_file);
    my %hash_fa;
    # read hg19 and store the results in a hash.
    %hash_fa = &readFA($fa_file);
    # read gene list and store the results in a hash.
    my @cancer_marker = &readList($cancer_marker_file);

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
        if($end - $start < 10){
            warn "$gene length is shorter than 20 bp.\n";
            next;
        }
        if($end - $start > 10000000){
            warn "$gene length is larger than 10M bp.\n";
            next;
        }

        # get the sequence of the gene + up and down stream 5kb seq.
        my $seq = &fetchSeq($chrm,$start-500,$end+500,\%hash_fa);
        $seq = uc($seq); # convert all seq letters to upper case (capitalized).

        # start from "+" strand and search candidate probe from flanking sequence of each CpG site.
        my $strand = "+";
        my @c_site = 40..(length($seq)-40); # is not a C_sie. but the end site of the first probe.
        my $msg = join("|",$gene,$cor,$strand)."\n"; 
        warn "$msg";

        # search candidate probe from flanking sequence of each CpG site.
        for(my $i = 1;$i<$#c_site;$i++){
            my $c_site = $c_site[$i];

            # if(($c_site[$i] - $c_site[$i-1]) > 40 or ($c_site[$i+1] - $c_site[$i]) > 40){
            #     next;
            # }

            my $g_site = $c_site + 1; # not actual g_site. start of second probe.
            # judge whether it is a qualified candidate probe.
            my $cpg_check = substr($seq,$c_site,2);
            if($cpg_check eq "CG"){
                # Skip CG site.
                next;
            }
            my ($hash_left_vote_ref,$hash_right_vote_ref,$left_vote,$right_vote) = &printCandidateProbe($c_site,$seq);
            if($hash_left_vote_ref){
                #print "Yes\n";
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
                        print ">";
                        print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";
                        print "L@";
                        print join("|",${$hash_left_vote{$left_key}}[1],${$hash_left_vote{$left_key}}[2],${$hash_left_vote{$left_key}}[3],${$hash_left_vote{$left_key}}[4]),"\n";                
                        print ${$hash_left_vote{$left_key}}[0],"\n";
                        print "R@";
                        print join("|",${$hash_right_vote{$right_key}}[1],${$hash_right_vote{$right_key}}[2],${$hash_right_vote{$right_key}}[3],${$hash_right_vote{$right_key}}[4]),"\n";                
                        print ${$hash_right_vote{$right_key}}[0],"\n";
                    }
                }
            }else{
                next;
            }
        }

        # repeat the procedure above for "-" strand
        $seq = reverseComplementDNA($seq);
        $strand = "-";
        @c_site = 40..(length($seq)-40);
        $msg = join("|",$gene,$cor,$strand)."\n"; 
        warn "$msg";
        for(my $i = 1;$i<$#c_site;$i++){
            my $c_site = $c_site[$i];

            my $g_site = $c_site + 1;
            my $cpg_check = substr($seq,$c_site,2);
            if($cpg_check eq "CG"){
                next;
            }
            my ($hash_left_vote_ref,$hash_right_vote_ref,$left_vote,$right_vote) = &printCandidateProbe($c_site,$seq);
            #print $hash_left_vote_ref,"\n";
            if($hash_left_vote_ref){
                #print "Yes\n";
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
                        print ">";
                        print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";
                        print "L@";
                        print join("|",${$hash_left_vote{$left_key}}[1],${$hash_left_vote{$left_key}}[2],${$hash_left_vote{$left_key}}[3],${$hash_left_vote{$left_key}}[4]),"\n";                
                        print ${$hash_left_vote{$left_key}}[0],"\n";
                        print "R@";
                        print join("|",${$hash_right_vote{$right_key}}[1],${$hash_right_vote{$right_key}}[2],${$hash_right_vote{$right_key}}[3],${$hash_right_vote{$right_key}}[4]),"\n";                
                        print ${$hash_right_vote{$right_key}}[0],"\n";
                    }
                }
            }else{
                next;
            }
        }


    }
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
    my $seq = shift;
    my $left_vote = 0;
    my $right_vote = 0;   

    my %hash_left_vote;
    my %hash_right_vote; 

    my $bs_seq = &bsConvert($seq);

    my $g_site = $c_site + 1;

    my $left_start;
    my $right_end;
    my $left_seq;
    my $left_bs_seq;
    my $right_seq;
    my $right_bs_seq;

    for(my $left_length = 40;$left_length > 19;$left_length--){
        if($c_site-$left_length + 1 < 0){
            $left_start = 0;
            $left_bs_seq = substr($bs_seq,$left_start,$c_site + 1);
            $left_seq = substr($seq,$left_start,$c_site + 1);

            my @vote = &criteriaVote($left_bs_seq);
            if($vote[0]){
                $left_vote++;
                #push @{$hash_left_vote{$left_vote}},$left_seq,@vote;
                ${$hash_left_vote{$left_vote}}[0] = $left_seq;
                @{$hash_left_vote{$left_vote}}[1..4] = @vote;
                last;

                #print "1\n";
                # print ">";
                # print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";                
                # print $left_seq,"\n";
                # print "+";
                # print join("|",@vote),"\n";
            }
            last;
        }
        $left_start = $c_site-$left_length + 1;
        $left_bs_seq = substr($bs_seq,$left_start,$left_length);
        $left_seq = substr($seq,$left_start,$left_length);
        my @vote = &criteriaVote($left_bs_seq);
        if($vote[0]){
            $left_vote++;
            #push @{$hash_left_vote{$left_vote}},$left_seq,@vote;
            ${$hash_left_vote{$left_vote}}[0] = $left_seq;
            @{$hash_left_vote{$left_vote}}[1..4] = @vote;
            last;
            # print ">";
            # print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";                
            # print $left_seq,"\n";
            # print "+";
            # print join("|",@vote),"\n";
        }
    }

    for(my $right_length = 40;$right_length > 19;$right_length--){
        if($g_site+$right_length - 1 > (length($seq)-1)){

            $right_end = length($seq)-1;
            $right_bs_seq = substr($bs_seq,$g_site,$right_end-$g_site+1);
            $right_seq = substr($seq,$g_site,$right_end-$g_site+1);

            my @vote = &criteriaVote($right_bs_seq);
            if($vote[0]){
                $right_vote++;
                #push @{$hash_right_vote{$right_vote}},$right_seq,@vote;
                ${$hash_right_vote{$right_vote}}[0] = $right_seq;
                @{$hash_right_vote{$right_vote}}[1..4] = @vote;
                last;

                # print ">";
                # print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";                
                # print $right_seq,"\n";
                # print "-";
                # print join("|",@vote),"\n";
            }
            last;
        }


        $right_bs_seq = substr($bs_seq,$g_site,$right_length);
        $right_seq = substr($seq,$g_site,$right_length);

        my @vote = &criteriaVote($right_bs_seq);
        if($vote[0]){
            $right_vote++;
            #push @{$hash_right_vote{$right_vote}},$right_seq,@vote;
            ${$hash_right_vote{$right_vote}}[0] = $right_seq;
            @{$hash_right_vote{$right_vote}}[1..4] = @vote;
            last;
            # print ">";
            # print join("|",$gene,$cor,$strand,$c_site,$g_site),"\n";                
            # print $right_seq,"\n";
            # print "-";
            # print join("|",@vote),"\n";
        }
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
    #my $bs_seq = &bsConvert($seq);
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
sub WAKTm (){
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

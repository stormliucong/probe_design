
# coding: utf-8

# In[1]:

# class perform sequence alignment using dynamic programming.

class DynamicProgrammer(object):
    int_to_char = {0:'A', 1:'C', 2:'G', 3:'T', 4:'N'}
    def __init__(self, seq1, seq2,indel,scoring):
        self.seq1 = seq1
        self.seq2 = seq2
        self.indel = indel
        self.scoring = scoring
        self.D = None

    def find_gobal_alignment(self):
        self.D = np.zeros((self.seq1.size+1, self.seq2.size+1), dtype=np.int16)
        self._compute_array()
        # self._print_sequences()
        return self._get_final_score()

    def _compute_array(self):
        for i in xrange(self.seq1.size+1):
            self.D[i,0] = i*self.indel
        for j in xrange(self.seq2.size+1):
            self.D[0,j] = j*self.indel
        for i in xrange(1, self.seq1.size+1):
            for j in xrange(1, self.seq2.size+1):
                self.D[i,j] = max(  self.D[i-1, j-1] + self._get_score(i, j),
                                    self.D[i-1, j] + self.indel,
                                    self.D[i, j-1] + self.indel)
    def _get_score(self, i, j):
        ''' The indexing is quite tricky because the matrix as one more row & column.
        That causes a shift between the matrix index and the sequence indices.
        Therefore, to obtain the correct nucleotide in the sequence, we must
        substract 1 to the matrix index. '''
        return self.scoring[self.seq1[i-1], self.seq2[j-1]]
    
    def _get_final_score(self):
        i = self.seq1.size
        j = self.seq2.size
        return self.D[i,j]
    
    def _get_aligned_pair(self, i, j):
        n1 = int_to_char[self.seq1[i-1]] if i>0 else '_'
        n2 = int_to_char[self.seq2[j-1]] if j>0 else '_'
        return (n1, n2)

    def _traceback(self):
        alignment= []
        i = self.seq1.size
        j = self.seq2.size
        while i >0 and j>0:
            if self.D[i-1, j-1] + self._get_score(i, j) == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, j))
                i -= 1
                j -= 1
            elif self.D[i-1, j] + self.indel == self.D[i,j]:
                alignment.append(self._get_aligned_pair(i, 0))
                i -= 1
            else:
                alignment.append(self._get_aligned_pair(0, j))
                j -= 1
        while i > 0:
            alignment.append(self._get_aligned_pair(i, 0))
            i -= 1
        while j > 0:
            alignment.append(self._get_aligned_pair(0, j))
            j -= 1
        alignment.reverse()
        return alignment
    
    def _print_sequences(self):
        pairs = self._traceback()
        top_seq = []
        bottom_seq = []
        for (b, t) in pairs:
            bottom_seq.append(b)
            top_seq.append(t)
        for n in top_seq:
            print n,
        print ' '
        for n in bottom_seq:
            print n,
        print ("")


# Above class is a class to perform dynamic programming for sequence alignment. Notice that the class have to used within barcode_compare function defined below. **barcode_compare** function is to compare two sequences using dynamic programming ,difflib or other methods. The scoring matrix and penalty score could be altered. By default, 
# * DP (dynamic programming) is to final optimal global alignment of two sequences.
#     * indel = -2
#     * match = 1
#     * mimatch = -1
#     * N = 0
# * difflib is a python module to compare two strings. N won't be considered specifically. Only uses default score in **difflib** module.
# * other is a method to compare two sequences base by base from 5' to 3'. Notice two sequences should be in the same length. Otherwise error will be raised.
#     * match = 1
#     * mimatch = -1
#     * N = 0
# 

# In[2]:

# function perform barcode comaprison and return the similarity score.

def barcode_compare (seq1, seq2, method = "DP"):
    '''
    Input: 
    seq1 is the a barcode sequence.
    seq2 is another barcode sequence.
    There are three methods to compare the similarity of two barcodes.
     - DP: dynamic programming
     - difflib: a python package to compare two string similarity
     - other: compare two barcode base by base from 5' to 3'
    Output:
    final_score describe the differences.
     - DP: similarity matrix based dynamic programming score.
     - difflib. raito represent the similarity.
     - other: score calculated base by base.
    
    '''
    
    char_to_int = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    
    if method == "DP":
        # make sure the sequence comapred is DNA sequences. (N is permitted).
        valid_dna = "ATCGN"
        is_DNA1 = all(i in valid_dna for i in seq1)
        is_DNA2 = all(i in valid_dna for i in seq1)
        assert(is_DNA1 & is_DNA2)

        
        # define penalty score.
        indel = -2 
        scoring = array([[1,-1,-1,-1,0],
                         [-1,1,-1,-1,0],
                         [-1,-1,1,-1,0],
                         [-1,-1,-1,1,0],
                         [0,0,0,0,0]])

        s1 = array([char_to_int[i] for i in seq1], dtype=np.int16)
        s2 = array([char_to_int[i] for i in seq2], dtype=np.int16)
        aligner = DynamicProgrammer(s1, s2,indel,scoring)
        final_score = aligner.find_gobal_alignment()
    elif method == "difflib":
        s = df.SequenceMatcher(lambda x: x == "N",seq1,seq2)
        final_score = s.ratio()
    else:
        assert(len(seq1) == len(seq2))
        final_score = 0
        match_score = 1
        N_score = 0
        mimatch_score = -1
        for i in range(len(seq1)):
            if seq1[i] != "N" and seq2[i] != "N":
                if seq1[i] == seq2[i]:
                    final_score += match_score
                else:
                    final_score += mimatch_score
            else:
                final_score += N_score
    
    return(final_score)


# Next we will show some examples to observe how to works for different barcode comparison. From the result we could see that (for barcode with length eight):
# * perfect match : 8/1/8
# * perfect match w N:  N: 7/0.875/7
# * one mismatch: 6/0.875/6
# * one indel: 5/0.933/NA 
# * one mismatch with N: 5/0.625/5
# * two mismatches: 4/0.625/4
# 
# ** Therefore, a single-end cut-off could be 6/0.875/6.**

# In[3]:

# # an example show barcode_compare()
# import difflib as df
# import numpy as np
# from numpy import array
# import re
    
# bcode1 = "ACTGGACG"
# bcode2 = "ACTGGACG"
# r1 = barcode_compare(bcode1,bcode2) # dynamic comparison.
# r2 = barcode_compare(bcode1,bcode2,method="difflib") # diff.ratio comparison.
# r3= barcode_compare(bcode1,bcode2,method="other") # custermized one-to-one comparison.

# print("Perfectly match r1:" + str(r1))
# print("Perfectly match r2:" + str(r2))
# print("Perfectly match r3:" + str(r3))


# bcode1 = "ACTGGACG"
# bcode2 = "ACTGNACG"
# r1 = barcode_compare(bcode1,bcode2)
# r2 = barcode_compare(bcode1,bcode2,method="difflib")
# r3= barcode_compare(bcode1,bcode2,method="other")

# print("Perfectly match but with 1N r1:" + str(r1))
# print("Perfectly match but with 1N r2:" + str(r2))
# print("Perfectly match but with 1N r3:" + str(r3))


# bcode1 = "ACTGGACG"
# bcode2 = "ACTGCACG"
# r1 = barcode_compare(bcode1,bcode2)
# r2 = barcode_compare(bcode1,bcode2,method="difflib")
# r3= barcode_compare(bcode1,bcode2,method="other")


# print("One mismatch r1:" + str(r1))
# print("One mismatch r2:" + str(r2))
# print("One mismatch r3:" + str(r3))



# bcode1 = "ACTGGACG"
# bcode2 = "ACTGACG"
# r1 = barcode_compare(bcode1,bcode2)
# r2 = barcode_compare(bcode1,bcode2,method="difflib")

# print("One indel r1:" + str(r1))
# print("One indel r2:" + str(r2))

# bcode1 = "ACTGGACG"
# bcode2 = "ACTGCANG"
# r1 = barcode_compare(bcode1,bcode2)
# r2 = barcode_compare(bcode1,bcode2,method="difflib")
# r3= barcode_compare(bcode1,bcode2,method="other")


# print("One mismatch with 1N r1:" + str(r1))
# print("One mismatch with 1N r2:" + str(r2))
# print("One mismatch with 1N r3:" + str(r3))


# bcode1 = "ACTGGACG"
# bcode2 = "ACTGCATG"
# r1 = barcode_compare(bcode1,bcode2)
# r2 = barcode_compare(bcode1,bcode2,method="difflib")
# r3= barcode_compare(bcode1,bcode2,method="other")
# print("Two mismatch r1:" + str(r1))
# print("Two mismatch r2:" + str(r2))
# print("Two mismatch r3:" + str(r3))


# Speed will be very different given different comparison methods. Let's test the computation time. From the results we could tell the running time relationship is 
#  ** DP = 10 $\times$ difflib $\times$ others**

# In[44]:

# this is a test to see the speed of the barcode_compare.
# import time
# bcode1 = "ACTGGACG"
# bcode2 = "ACTGNACG"
# methods = ["DP","difflib","other"]
# compare_times = [10,100,1000,10000]
# for h in range(len(compare_times)):
#     compare_time = compare_times[h]
#     for i in range(len(methods)):
#         method = methods[i]
#         tic = time.time()
#         for j in range(compare_time):
#             r = barcode_compare(bcode1,bcode2,method=method)
#         tac = time.time()
#         delta = tac - tic
#         print ("computation time for {} + {} is {}s".format(compare_time,method,delta))
# print (barcode_count)


# Our algorithm will cluster the barcode together if one barcode is similar to another by certain cut-off. The greedy algorithm is designed to cluster the barcode with lower count into barcode with higher count. There are some technical details should be noticed:
# * for each probe barcode is ordered accroding to the count both asc and descly.
# * from the lowest count barcode, barcode will compare to from the highest count barcode to decide whether to merge.
# * speed will vary in different methods. The speed of greed_cluster is $O(N^2)$
# * to reduce the computational burden, we introduce $desc\_threshold$ and $asc\_threshold$. (See comment for the meaning)
# * $score1\_threshold$ and $score2\_threshold$ should be modified according to different barcode_compare_methods

# In[66]:

# function to perform greedy cluster.
def greedy_cluster (barcode_count, barcode_compare_method = "other", score1_threshold = 6, score2_threshold = 6, desc_threshold = 99999, asc_threshold = 99999):
    '''
    Input: 
    barcode_count is a dictionary. key is the barcode pair, and
    value is the count of this pair of barcode on the probe.
    
    e.g. barcode_count = {('ACTGGACG', 'ACTGCGTT'): 100, ('ACTGGACG', 'ACTGGACGC'): 5,('ACNGGACG','TCATGACG'):2}
    
    barcode_compare_method is the method to calculate pair-wise barcode similarity.

    score1_threshold is 5'end comparison lower bound.
    score2_threshold is 3'end comparison lower bound.
    
    desc_threshold is the number of candidate barcodes to be merged into. - i.e. added.
    asc_threshold is the candidate of candidate barcodes to be merged from - i.e. deleted.
    
    Output: 
    updated barcode_count.
    '''
    
    from collections import OrderedDict

    barcode_count_asc = OrderedDict(sorted(barcode_count.items(), key=lambda kv: kv[1]))
    barcode_count_desc = OrderedDict(sorted(barcode_count.items(), key=lambda kv: -kv[1]))

    # define the candidate query set need merged.
    q_pair_set = barcode_count_asc.keys()
    
    for i in range(min(asc_threshold, len(q_pair_set))):
                
        q5 = q_pair_set[i][0]
        q3 = q_pair_set[i][1]
        count = barcode_count_asc[q_pair_set[i]]
        
        # do not compare with itself
        del barcode_count_desc[q_pair_set[i]]
        
        # the reference dict is not empty.
        if barcode_count_desc:
            
            # define the candidate reference set could be merged into.
            r_pair_set = barcode_count_desc.keys()

            for j in range(min(desc_threshold, len(r_pair_set))):
                r5 = r_pair_set[j][0]
                r3 = r_pair_set[j][1]
                
                score1 = barcode_compare(q5,r5,method=barcode_compare_method)
                score2 = barcode_compare(q3,r3,method=barcode_compare_method)        
                
                if score1 >= score1_threshold and score2 >= score2_threshold:
                    
                    barcode_count[r_pair_set[j]] += count
                    del barcode_count[q_pair_set[i]]
                    
                    break
    
    return (barcode_count)
    


# In[67]:

# this is a test to see the speed of the greedy_cluster.
# import time
# Ns = [10] # barcode_number
# for h in range(len(Ns)):
#     N = Ns[h]
#     tic = time.time()
#     for i in range(N*N/2): # O(N^2) 
#         barcode_count_orig = {('GCCTCTAC', 'AGAGGGCA'): 2, ('TAATGGCC', 'AAAGCTCG'): 2, ('GCCTCTAC', 'AGAGGTCA'): 2, ('GATAAGCT', 'CAACCGAC'): 2, ('CTGCGTCC', 'GGACCCGA'): 4, ('GTACTGGA', 'GACACCAA'): 2, ('ACCCGAGG', 'TACTACAC'): 2}
#         barcode_count_new = greedy_cluster(barcode_count_orig)
#     tac = time.time()
#     delta = tac - tic
#     print ("computation time for {} is {}s".format(N,delta))
# print (barcode_count_new)


# From the test above we could see the computation could be huge - ** $O(N^2)$ **
# for N = 1000, it takes ~1min to compute.
# Since in each dataset, we are expected to have 100,000,000 reads for one sample, and thus expect to have ~100,000 reads for each probe. The computing time could be larger than 1 year. ( The true time will be much less than that due to reads are not evenly distributed across barcodes).
# 
# ** Threfore, in practce, we have to optimize it **.
# * add a cut-off desc_threshold: only query top ranked barcode_pair in barcode_count_desc dictionary. (implmented already)
# * add a cut-off asc_threshold: only query top ranked barcode_pair in barcode_count_asc dictionary. (implmented already)
# * optimize the compare algorithm. (implmented already)
#     * DP = 10  ××  difflib  ××  others

# Next we implement read fastq and read sam file. 
# 
# Below is an example for a ** SAM ** line. SAM is ordered as fastq reads order.
# 
# ** TPNB500121:117:H7CLYBGX2:1:11101:17864:1043     73      MLH1_3_mc       1       42      30M     =       1       0
#        ACGCCAAAATATCGTTCGTAACAAAAATTA  AEEAEEEEEEEEEEEEEEEEEEEEAEEEEA  AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:30 YT:Z:UP**
# 
# The raw fastq seq for this reads is 
# 
# ** TACTANTAACGCCAAAATATCGTTCGTAACAAAAATTA **
# 
# For this reads the barcode is ** TACTANTA **

# In[8]:

# function to read and parse fastq
def read_fastq (fq1,fq2,barcode_length = 8):
    '''
    Input: 
    a pair of fastq file fq1 and fq2
    length of barcode
    Output:
    a dictionary key is reads id and value is [barcode1,barcode2]
    '''
    import gzip
    from Bio import SeqIO
    barcode_dict = {}
    for fq in [fq1,fq2]:
        with gzip.open(fq) as f:
            records = SeqIO.parse(f, "fastq")
            for record in records:
                seq_id = record.id
                seq_barcode = str(record.seq)[:barcode_length]
                if seq_id in barcode_dict.keys():
                    barcode_dict[seq_id].append(seq_barcode)
                else:
                    barcode_dict[seq_id] = []
                    barcode_dict[seq_id].append(seq_barcode)
                    
    return barcode_dict


# In[9]:

# an example to show the function above.
# fq1 = "/Users/cong/Desktop/barcode/SEQ4416_1_40000.fastq.gz"
# fq2 = "/Users/cong/Desktop/barcode/SEQ4416_2_40000.fastq.gz"
# barcode_dict = read_fastq(fq1,fq2,8)
# first2pairs = {k: barcode_dict[k] for k in barcode_dict.keys()[:2]}
# print(first2pairs)


# ** read_fastq ** function is to read and parse fastq data. The output is a dictionary where key is the reads id and value is a list of 5' and 3' barcode. Something have to be pay attention.
# * The RAM request is huge to store the whole fastq in memory.
# * Biopython should be installed.
# * What to do with NNNNNNNN?

# In[103]:

# function to read and parse SAM
def read_bam (bam_file,barcode_dict, output = 'everything'):
    '''
    Input:
    bam_file. only bam file. SAM file does not work
    barcode_dict obtained from read_fastq()
    output_barcode: Should it output barcode or reads count.
    barcode_details: Should it ouput a summary of barcode count for each probe or details.
    Output:
    bc_dict: a dictionary containning the count of barcode for each probe after greed_clustered.
    '''
    # change it from sam to bam.
#     if filename.split(".")[-1] == "sam":
#         pysam.sort("-o", bam_file+"sam", bam_file)
    
    import pysam
  
    # creat an index for bam file for the purposes of random access 

    pysam.index(bam_file) 
    
    f = pysam.AlignmentFile(bam_file,'rb')
    
    barcode_detail_dict = {}
    barcode_detail_dict["pm"] = {}
    barcode_detail_dict["mm"] = {}
    barcode_detail_dict["mu"] = {}
    
    uniq_barcode_dict = {}
    uniq_barcode_dict["pm"] = {}
    uniq_barcode_dict["mm"] = {}
    uniq_barcode_dict["mu"] = {}
    
    reads_dict = {}
    reads_dict["pm"] = {}
    reads_dict["mm"] = {}
    reads_dict["mu"] = {}
    
    for probe in (f.references):
        barcode_count_pm = {} # both ends mapped to the same ref.
        barcode_count_mm = {} # mate maps to different ref.
        barcode_count_mu = {} # mate unmapped.
        reads = f.fetch(probe)
        sam = []
        for read in reads:
            barcode = barcode_dict[read.query_name]
            barcode = tuple(barcode) # list is not hashaable.
# the reads status could be the following.
# see http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.aligned_pairs for details.
#             read.mate_is_reverse
#             read.mate_is_unmapped
#             read.is_unmapped
#             read.is_duplicate
#             read.is_paired
#             read.is_proper_pair
#             read.is_qcfail
#             read.is_read1
#             read.is_read2
#             read.is_reverse
#             read.is_secondary
#            sam.append(str(read))

            if read.is_proper_pair:
                if barcode in barcode_count_pm.keys():
                    barcode_count_pm[barcode] += 1
                else:
                    barcode_count_pm[barcode] = 1
                
            else:
                if not (read.is_unmapped):
                    if read.mate_is_unmapped:
                        if barcode in barcode_count_mu.keys():
                            barcode_count_mu[barcode] += 1
                        else:
                            barcode_count_mu[barcode] = 1
                    else:
                        if barcode in barcode_count_mm.keys():
                            barcode_count_mm[barcode] += 1
                        else:
                            barcode_count_mm[barcode] = 1
                else:
                    next # don't count unmapped reads since it already counted by its mate.
        # print len(barcode_count.keys())
        barcode_detail_dict['pm'][probe] = greedy_cluster(barcode_count=barcode_count_pm,asc_threshold=50,desc_threshold=50)
        barcode_detail_dict['mm'][probe] = greedy_cluster(barcode_count=barcode_count_mm,asc_threshold=50,desc_threshold=50)
        barcode_detail_dict['mu'][probe] = greedy_cluster(barcode_count=barcode_count_mu,asc_threshold=50,desc_threshold=50)
        
        uniq_barcode_dict['pm'][probe] = len(barcode_detail_dict['pm'][probe].keys())
        uniq_barcode_dict['mm'][probe] = len(barcode_detail_dict['mm'][probe].keys())
        uniq_barcode_dict['mu'][probe] = len(barcode_detail_dict['mu'][probe].keys())
        
        reads_dict["pm"][probe] = sum(barcode_detail_dict['pm'][probe].values())/2 # a pair will only count once. 
        reads_dict["mm"][probe] = sum(barcode_detail_dict['mm'][probe].values())
        reads_dict["mu"][probe] = sum(barcode_detail_dict['mu'][probe].values())

        
#         reads_dict[probe] = sam
    if output == 'barcode count':
        return barcode_detail_dict
    elif output == "read count":
        return reads_dict
    elif output == "barcode detail":
        return uniq_barcode_dict
    else:
        return barcode_detail_dict,reads_dict,uniq_barcode_dict





# In[105]:

# a program to test the read_bam.
# import time
# import pysam

# tic = time.time()
# sam_file = "/Users/cong/Desktop/barcode/SEQ4416_10112.sam"
# bam_file = "/Users/cong/Desktop/barcode/SEQ4416_10112.bam"
# # pysam.sort("-o", bam_file, sam_file)
# barcode_detail_dict, reads_dict, uniq_barcode_dict = read_bam(bam_file,barcode_dict)

# # print reads_dict['MAL_1_mc']
# tac = time.time()
# print (tac - tic)


# $read\_bam$ is a function to read bam and count the barcode count for each probe. it could 
# * return a dictionary containning count for each barcodes and each probe
# * OR return a dictionary containing total barcode count for each probe.
# * OR return a dictionary containing total reads count for each probe.
# 
# Three types of reads count are returned.
# * pm: read is mapped in this probe in a proper pair (a pair will only count once)
# * mm: one end of read is mapped in this probe, another in another or discordordantly
# * mu: one end of read is mapped in this probe, another unmapped
# 
# ** pysam module is required to run the code.**
# 

# Finally, we implement main function to execute the barcode count program.

# In[106]:

def main(bam_file,fq1,fq2):
    
    import time
    from collections import OrderedDict
    import difflib as df
    import numpy as np
    from numpy import array
    import re
    import csv
    import pandas as pd
    
    barcode_dict = read_fastq(fq1,fq2,8)
    barcode_detail_dict, reads_dict, uniq_barcode_dict = read_bam(bam_file,barcode_dict)
    
#     count_df = pd.DataFrame(count)
#     count_df.to_csv(output)
#     with open(output, 'wb') as csv_file:
#         csv_writer = csv.writer(csv_file)
#         for key, value in count.items():
#             csv_writer.writerow([key, value])
        
    return barcode_detail_dict, reads_dict, uniq_barcode_dict


# In[91]:

# a test program to see whether it returns the correct answer.

# import pysam
# f = pysam.AlignmentFile(bam_file,'rb')
# reads = f.fetch('CDH13_2_mc')
# barcode_dict = read_fastq(fq1,fq2,8)
# for read in reads:
#     print ('{}\t{}\t{}\t{}\t{}\t{}').format(read.query_name,barcode_dict[read.query_name],read.reference_name,read.is_proper_pair,read.is_unmapped,read.mate_is_unmapped)


# # In[107]:

import pandas as pd
import sys
# sam_file = "/Users/cong/Desktop/barcode/SEQ4416_10112.sam"
# bam_file = "/Users/cong/Desktop/barcode/SEQ4416_10112.bam"
bam_file = sys.argv[1]
# fq1 = "/Users/cong/Desktop/barcode/SEQ4416_1_40000.fastq.gz"
fq1 = sys.argv[2]
# fq2 = "/Users/cong/Desktop/barcode/SEQ4416_2_40000.fastq.gz"
fq2 = sys.argv[3]
# output_read = "/Users/cong/Desktop/barcode/SEQ4416_40000_read_count.csv"
output_read = sys.argv[4]
# output_barcode = "/Users/cong/Desktop/barcode/SEQ4416_40000_barcode_count.csv"
output_barcode = sys.argv[5]

barcode_detail_dict, reads_dict, uniq_barcode_dict = main(bam_file=bam_file,fq1=fq1,fq2=fq2)
read_df = pd.DataFrame(reads_dict)
read_df.to_csv(output_read)
barcode_df = pd.DataFrame(uniq_barcode_dict)
barcode_df.to_csv(output_barcode)

# # In[108]:

# print ("barcode details.")
# print ("----------------")
# print ("proper pair: {}").format(barcode_detail_dict['pm']['CDH13_2_mc'])
# print ("unproper pair: {}").format(barcode_detail_dict['mm']['CDH13_2_mc'])
# print ("single map: {}").format(barcode_detail_dict['mu']['CDH13_2_mc'])
# print ("barcode count.")
# print ("----------------")
# print ("proper pair: {}").format(uniq_barcode_dict['pm']['CDH13_2_mc'])
# print ("unproper pair: {}").format(uniq_barcode_dict['mm']['CDH13_2_mc'])
# print ("single map: {}").format(uniq_barcode_dict['mu']['CDH13_2_mc'])
# print ("reads count.")
# print ("----------------")
# print ("proper pair: {}").format(reads_dict['pm']['CDH13_2_mc'])
# print ("unproper pair: {}").format(reads_dict['mm']['CDH13_2_mc'])
# print ("single map: {}").format(reads_dict['mu']['CDH13_2_mc'])


# From the manual check we could find that
# * **TPNB500121:117:H7CLYBGX2:1:11101:25086:1858** and **TPNB500121:117:H7CLYBGX2:1:11101:18196:4543** only have one mismatch on barcode so that they were merged together same barcode
# The results are manually checked correct.
# 

# In[ ]:




# Issues:
# 
# * manual check in small samples and compare with other tools (i.e. feature count) 
# * test one sample in server.
# * discuss w/ groups for next steps.

# In[ ]:




# In[ ]:




# pm: 
# * TPNB500121:117:H7CLYBGX2:1:11101:12199:1424
# * TPNB500121:117:H7CLYBGX2:1:11101:6005:1606
# * TPNB500121:117:H7CLYBGX2:1:11101:25086:1858
# * TPNB500121:117:H7CLYBGX2:1:11101:20480:2932
# * TPNB500121:117:H7CLYBGX2:1:11101:6782:3708
# * TPNB500121:117:H7CLYBGX2:1:11101:25630:4174
# * TPNB500121:117:H7CLYBGX2:1:11101:16199:4357
# * TPNB500121:117:H7CLYBGX2:1:11101:18196:4543
# 
# 
# 
# mu:
# * TPNB500121:117:H7CLYBGX2:1:11101:12066:2614
# * TPNB500121:117:H7CLYBGX2:1:11101:11492:4570
# * TPNB500121:117:H7CLYBGX2:1:11101:11306:1588
# * TPNB500121:117:H7CLYBGX2:1:11101:15236:4264
# * TPNB500121:117:H7CLYBGX2:1:11101:1873:3252
# * TPNB500121:117:H7CLYBGX2:1:11101:4557:5216
# 
# mm:
# * TPNB500121:117:H7CLYBGX2:1:11101:12861:1458 - 145
# * TPNB500121:117:H7CLYBGX2:1:11101:6081:2324
# * TPNB500121:117:H7CLYBGX2:1:11101:3747:2948
# 

# In[ ]:





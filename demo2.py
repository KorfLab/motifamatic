import argparse
import random
import sys
import pwm
import mm
import motif_finder

# Test anrwdistance

#######
# CLI #
#######

parser = argparse.ArgumentParser(
	description='Motif embedding and finding demonstration program')
parser.add_argument('-n', type=int, metavar='<int>', required=False,
	default=50, help = 'number of sequences to generate [%(default)i]')
parser.add_argument('-s', type=int, metavar='<n>', required=False,
	default=60, help = 'size of sequences [%(default)i]')
parser.add_argument('-p', type=float, metavar='<float>', required=False,
	default=0.03, help = 'probability of motif [%(default).2f]')
parser.add_argument('-m', type=str, metavar='<motif>', required=False,
	default='AASTT', help = 'string representation of motif [%(default)s]')
parser.add_argument('-d', action='store_true', help='double-stranded')
parser.add_argument('-r', action='store_true', help='fix random seed')
arg = parser.parse_args()

########
# INIT #
########

if arg.r: random.seed(1)

strand = '+'
if arg.d: strand = '='

################
# Create a PWM #
################
pwm = pwm.string2pwm(arg.m)
print(pwm.pwm_file())

#######################
# Create training set #
#######################
"""
seqs   = [] # synthetic sequences
locs   = [] # locations of synthetic motifs in each sequence (not used yet)
for i in range(arg.n):
	seq, loc = motif_finder.motifembedder(pwm, arg.p, arg.s, strand=strand)
	seqs.append(seq)
	locs.append(loc)
print(seqs)
"""

dtest = ['AGGATCTACGTTCGAACGCATGATACTAGGATAAACCCCGGAGCTTATTATATTAGCAGT', 'AGATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGG', 'ATAGAATGTTGGACTCGCGATCCAAACTTCCCGGCCGACGGATGATGAGTCACGGGATCC']

dtest2 = ['AGGATCTACGTTCGATGAACGCATGATACTAGAAACCCCGGAGCTTATTATATTAGCAGT', 'AGATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGG', 'ATAGAGATATGGATGATTTGGACTCGCGATCCAAACTTCCCGGCCGACGGAGTCACGGCC']

dtest3 = ['AGCTACGTTCGAACGCATACTAGAAACCCCGGAGCTTATTATAGATTTGATAGCGATAGT', 'ATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGATG', 'ATAGAATGTTGGACTCGCCCAAACTTCCCGGCCGACGGATGATGAGTCAGATCGGGATCC']

###########################
# Create background model #
###########################
bkg1 = mm.MM(dtest, order=0)
print(bkg1.mm_file())
bkg2 = mm.MM(dtest2, order=0)
bkg3 = mm.MM(dtest3, order=0)

##########################
# Find kmer-based motifs #
##########################
found_anr = motif_finder.kmer_finder(dtest, bkg1, motif_finder.anr, len(arg.m))
found2_anr = motif_finder.kmer_finder(dtest2, bkg2, motif_finder.anr, len(arg.m))
found3_anr = motif_finder.kmer_finder(dtest3, bkg3, motif_finder.anr, len(arg.m))
found_awd = motif_finder.kmer_finder(dtest, bkg1, motif_finder.anrwdistance, len(arg.m))
found2_awd = motif_finder.kmer_finder(dtest2, bkg2, motif_finder.anrwdistance, len(arg.m))
found3_awd = motif_finder.kmer_finder(dtest3, bkg3, motif_finder.anrwdistance, len(arg.m))
print("anr:")
print(found_anr)
print(found2_anr)
print(found3_anr)
print('Random (kmer):')
print(found_awd)
print('Close to start (kmer):')
print(found2_awd)
print('Far from start (kmer):')
print(found3_awd)

###########################
# Find regex-based motifs #
###########################
found_anr = motif_finder.regex_finder(dtest, bkg1, motif_finder.anr, len(arg.m))
found2_anr = motif_finder.regex_finder(dtest2, bkg2, motif_finder.anr, len(arg.m))
found3_anr = motif_finder.regex_finder(dtest3, bkg3, motif_finder.anr, len(arg.m))
found_awd = motif_finder.regex_finder(dtest, bkg1, motif_finder.anrwdistance, len(arg.m))
found2_awd = motif_finder.regex_finder(dtest2, bkg2, motif_finder.anrwdistance, len(arg.m))
found3_awd = motif_finder.regex_finder(dtest3, bkg3, motif_finder.anrwdistance, len(arg.m))
print("anr:")
print(found_anr)
print(found2_anr)
print(found3_anr)
print('Random (kmer):')
print(found_awd)
print('Close to start (kmer):')
print(found2_awd)
print('Far from start (kmer):')
print(found3_awd)

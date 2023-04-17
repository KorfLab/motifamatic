import argparse
import random
import sys
import pwm
import mm
import motif_finder

# Test distributions with anr (GAT motif)
# Example command line: python3 demo2.py -m GAT -t linear

#######
# CLI #
#######

parser = argparse.ArgumentParser(
	description='Motif finding with anr distributions test program')
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
parser.add_argument('-t', type=str, metavar='<string>', required=False,
	default='none', help = 'distribution of anr: none, linear, quadratic, uniform [%(default)s]')
	
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

# Using GAT as motif
random = ['AGGATCTACGTTCGAACGCATGATACTAGGATAAACCCCGGAGCTTATTATATTAGCAGT', 'AGATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGG', 'ATAGAATGTTGGACTCGCGATCCAAACTTCCCGGCCGACGGATGATGAGTCACGGGATCC']

close = ['AGGATCTACGTTCGATGAACGCATGATACTAGAAACCCCGGAGCTTATTATATTAGCAGT', 'AGATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGG', 'ATAGAGATATGGATGATTTGGACTCGCGATCCAAACTTCCCGGCCGACGGAGTCACGGCC']

far = ['AGCTACGTTCGAACGCATACTAGAAACCCCGGAGCTTATTATAGATTTGATAGCGATAGT', 'ATCCTTGGCTTTGGTAGAGCTTTGCACGCGGCGTTAAAAGGTTGGTGCCCACGTAAGATG', 'ATAGAATGTTGGACTCGCCCAAACTTCCCGGCCGACGGATGATGAGTCAGATCGGGATCC']

###########################
# Create background model #
###########################
bkg1 = mm.MM(random, order=0)
print(bkg1.mm_file())
bkg2 = mm.MM(close, order=0)
bkg3 = mm.MM(far, order=0)

##########################
# Find kmer-based motifs #
##########################
random_anr = motif_finder.kmer_finder(random, bkg1, motif_finder.anr_general, len(arg.m), distr="none")

close_anr = motif_finder.kmer_finder(close, bkg2, motif_finder.anr_general, len(arg.m), distr="none")

far_anr = motif_finder.kmer_finder(far, bkg3, motif_finder.anr_general, len(arg.m), distr="none")

random_distr = motif_finder.kmer_finder(random, bkg1, motif_finder.anr_general, len(arg.m), distr=arg.t)

close_distr = motif_finder.kmer_finder(close, bkg2, motif_finder.anr_general, len(arg.m), distr=arg.t)

far_distr = motif_finder.kmer_finder(far, bkg3, motif_finder.anr_general, len(arg.m), distr=arg.t)

print("anr:")
print(random_anr)
print(close_anr)
print(far_anr)
print('Random (kmer):')
print(random_distr)
print('Close to start (kmer):')
print(close_distr)
print('Far from start (kmer):')
print(far_distr)

###########################
# Find regex-based motifs #
###########################
random_anr = motif_finder.regex_finder(random, bkg1, motif_finder.anr_general, len(arg.m), distr="none")

close_anr = motif_finder.regex_finder(close, bkg2, motif_finder.anr_general, len(arg.m), distr="none")

far_anr = motif_finder.regex_finder(far, bkg3, motif_finder.anr_general, len(arg.m), distr="none")

random_distr = motif_finder.regex_finder(random, bkg1, motif_finder.anr_general, len(arg.m), distr=arg.t)

close_distr = motif_finder.regex_finder(close, bkg2, motif_finder.anr_general, len(arg.m), distr=arg.t)

far_distr = motif_finder.regex_finder(far, bkg3, motif_finder.anr_general, len(arg.m), distr=arg.t)

print("anr:")
print(random_anr)
print(close_anr)
print(far_anr)
print('Random (regex):')
print(random_distr)
print('Close to start (regex):')
print(close_distr)
print('Far from start (regex):')
print(far_distr)

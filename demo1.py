import argparse
import random
import sys

import pwm
import mm
import motif_finder

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
seqs = [] # synthetic sequences
locs = [] # locations of synthetic motifs in each sequence (not used yet)
for i in range(arg.n):
	seq, loc = motif_finder.motifembedder(pwm, arg.p, arg.s, strand=strand)
	seqs.append(seq)
	locs.append(loc)

###########################
# Create background model #
###########################
bkg0 = mm.MM(seqs, order=0)
print(bkg0.mm_file())

##########################
# Find kmer-based motifs #
##########################
found = motif_finder.kmer_finder(seqs, bkg0, motif_finder.anr, len(arg.m))
print(found)

###########################
# Find regex-based motifs #
###########################
found = motif_finder.regex_finder(seqs, bkg0, motif_finder.anr, len(arg.m))
print(found)


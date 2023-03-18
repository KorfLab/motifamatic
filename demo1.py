import argparse
import random
import motiflib


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
	default=0.02, help = 'probability of motif [%(default).2f]')
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

#pwm = motiflib.string2pwm(arg.m)

#UNFINISHED THOUGHTS

#seqs = [] # synthetic sequences
#locs = [] # locations of synthetic motifs in each sequence
#for i in range(arg.n):
#	seq, loc = motiflib.motifembedder(pwm, arg.p, arg.s, strand=strand)
#	seqs.append(seq)
#	locs.append(loc) # not sure we're using these here

train = ['ACGT']
bkgd = motiflib.MM(train, order=0)
#print(bkgd.seq_prob('A'))

rm = motiflib.regex2pwm('A[AC]')
print(rm.pwm_file())

#sm = motiflib.string2pwm('N')
#print(bkgd.pwm_prob(sm))

#s = mm.generate(10, pre='')
#print(s, len(s), s.count('A'))
#print(mm.seq_prob('AAA'))
#print(mm.re_prob('A[AT]A'))




"""
for mname in motif_model:
	finder = motif_model[mname]
	best_motifs = finder(seqs, pwm.length, bkgd, n=10)
	for motif in best_motifs:




pwm model: kmer, regex, dpwm
score model: nseqs, nmotifs
"""


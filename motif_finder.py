import gzip
import math
import itertools
import random
import sys
import re
import pwm
import mm
import tools

#################
# Motif Finding #
#################

def motifembedder(pwm, p, size, choice='ACGT', strand='='):
	seq = ''
	locs = []
	while len(seq) < size:
		if len(seq) > 0 and len(seq) < size - pwm.length and random.random() < p:
			kmer = pwm.generate()
			loc = len(seq)
			if   strand == '=': s = random.choice('+-')
			elif strand == '+': s = '+'
			else:               s = '-'


			if s == '+':
				seq += kmer
				locs.append(loc)
			else:
				seq += tools.anti(kmer)
				locs.append(-loc)
		else:
			seq += random.choice(choice)

	return seq, locs

"""

A scoring function
	requires
		seqs: sequences in case it needs lengths
		locations: where on each sequence the motif is found
		exp: expected probability of motif (given pwm and bkgd model)

	returns
		score


"""

def zoops(seqs, locs, exp):
	n = 0
	x = 0
	for seq, loc in zip(seqs, locs):
		if len(loc) > 0: n += 1
		x += exp * len(seq)

	if n == 0: return 0 # really?

	return math.log2(n / x)


def anr(seqs, locs, exp):
	n = 0
	x = 0
	for seq, loc in zip(seqs, locs):
		n += len(loc)
		x += exp * len(seq)

	if n == 0: return 0 # really?

	return math.log2(n / x)
	
def anrwdistance(seqs, locs, exp):
	n = 0 
	x = 0
	for seq, loc in zip(seqs, locs):
		for l in loc:
			n += (len(seq)-l)/(len(seq))
		x += exp * len(seq)
	
	if n == 0: return 0
	return math.log2(n /x)

# anr without distance
def dn(loc):
	# scoring equation
	score = len(loc)
	return score

# linear distribution for anr with distance
def dl(loc, dist):
	# scoring equation
	if (loc >= dist): score = 0
	else:             score = (dist - loc) / (dist)
	return score
# quadratic distribution for anr with distance
def dq(loc, dist):
	# scoring equation
	if (loc >= dist): score = 0
	else:             score = (1 / dist)(loc + dist)((dist - loc) / dist)
	return score

# uniform distribution for anr with distance
def du(loc, dist):
	# scoring equation
	if (loc >= dist): score = 0
	else:             score = (1 / dist)
	return score

def anr_general(seqs, locs, exp, distr, d):
	n = 0
	x = 0
	# distr is probability distribution for scoring
	# default is 'none' meaning anr without distance 
	if   distr == 'none':      pd = dn      
	elif distr == 'linear':    pd = dl
	elif distr == 'quadratic': pd = dq
	elif distr == 'uniform':   pd = du
	
	for seq, loc in zip(seqs, locs):
		# d is the distance after which scores = 0
		if d == 0: dist = len(seq)
		else:      dist = d
		# scoring with anr no distance
		if distr == "none":
			n += pd(loc)
		# scoring anr with distance
		else: 
			for l in loc:
				n += pd(l, dist)
		x += exp * dist
	if (n == 0): return 0
	return math.log2(n/x)

"""
A motif-finder
	requires
		seqs: a list of sequences of arbitrary length
		bkgd: an n-th order background model
		sfunc: scoring function
		k: some length

	options:
		stranded
		filtering heuristics (entropy, expectation)
		threshold score?

	returns
		motifs with some minimum score?
		best n motifs?

Given a background model, what is the prob of
	kmer: ACG
	regex: A[S]G
	dpwm: AcG
	pwm: complete pwm
	gpwm: gappable pwm


"""

def kmer_finder(seqs, bkgd, func, k, n=10, distr="none", d=0):

	# get all of the kmers present (rather than all possible)
	kmers = {}
	for seq in seqs:
		for i in range(len(seq) -k +1):
			kmers[ seq[i:i+k] ] = True

	# calculate scores of all kmers, and keep the good ones
	keep = []
	for kmer in kmers:
		locs = []
		for seq in seqs:
			pos = []
			off = 0
			while True:
				idx = seq.find(kmer, off)
				if idx == -1: break
				pos.append(idx)
				off = idx + len(kmer)
			locs.append(pos)

		score = func(seqs, locs, bkgd.seq_prob(kmer), distr, d)
		keep.append( (kmer, score) )

		# prevent the list from growing too much
		if len(keep) > 1000:
			keep = sorted(keep, key=lambda t: t[1], reverse=True)
			keep = keep[:n]

	# final sort-n-trim
	keep = sorted(keep, key=lambda t: t[1], reverse=True)
	return keep[:n]

XNT = {
	'A': 0.25,
	'C': 0.25,
	'G': 0.25,
	'T': 0.25,
	'R': 0.50,
	'Y': 0.50,
	'M': 0.50,
	'K': 0.50,
	'W': 0.50,
	'S': 0.50,
	'B': 0.75,
	'D': 0.75,
	'H': 0.75,
	'V': 0.75,
	'N': 1.00,
}

def regex_finder(seqs, bkgd, func, k, n=10, x=0.35, alph='ACGTRYMKWSN',
                 distr="none", d=0):

	keep = []
	for t in itertools.product(alph, repeat=k):
		s = ''.join(t)
		p = 1.0
		for letter in s: p *= XNT[letter]
		if p > x ** len(s): continue

		regex = ''
		for letter in t: regex += tools.NT2RE[letter]

		locs = []
		for seq in seqs:
			pos = []
			for m in re.finditer(regex, seq):
				pos.append(m.span()[0])
			locs.append(pos)

		score = func(seqs, locs, bkgd.re_prob(regex), distr, d)
		keep.append( (regex, score) )

		# prevent the list from growing too much
		if len(keep) > 1000:
			keep = sorted(keep, key=lambda t: t[1], reverse=True)
			keep = keep[:n]

	# final sort-n-trim
	keep = sorted(keep, key=lambda t: t[1], reverse=True)
	return keep[:n]

def dpwm_finder(seqs, bkgd, func, k, n=10, alph='ACGTRYMKWSN'):
	pass


def motiffinder(seqs, k):
	freqs = {}
	total = 0
	for seq in seqs:
		seq = seq.upper()
		for i in range(len(seq)-k+1):
			kmer = seq[i:i+k]
			if kmer not in freqs:
				freqs[kmer] = 0
			freqs[kmer] += 1
	for kmer in freqs:
		print(kmer, freqs[kmer])




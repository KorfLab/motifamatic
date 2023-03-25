import gzip
import math

####################
# Global Constants #
####################

NT2RE = {
	'A': 'A',
	'C': 'C',
	'G': 'G',
	'T': 'T',
	'R': '[AG]',
	'Y': '[CT]',
	'M': '[AC]',
	'K': '[GT]',
	'W': '[AT]',
	'S': '[CG]',
	'B': '[CGT]',
	'D': '[AGT]',
	'H': '[ACT]',
	'V': '[ACG]',
	'N': '[ACGT]',
}

##################
# File Utilities #
##################

def get_filepointer(thing):
	if (type(thing) == str):
		if   thing == '-':          return sys.stdin
		elif thing.endswith('.gz'): return gzip.open(thing, 'rt')
		else:                       return open(thing)
	elif (type(thing)) == io.StringIO: return thing
	elif (type(thing)) == io.TextIOWrapper: return thing
	else: raise ValueError('unknown thing')

def read_fasta(input):

	fp = get_filepointer(input)

	name = None
	seqs = []

	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				yield(name, ''.join(seqs))
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)

	yield(name, ''.join(seqs))
	fp.close()
	
##################
# Math Utilities #
##################

def check_probs(vals, tol=1e-3):
	for val in vals: assert(val >= 0 and val <= 1)
	assert(math.isclose(1.0, sum(vals), abs_tol=tol))

def entropy(vals, check=True):
	if check: check_probs(vals)
	h = 0
	for val in vals:
		if val != 0: h -= val * math.log2(val)
	return h

def dl1(ps, qs, check=True):
	# Manhattan, Taxicab, City Block
	if check:
		check_probs(ps)
		check_probs(qs)

	d = 0
	for p, q, in zip(ps, qs):
		d += abs(p - q)
	return d

def dl2(ps, qs, check=True):
	# Euclidean
	if check:
		check_probs(ps)
		check_probs(qs)

	d = 0
	for p, q in zip(ps, qs):
		d += (p - q) ** 2
	return d ** 0.5

def dkl(ps, qs, check=True):
	# Kullback-Leibler - not recommended
	if check:
		check_probs(ps)
		check_probs(qs)
		for p in ps: assert(p != 0)
		for q in qs: assert(q != 0)

	d = 0
	for p, q in zip(ps, qs):
		d += p * math.log2(p/q)
	return d

######################
# Sequence Utilities #
######################

def anti(seq):
	anti = ''
	for nt in seq[::-1]:
		if   nt == 'A': anti += 'T'
		elif nt == 'C': anti += 'G'
		elif nt == 'G': anti += 'C'
		elif nt == 'T': anti += 'A'
		elif nt == 'R': anti += 'Y'
		elif nt == 'Y': anti += 'R'
		elif nt == 'M': anti += 'K'
		elif nt == 'K': anti += 'M'
		elif nt == 'W': anti += 'W'
		elif nt == 'S': anti += 'S'
		elif nt == 'B': anti += 'V'
		elif nt == 'D': anti += 'H'
		elif nt == 'H': anti += 'D'
		elif nt == 'V': anti += 'B'
		elif nt == 'N': anti += 'N'
		elif nt == 'a': anti += 't'
		elif nt == 'c': anti += 'g'
		elif nt == 'g': anti += 'c'
		elif nt == 't': anti += 'a'
		elif nt == 'r': anti += 'y'
		elif nt == 'y': anti += 'r'
		elif nt == 'm': anti += 'k'
		elif nt == 'k': anti += 'm'
		elif nt == 'w': anti += 'w'
		elif nt == 's': anti += 's'
		elif nt == 'b': anti += 'v'
		elif nt == 'd': anti += 'h'
		elif nt == 'h': anti += 'd'
		elif nt == 'v': anti += 'b'
		elif nt == 'n': anti += 'n'
		else: raise
	return anti

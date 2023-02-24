import gzip
import math
import sys

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

###############
# Motif Class #
###############

class PWM:
	"""Class representing a nucleotide position weight matrix."""

	def __init__(self, init, name=None, source=None):
		"""
		Attributes
		----------
		+ name     `str`   object of type `sequence.DNA`
		+ source   `str`   origin of data, if known
		+ pwm      `list`  a list of `dict` {ACGT}
		+ length   `int`   length of PWM
		+ entropy  `float` sum entropy (or actually 2-H)


		Initialization
		--------------
		+ fasta    `file`  read a fasta file of sequences
		+ seqs     `list`  a list of sequences (strings)
		+ pwm      `list`  a list of dictionaries
		"""

		self.name = name
		self.source = source
		self.pwm = None
		self.length = None
		self.entropy = None

		if   type(init) == str: self._from_fasta(init)
		elif type(init) == list:
			if   type(init[0]) == dict: self._from_pwm(init)
			elif type(init[0]) == str: self._from_seqs(init)
			else: raise ValueError('unknown init type')
		else: raise ValueError('unknown init type')

		# make sure the PWM columns sum close to 1.0
		valid = True
		for c in self.pwm:
			if not math.isclose(1.0, sum(c.values()), abs_tol=1e-3):
				valid = False
				break
		if not valid: raise ValueError(f'sum != 1 {c}')

		self.length = len(self.pwm)
		self.entropy = 0
		for c in self.pwm: self.entropy += entropy(c.values())

	def _from_fasta(self, filename):
		seqs = []
		for name, seq in read_fasta(filename): seqs.append(seq)
		self._from_seqs(seqs)

	def _from_seqs(self, seqs):
		count = []
		for j in range(len(seqs[0])):
			count.append({'A': 0, 'C': 0, 'G': 0, 'T': 0})
		for i in range(len(seqs)):
			for j in range(len(seqs[i])):
				count[j][seqs[i][j]] += 1
		freq = []
		for i in range(len(count)):
			f = {}
			for nt in count[i]: f[nt] = count[i][nt] / len(seqs)
			freq.append(f)
		self.pwm = freq

	def _from_pwm(self, pwm):
		self.pwm = pwm

	def __str__(self):
		lines = []
		lines.append(f'% PWM {self.name} {self.length}')
		nts = 'ACGT'
		for c in self.pwm:
			vals = []
			for nt in c: vals.append(f'{c[nt]:.4f}')
			lines.append(' '.join(vals))
		return '\n'.join(lines)

	def svg(self):

		color = {
			'A': 'fill="red"',
			'C': 'fill="blue"',
			'G': 'fill="green"',
			'T': 'fill="orange"',
		}

		W = len(self.pwm) * 30
		H = 100

		# header
		style = '<style>.sm {font:10px sans-serif;} .lg {font:70px sans-serif;}</style>'
		xmlns = 'xmlns="http://www.w3.org/2000/svg"'
		vbox = f'viewBox="0 0 {W} {H}'

		svg = []
		svg.append(f'<svg {vbox} {xmlns}>')
		svg.append(style)

		# y-axis
		sk = 'stroke="black"'
		svg.append(f'<line x1="{0}" y1="{0}" x2="{0}" y2="{H}" {sk}/>')

		ty = [0, 25, 50, 75, 100]
		x1 = 0
		x2 = 2
		sk = 'stroke="black"'
		for y in ty:
			y1 = y
			y2 = y
			svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')

		# letters
		for i, col in enumerate(self.pwm):
			ys = 2- entropy(col)
			if ys == 0: continue
			yoff = 0
			xoff = i * 46 + 30
			for nt, p in sorted(col.items(), key=lambda item: item[1]):
				if p == 0: continue
				yp = p * ys
				c = color[nt]
				a = 'text-anchor="middle"'
				t = f'transform="scale(1, {yp:.3f})"'
				s = 'class="lg"'
				x = f'x="{xoff}"'
				y0 = H / yp
				yd = 0.5 * yoff * H / yp
				y = f'y="{y0 - yd}"'
				svg.append(f'<g {c} {t}><text {a} {s} {x} {y}>{nt}</text></g>')
				yoff += yp

		# footer
		svg.append('</svg>\n')

		return '\n'.join(svg)

#################################
# Motif Generating Constructors #
#################################

def read_pwm_file(input):

	fp = get_filepointer(input)
	for line in fp:
		if line.startswith('% PWM'):
			pwm = []
			f = line.split()
			name = f[2]
			length = int(f[3])
			for i in range(length):
				line = fp.readline()
				f = line.split()
				A = float(f[0])
				C = float(f[1])
				G = float(f[2])
				T = float(f[3])
				tot = A + C + G + T
				pwm.append({'A': A/tot, 'C': C/tot, 'G': G/tot, 'T': T/tot})
			yield PWM(pwm, name=name)
	fp.close()

def read_transfac(input):

	fp = get_filepointer(input)
	for line in fp:
		# AC = accession number, ID = identifier, NA = name
		if line.startswith('ID'):
			f    = line.split()
			name = f[1]
			pwm  = []
		if line[0].isdigit():
			while line[0].isdigit():
				f = line.split()
				A = float(f[1])
				C = float(f[2])
				G = float(f[3])
				T = float(f[4])
				tot = A + C + G + T
				pwm.append({'A': A/tot, 'C': C/tot, 'G': G/tot, 'T': T/tot})
				line = fp.readline()
			yield PWM(pwm, name=name, source='transfac')
	fp.close()

def _get_count_jaspar(fp):
	line = fp.readline()
	counts = []
	f = line.split()
	for val in f[2:-1]:
		counts.append(float(val))
	return counts

def read_jaspar(input):

	fp = get_filepointer(input)
	while True:
		defline = fp.readline()
		if defline == '': break
		words = defline.split()
		na = _get_count_jaspar(fp)
		nc = _get_count_jaspar(fp)
		ng = _get_count_jaspar(fp)
		nt = _get_count_jaspar(fp)
		pwm = []
		for a, c, g, t in zip(na, nc, ng, nt):
			pwm.append({
				'A': a / (a + c + g + t),
				'C': c / (a + c + g + t),
				'G': g / (a + c + g + t),
				'T': t / (a + c + g + t),
			})
		yield PWM(pwm, name=words[1], source='jaspar')
	fp.close()

###################
# Motif Utilities #
###################

def _align_pwms(m1, m2):
	alignments = []
	if(len(m1.pwm) > len(m2.pwm)):
		for i in range(len(m1.pwm) - len(m2.pwm) + 1):
			pos = zip(m1.pwm[i: i + len(m2.pwm)], m2.pwm)
			alignments.append(pos)
	elif(len(m1.pwm) < len(m2.pwm)):
		for i in range(len(m2.pwm) - len(m1.pwm) + 1):
			pos = zip(m1.pwm, m2.pwm[i: i + len(m1.pwm)])
			alignments.append(pos)
	yield alignments

def motif_distance(m1, m2, method='taxi'):
	# PWMs of same length
	if (len(m1.pwm) == len(m2.pwm)):
		dist = 0
		for pos_m1, pos_m2 in zip(m1.pwm, m2.pwm):
			for nt in pos_m1:
				if (method == 'kl'):
					if (pos_m1[nt] != 0 and pos_m2[nt] != 0):
						dist += pos_m1[nt] * math.log2(pos_m1[nt] / pos_m2[nt])
				else:
						dist += abs(pos_m1[nt] - pos_m2[nt])
	# PWMs of different length
	elif (len(m1.pwm) != len(m2.pwm)):
		dist = 200
		for windows in _align_pwms(m1, m2):
			for i in range(len(windows)):
				d = 0
				for pos_m1, pos_m2 in windows[i]:
					for nt in pos_m1:
						if (method == 'kl'):
							if (pos_m1[nt] != 0 and pos_m2[nt] != 0):
								d += pos_m1[nt] * math.log2(pos_m1[nt]/pos_m2[nt])
						else:
								d += abs(pos_m1[nt] - pos_m2[nt])
				if (d < dist): 	dist = d
	return dist


def cmp_motifs(m1, m2, method='taxi'):

	if   method == 'taxi':   dfunc = dl1
	elif method == 'euclid': dfunc = dl2
	elif method == 'dkl':    dfunc = dkl
	else: raise ValueError('unknown method type')

	if (m1.length < m2.length): (m1, m2) = (m2, m1)

	dmin = None
	for i in range(m1.length - m2.length + 1):
		d = 0
		for c1, c2 in zip(m1.pwm[i:], m2.pwm):
			d += dfunc(c1.values(), c2.values(), check=False)
		if dmin is None or d < dmin: dmin = d

	return dmin

def align(m1, m2, gap=-2):
    match = 3
    scores = [[0.0]*(m1.length+1) for _ in range(m2.length+1)]
    trace = [['-']*(m1.length+1)for _ in range(m2.length+1)]
    for i in range(1, m1.length+1):
        trace[0][i] = 'L'
    for i in range(1, m2.length+1):
        trace[i][0] = 'U'

    maxscore = 0
    maxi = 0
    maxj = 0
    for i in range(1, m1.length + 1):
        for j in range(1, m2.length + 1):
            dist = dl1(m1.pwm[i-1].values(), m2.pwm[j-1].values())
            score = 0.0
            if dist != 2.0:
                score = match*(2-dist)

            left = scores[j-1][i] + gap
            top = scores[j][i-1] + gap
            di = scores[j-1][i-1] + score

            if di < 0 and left < 0 and top <0:
                scores[j][i] = 0.0
            elif di > top and di > left:
                if di > maxscore:
                    maxscore = di
                    maxi = i
                    maxj = j
                scores[j][i] = di
                trace[j][i] = 'D'
            elif top > left:
                scores[j][i] = top
                trace[j][i] = 'U'
            elif left > top:
                scores[j][i] = left
                trace[j][i] = 'L'

    seq = ""
    que = ""
    j = maxj
    i = maxi
    totalscore = 0.0
    while True:
        totalscore += scores[j][i]
        if scores[j][i] == 0:
            break
        if trace[j][i] == 'U':
            seq += f'{i}'
            que += '-'
            i -= 1
        elif trace[j][i] == 'L':
            que += f'{j}'
            seq += '-'
            j -= 1
        elif trace[j][i] == 'D':
            que += f'{j}'
            seq += f'{i}'
            i -= 1
            j -= 1

    #prints out location in sequence rather than nts
    print(seq[::-1])
    print(que[::-1])
    print(f'Score: {totalscore}')
    pass

###########
# Testing #
###########

if __name__ == '__main__':
	import io
	
	fasta_file = """\
>s1
AAATGC
>s2
AAATGC
>s3
ACATGC
>s4
ACATGC
>s5
AGCTGC
>s6
AGCTGC
>s7
ATGTAA
>s8
ATGGAG
"""

	pwm_file = """\
% PWM test1 4
1.0 0.0 0.0 0.0
0.0 1.0 0.0 0.0
0.0 0.0 1.0 0.0
0.0 0.0 0.0 1.0
% PWM test2 3
1.0 0.0 0.0 0.0
0.4 0.3 0.2 0.1
0.1 0.2 0.3 0.4
"""

	transfac_file = """\
AC  M00001
ID  V$MYOD_01
P0      A      C      G      T
01      1      2      2      0      S
02      2      1      2      0      R
03      3      0      1      1      A
04      0      5      0      0      C
//
ID any_old_name_for_motif_1
BF species_name_for_motif_1
P0      A      C      G      T
01      1      2      2      0      S
02      2      1      2      0      R
03      3      0      1      1      A
04      0      5      0      0      C
05      5      0      0      0      A
06      0      0      4      1      G
07      0      1      4      0      G
08      0      0      0      5      T
09      0      0      5      0      G
10      0      1      2      2      K
11      0      2      0      3      Y
12      1      0      3      1      G
//
"""

	jaspar_file = """\
>MA0001.1	AGL3
A  [     0      3     79     40     66     48     65     11     65      0 ]
C  [    94     75      4      3      1      2      5      2      3      3 ]
G  [     1      0      3      4      1      0      5      3     28     88 ]
T  [     2     19     11     50     29     47     22     81      1      6 ]
>MA0002.1	RUNX1
A  [    10     12      4      1      2      2      0      0      0      8     13 ]
C  [     2      2      7      1      0      8      0      0      1      2      2 ]
G  [     3      1      1      0     23      0     26     26      0      0      4 ]
T  [    11     11     14     24      1     16      0      0     25     16      7 ]
"""

	ff = io.StringIO(fasta_file)
	pf = io.StringIO(pwm_file)
	tf = io.StringIO(transfac_file)
	jf = io.StringIO(jaspar_file)
	
	# create a motif from a fasta file
	print('\nFASTA file')
	seqs = [seq for name, seq in read_fasta(ff)]
	m = PWM(seqs, name='testfasta', source='motiflib')
	print(m.name, m.source, m.length, m.entropy)
	print(m)

	# generate motifs objects by reading motif files
	print('\nPWM file')
	for m in read_pwm_file(pf): print(m)
	
	print('\nJASPAR file')
	for m in read_jaspar(jf): print(m)
	
	print('\nTRANSFAC file')
	for m in read_transfac(tf): print(m)

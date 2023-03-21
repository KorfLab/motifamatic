import math
import random
import re
import tools


###################
# PWM/Motif Class #
###################

class PWM:
	"""Class representing a nucleotide position weight matrix."""

	def __init__(self, seqs=None, pwm=None, string=None, name=None, source=None):
		"""
		Initialization (choose one of the options below)
		--------------
		+ string   `str`   an optional motif-like string (e.g. ttNGATYTG)
		+ seqs     `list`  an optional list of sequences (strings)
		+ pwm      `list`  an optional list of dictionaries

		Attributes
		----------
		+ name     `str`   name, if given
		+ source   `str`   origin of data, if known
		+ pwm      `list`  a list of `dict` {ACGT}
		+ length   `int`   length of PWM
		+ entropy  `float` sum entropy (or actually 2-H)

		Methods
		-------
		+ prob(seq)        probability of generating sequence
		+ generate()       an exemplar sequence
		+ string()         a string representation
		+ pwm_file()       a PWM file representation
		+ svg()            an SVG representation
		"""

		assert(seqs is not None or pwm is not None or string is not None)

		self.name = name
		self.source = source
		self.pwm = None
		self.length = None
		self.entropy = None

		# initializers
		if   pwm:    self._from_pwm(pwm)
		elif seqs:   self._from_seqs(seqs)
		elif string: self._from_string(string)
		else: raise ValueError('no initizer')

		# make sure the PWM columns sum close to 1.0
		valid = True
		for c in self.pwm:
			if not math.isclose(1.0, sum(c.values()), abs_tol=1e-3):
				valid = False
				break
		if not valid: raise ValueError(f'sum != 1 {c}')

		self.length = len(self.pwm)
		self.entropy = 0
		for c in self.pwm: self.entropy += tools.entropy(c.values())

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

	def _from_string(self, string):
		p = string2pwm(string)
		self.pwm = p.pwm

	def __str__(self, probs=[]):
		return pwm2string(self, probs=probs)

	def prob(self, seq):
		assert(self.length == len(seq))
		p = 1.0
		for i, nt in enumerate(seq):
			p *= self.pwm[i][nt]
		return p

	def generate(self):
		seq = ""
		for i in range(self.length):
			r = random.random()
			if r < self.pwm[i]['A']:
				seq += "A"
			elif r < self.pwm[i]['A']+self.pwm[i]['C']:
				seq += "C"
			elif r < self.pwm[i]['A'] + self.pwm[i]['C'] + self.pwm[i]['G']:
				seq += "G"
			else:
				seq += "T"
		return seq

	def string(self, probs=[]):
		return pwm2string(self, probs=probs)

	def pwm_file(self):
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
			ys = 2- tools.entropy(col)
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
			yield PWM(pwm=pwm, name=name)
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
			yield PWM(pwm=pwm, name=name, source='transfac')
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
		yield PWM(pwm=pwm, name=words[1], source='jaspar')
	fp.close()

def random_motif(length, name=None, source=None):
	pwm = []
	for i in range(length):
		a = random.random()
		c = random.random()
		g = random.random()
		t = random.random()
		tot = a + c + g + t
		a /= tot
		c /= tot
		g /= tot
		t /= tot
		pwm.append({'A': a, 'C': c, 'G': g, 'T': t})
	return PWM(pwm=pwm, name=name, source=source)

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
                score = match+(1-dist)

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

################################
# Regular Expressions and PWMs #
################################

def regex2pwm(regex, name=None, source=None):
	""" old implementation
	pwm = []
	positions = []
	i = 0
	while (i < len(regex)):
		if (regex[i] == '['):
			x = re.search(']', regex[i:]).start() + i
			positions.append(regex[i+1:x])
			i = x
		elif (regex[i] != '[' and regex[i] != ']'):
			positions.append(regex[i])
		i += 1
	for pos in positions:
		probs = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		for nt in pos:
			if nt not in probs: raise ValueError(f'letter {nt} not allowed')
			p = 1 / len(pos)
			probs[nt] = p
		pwm.append(probs)
	"""

	pwm = []
	pat = '([ACGT])|\[([ACGT]+)\]'
	for m in re.finditer(pat, regex):
		prob = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
		if   m.group(1): nts = m.group(1)
		elif m.group(2): nts = m.group(2)
		else: raise Exception("unexpected letter or pattern")
		for nt in nts: prob[nt] = 1/len(nts)
		pwm.append(prob)
	return PWM(pwm=pwm, name=name, source=source)


def pwm2regex(pwm):
	regex = ''
	for letter in pwm2string(pwm):
		if letter not in tools.NT2RE: regex += '[ACGT]'
		else:                   regex += tools.NT2RE[letter]
	return regex

########################################
# Discretized Nucleotide Probabilities #
########################################

def dnp_table(probs=[0.97, 0.49, 0.33, 0.7, 0.4, 0.3]):
	c1, c2, c3, l1, l2, l3 = probs
	assert(c1 > 0.25)
	assert(c2 < 0.5)
	assert(c3 < 1/3)
	assert(l1 > 0.25)
	assert(l2 < 0.5)
	assert(l3 < 1/3)
	d1 = (1 - c1) / 3
	d2 = 0.5 - c2
	d3 = 1 - c3 * 3
	m1 = (1 - l1) / 3
	m2 = 0.5 - l2
	m3 = 1 - l3 * 3

	return {
		'A': {'A': c1, 'C': d1, 'G': d1, 'T': d1},
		'C': {'A': d1, 'C': c1, 'G': d1, 'T': d1},
		'G': {'A': d1, 'C': d1, 'G': c1, 'T': d1},
		'T': {'A': d1, 'C': d1, 'G': d1, 'T': c1},
		'R': {'A': c2, 'C': d2, 'G': c2, 'T': d2},
		'Y': {'A': d2, 'C': c2, 'G': d2, 'T': c2},
		'M': {'A': c2, 'C': c2, 'G': d2, 'T': d2},
		'K': {'A': d2, 'C': d2, 'G': c2, 'T': c2},
		'W': {'A': c2, 'C': d2, 'G': d2, 'T': c2},
		'S': {'A': d2, 'C': c2, 'G': c2, 'T': d2},
		'B': {'A': d3, 'C': c3, 'G': c3, 'T': c3},
		'D': {'A': c3, 'C': d3, 'G': c3, 'T': c3},
		'H': {'A': c3, 'C': c3, 'G': d3, 'T': c3},
		'V': {'A': c3, 'C': c3, 'G': c3, 'T': d3},
		'a': {'A': l1, 'C': m1, 'G': m1, 'T': m1},
		'c': {'A': m1, 'C': l1, 'G': m1, 'T': m1},
		'g': {'A': m1, 'C': m1, 'G': l1, 'T': m1},
		't': {'A': m1, 'C': m1, 'G': m1, 'T': l1},
		'r': {'A': l2, 'C': m2, 'G': l2, 'T': m2},
		'y': {'A': m2, 'C': l2, 'G': m2, 'T': l2},
		'm': {'A': l2, 'C': l2, 'G': m2, 'T': m2},
		'k': {'A': m2, 'C': m2, 'G': l2, 'T': l2},
		'w': {'A': l2, 'C': m2, 'G': m2, 'T': l2},
		's': {'A': m2, 'C': l2, 'G': l2, 'T': m2},
		'b': {'A': m3, 'C': l3, 'G': l3, 'T': l3},
		'd': {'A': l3, 'C': m3, 'G': l3, 'T': l3},
		'h': {'A': l3, 'C': l3, 'G': m3, 'T': l3},
		'v': {'A': l3, 'C': l3, 'G': l3, 'T': m3},
		'N': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
		'n': {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
	}

def string2pwm(string, probs=[], name=None, source=None):
	if   len(probs) == 6:  t = dnp_table(probs)
	elif len(probs) == 0:  t = dnp_table()
	else: raise ValueError('requires 6 arguments')

	pwm = []
	for nt in string:
		if nt not in t: raise ValueError(f'letter {nt} not allowed')
		pwm.append(t[nt])
	return PWM(pwm=pwm, name=name, source=source)

def pwm2string(pwm, probs=[]):
	if   len(probs) == 6:  t = dnp_table(probs)
	elif len(probs) == 0:  t = dnp_table()
	else: raise ValueError('requires 6 arguments')

	s = ''
	for c in pwm.pwm:
		dmin = None
		best = None
		for nt in t:
			d = dl1(c.values(), t[nt].values())
			if dmin is None or d < dmin:
				dmin = d
				best = nt
		s += best

	return s

#################
# pHMM with PWM #
#################
file_gen = tools.read_fasta(input)

def states(file_gen):
    # list of match states
    bool_states = []
    for line in file_gen:
        temp = []
        for i in line[1]:
            if i != "-": temp.append(True)
            else:temp.append(False)
        bool_states.append(temp)
    marked = []
    for i in range(len(bool_states[0])):
        count = 0
        for j in range(len(bool_states)):
            if bool_states[j][i] == True: count += 1
        if count > 0.5*len(bool_states): marked.append("M")
        else: marked.append("I")
    states = []
    for row in bool_states:
        temp = []
        for ind_col, col in enumerate(row):
            if col == True and marked[ind_col] == "M": temp.append("M")
            elif col == False and marked[ind_col] == "M": temp.append("D")
            elif col == True and marked[ind_col] == "I": temp.append("I")
            else: temp.append("-")
        states.append(temp)
    return states, marked

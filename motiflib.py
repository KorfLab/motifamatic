import gzip
import math
import sys
import json


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
			print(c)
			if not math.isclose(1.0, sum(c.values()), abs_tol=1e-3):
				valid = False
				break
		if not valid: raise ValueError(f'sum != 1 {c}')

		self.length = len(self.pwm)
		self.entropy = 0
		for c in self.pwm: self.entropy += entropy(c)

	def _from_fasta(self, filename):
		seqs = []
		for name, seq in read_fasta(filename): seqs.append(seq)
		self._from_seqs(seqs)
	def _from_seqs(self, seqs):
		self.pwm = seqs2motif(seqs)
	def _from_pwm(self, pwm):
		self.pwm = pwm

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
			print(i, ys, file=sys.stderr)
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


def read_fasta(filename):

	if   filename == '-':          fp = sys.stdin
	elif filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	name = None
	seqs = []

	while True:
		line = fp.readline();
		if line == '': break
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

def seqs2motif(seqs):

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

	return freq

def entropy(col):
	h = 0
	for nt in col:
		if col[nt] != 0: h -= col[nt] * math.log2(col[nt])
	return h



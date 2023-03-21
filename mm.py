import itertools
import re
import random

#######################
# MM/Background Class #
#######################

class MM:
	"""Class representing a Markov model of nucleotide probabilities."""

	def __init__(self, seqs, order=0, pseudo=1.0, name=None):
		"""
		Initialization
		--------------
		+ seqs     `list`  a list of nt sequences (strings)
		+ order    `int`   integer >= 0, default 0
		+ name     `str`   optional name, default None
		+ pseudo   `float` optional pseudocount, default 1.0

		Attributes
		----------
		+ name     `str`
		+ order    `int`
		+ mm       `dict`  [context dict][nt dict] or [dict] if order == 0

		Methods
		-------
		+ seq_prob(seq)    probability of generating sequence
		+ re_prob(re)      probability of generating regex
		+ pwm_prob(pwm)    probability of generating pwm
		+ generate(len)    generate a sequence of some length
		+ mm_file()        file representation
		"""

		self.name = name
		self.order = order
		self.mm = {}

		# init
		counts = {}
		if order == 0:
			counts = {'A':pseudo, 'C':pseudo, 'G':pseudo, 'T':pseudo}
		else:
			for nts in itertools.product('ACGT', repeat=order):
				ctx = ''.join(nts)
				counts[ctx] = {'A':pseudo, 'C':pseudo, 'G':pseudo, 'T':pseudo}

		# add counts
		for seq in seqs:
			for i in range(order, len(seq) - order):
				if order == 0:
					nt = seq[i]
					counts[nt] += 1
				else:
					ctx = seq[i-order:i]
					nt = seq[i]
					counts[ctx][nt] += 1

		# assign frequencies
		if order == 0:
			total = sum(counts.values())
			for nt in counts: self.mm[nt] = counts[nt] / total
		else:
			for ctx in counts:
				total = sum(counts[ctx].values())
				if total == 0:
					raise Exception("not enough observations for context")
				self.mm[ctx] = {}
				for nt in counts[ctx]:
					self.mm[ctx][nt] = counts[ctx][nt] / total

	def seq_prob(self, seq):
		if self.order == 0:
			p = 1.0
			for nt in seq: p *= self.mm[nt]
			return p

		assert(len(seq) > self.order)
		p = 1.0
		for i in range(self.order, len(seq)):
			ctx = seq[i-self.order:i]
			nt = seq[i]
			p *= self.mm[ctx][nt]
		return p

	def re_prob(self, regex):
		# there must be a more efficient way to do this

		# create all the sequences from the regex
		seqs = ['']
		pat = '([ACGT])|\[([ACGT]+)\]'
		for m in re.finditer(pat, regex):
			if   m.group(1): nts = m.group(1)
			elif m.group(2): nts = m.group(2)
			else: raise Exception("unexpected letter or pattern")
			newseqs = []
			for seq in seqs:
				for nt in nts:
					newseqs.append(seq + nt)
			seqs = newseqs

		# sum the expected sequence probabilities
		p = 0.0
		for seq in seqs: p += self.seq_prob(seq)

		return p

	def pwm_prob(self, pwm, x=0.5):
		# this is brute force, there must be a better way

		threshold = x ** pwm.length
		n = 0
		for nts in itertools.product('ACGT', repeat=pwm.length):
			seq = ''.join(nts)
			p = pwm.prob(seq)
		#	q = self.seq_prob(seq)
		#	print(seq, p, q)
			if p > threshold: n += 1

		return n / 4 ** pwm.length

	def generate(self, n, pre='', marg=[0.25, 0.25, 0.25, 0.25]):
		if self.order == 0:
			s = ''
			for i in range(n):
				s += random.choices('ACGT', weights=self.mm.values())[0]
			return s

		s = pre
		while len(s) < self.order:
			s += random.choices('ACGT', weights=marg)[0]

		for i in range(len(s), n):
			ctx = s[i-self.order:i]
			s += random.choices('ACGT', weights=self.mm[ctx].values())[0]
		return s

	def mm_file(self):
		lines = []
		lines.append(f'% MM {self.name} {4**(self.order+1)}')
		if self.order == 0:
			for nt in 'ACGT':
				lines.append(f'{nt} {self.mm[nt]:.6f}')
		else:
			for ctx in self.mm:
				for nt in 'ACGT':
					lines.append(f'{ctx}{nt} {self.mm[ctx][nt]:.6f}')
				lines.append('')
		return '\n'.join(lines)

import pwm
import re 
import itertools
import motif_finder

string2re = {
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

def pwm2regex(pwm):
	regex = ''
	for letter in pwm.pwm2string(pwm):
		if letter not in string2re: regex += '[ACGT]'
		else: 						regex += string2re[letter]
	return regex

pwm = pwm.string2pwm('ACRTT')
#print(pwm.pwm_file())
#print(pwm.string())
#print(pwm)

#print(pwm2regex(pwm))

seqs = [seq for seq in motif_finder.motifembedder(pwm, 0.05, 50)]
print(seqs)
# expected probabilities:
exp = {
	'A': 0.25,
	'C': 0.25,
	'G': 0.25,
	'T': 0.25,
	'R': 0.5,
	'Y': 0.5,
	'M': 0.5,
	'K': 0.5,
	'W': 0.5,
	'S': 0.5,
	'B': 0.75,
	'D': 0.75,
	'H': 0.75,
	'V': 0.75,
	'N': 1.0,
}
# find patterns using regex and associated probabilities
for t in itertools.product('ACGTRYMKWSN', repeat=5):
	s = ''.join(t)
	p = 1
	for letter in s: 
		p *= exp[letter]
	if p > 0.3 ** len(s): continue
	regex = ''
	for letter in t: 
		regex += string2re[letter]
	count = 0
	for seq, pos in seqs:
		seq = seq.upper()
		#print(regex, seq)
		count += len(re.findall(regex, seq))
	print(regex, count, p, count/p)
		#for m in re.finditer(regex, seq): 
		#	print(m.span(), m.group())
	



import gzip
import math

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

def yscale(col):
	h = 0
	for nt in col:
		if col[nt] != 0: h -= col[nt] * math.log2(col[nt])
	return 2 - h


def motif2svg(motif):

	x = 100 + len(motif) * 20
	y = 100

	style = '<style>.sm {font:10px sans-serif;} .lg {font:66px sans-serif;}</style>'
	xmlns = 'xmlns="http://www.w3.org/2000/svg"'
	vbox = f'viewBox="0 0 {x} {y}'
	
	svg = []
	svg.append(f'<svg {vbox} {xmlns}>')
	svg.append(style)
	
	# y-axis
	x1 = 0
	y1 = 0
	x2 = 0
	y2 = 100
	sk = 'stroke="black"'
	svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')

	
	
	# x-axis
	x1 = -1
	y1 = 101
	x2 = 100
	y2 = 101
	sk = 'stroke="black"'
	svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')
	
	
	"""
	# bits label
	x1 = -68
	y1 = -15
	svg.append(f'<g transform="rotate(-90)"> <text x="{x1}" y="{y1}" class="sm">Bits</text> </g>')
	
	# bits ticks and values
	tn = ['2.0', '1.5', '1.0', '0.5', '0.0']
	ty = [0, 25, 50, 75, 100]
	x1 = gt -3
	x2 = gt
	sk = 'stroke="black"'
	for n, y in zip(tn, ty):
		y1 = y + gt
		y2 = y + gt
		svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')
		
	for n, y in zip(tn, ty):
		y1 = y + gt + 3
		x1 = -8
		svg.append(f'<text x="{x1}" y="{y1}" class="sm">{n}</text>')
	
	for col in motif:
		ys = yscale(col)
		
	
	#svg.append('<g fill="red" transform="scale(1.0, 1.0)"> <text x="10" y="109" class="lg">ACGT</text> </g>')
	svg.append('<g fill="red"> <text x="10" y="109" class="lg">ACGT</text> </g>')
	#svg.append('<g fill="red" transform="scale(1, 0.5)"> <text x="20" y="35" class="lg">ACGT</text> </g>')
	#svg.append('<g fill="red" transform="scale(1, 1.0)"> <text x="20" y="35" class="lg">ACGT</text> </g>')
	"""
	
	svg.append('</svg>\n')
		
	return '\n'.join(svg)
import gzip
import math
import sys
import json

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

	color = {
		'A': 'fill="red"',
		'C': 'fill="blue"',
		'G': 'fill="green"',
		'T': 'fill="orange"',
	}

	G = 20
	W = G + len(motif) * 30
	H = 100
	
	# header
	style = '<style>.sm {font:10px sans-serif;} .lg {font:70px sans-serif;}</style>'
	xmlns = 'xmlns="http://www.w3.org/2000/svg"'
	vbox = f'viewBox="0 0 {W} {H}'
	
	svg = []
	svg.append(f'<svg {vbox} {xmlns}>')
	svg.append(style)
	
	# y-axis
	x1 = G -2
	y1 = G
	x2 = G -2
	y2 = H + G
	sk = 'stroke="black"'
	svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')
	
	tn = ['2.0', '1.5', '1.0', '0.5', '0.0']
	ty = [0, 25, 50, 75, 100]
	x1 = G -5
	x2 = G -2
	sk = 'stroke="black"'
	for n, y in zip(tn, ty):
		y1 = y + G
		y2 = y + G
		svg.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" {sk}/>')
	
	for n, y in zip(tn, ty):
		y1 = y + 3 + G
		x1 = 0
		svg.append(f'<text x="{x1}" y="{y1}" class="sm">{n}</text>')
	
	# letters
	for i, col in enumerate(motif):
		ys = yscale(col)
		yoff = 0 # percent of pixels already used
		xoff = i * 45 + G + 25
		for nt, p in sorted(col.items(), key=lambda item: item[1]):
			if p == 0: continue			
			yp = p * ys # this nt proportion of scale
			c = color[nt]
			a = 'text-anchor="middle"'
			t = f'transform="scale(1, {yp:.3f})"'
			s = 'class="lg"'
			x = f'x="{xoff}"'
			y0 = (H+G) / (p * ys) # starting position
			yd = H * yoff * p / yp # y-delta, not quite correct
			y = f'y="{y0-yd}"'
			
			print(ys, nt, p, file=sys.stderr)
			
			svg.append(f'<g {c} {t}><text {a} {s} {x} {y}>{nt}</text></g>')
			yoff += p
	
	# footer
	svg.append('</svg>\n')
	
	return '\n'.join(svg)
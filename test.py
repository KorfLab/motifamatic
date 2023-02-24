import sys
import motiflib
import json


m = motiflib.PWM('test.fasta', name='from-fasta', source='fabricated')
print(m.name)
print(m.source)
print(m.length)
print(m.entropy)
print(m)
print(json.dumps(m.pwm, indent=4))

motifs = []
for m in motiflib.read_pwm_file('test.pwm'):
	print(m)
	motifs.append(m)

for m in motiflib.read_transfac('test.transfac'):
	print(m)

for m in motiflib.read_jaspar('JASPAR2022_CORE_redundant_pfms_jaspar.txt.gz'):
	print(m)
	break

print(motiflib.cmp_motifs(motifs[0], motifs[1], method='taxi'))
print(motiflib.cmp_motifs(motifs[0], motifs[1], method='euclid'))
try:
	print(motiflib.cmp_motifs(motifs[0], motifs[1], method='dkl'))
except:
	raise

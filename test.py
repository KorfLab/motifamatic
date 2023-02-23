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

print('----')
print(motifs[0])
print(motifs[1])

print(motiflib.cmp_motifs(motifs[1], motifs[0], method='dkl'))
print(motiflib.motif_distance(motifs[1], motifs[0], method='kl'))
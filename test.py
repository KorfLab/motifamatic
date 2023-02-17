import sys
import motiflib
import json

m = motiflib.PWM('test.fasta', name='from-fasta', source='fabricated')
print(m.name)
print(m.source)
print(m.length)
print(m.entropy)
print(m)

for m in motiflib.read_pwm_file('test.pwm'):
	print(m)
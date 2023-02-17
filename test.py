import sys
import motiflib
import json

m = motiflib.PWM('motif.fasta', name='test', source='fabricated')
print(m.name)
print(m.source)
print(m.length)
print(m.entropy)
print(m)

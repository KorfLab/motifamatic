import math
import pytest
import pwm

def test_random_motif():
	m = pwm.random_motif(5, name='random', source='random')
	assert(m.name == 'random')
	assert(m.source == 'random')
	assert(m.length == 5)

def test_from_seqs():
	seqs = [seq for name, seq in motiflib.read_fasta('tests/example.fasta')]
	m = pwm.PWM(seqs=seqs)
	assert(math.isclose(m.pwm[0]['A'], 1.0))
	assert(math.isclose(m.pwm[1]['A'], 0.25))

def test_from_pwm():
	pwm = [
		{'A': 1.00, 'C': 0.00, 'G': 0.00, 'T': 0.00},
		{'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
	]
	m = motiflib.PWM(pwm=pwm)
	assert(math.isclose(m.pwm[0]['A'], 1.0))
	assert(math.isclose(m.pwm[1]['A'], 0.25))

def test_error_from_nothing():
	with pytest.raises(AssertionError):
		m = pwm.PWM()

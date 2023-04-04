import pwm
import pytest
import math
import regex

def test_regex2pwm_ok():
	motif = '[ACGT]A'
	a = pwm.regex2pwm(motif)
	assert(math.isclose(a.pwm[0]['A'], 0.25))
	assert(math.isclose(a.pwm[0]['C'], 0.25))
	assert(math.isclose(a.pwm[0]['G'], 0.25))
	assert(math.isclose(a.pwm[0]['T'], 0.25))
	assert(math.isclose(a.pwm[1]['A'], 1.0))
	assert(math.isclose(a.pwm[1]['C'], 0.0))
	assert(math.isclose(a.pwm[1]['G'], 0.0))
	assert(math.isclose(a.pwm[1]['T'], 0.0))

def test_regex2pwm_invalid_nt():
	motif = 'BA[CG]'
	with pytest.raises(ValueError):
		pwm.regex2pwm(motif)

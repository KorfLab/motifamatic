import math
import pytest
import sys
import pwm

myprobs = [0.97, 0.49, 0.33, 0.70, 0.40, 0.30]

def test_string2pwm():

	m = pwm.string2pwm('ACGTRYMKWSBDHVNacgtrymkwsbdhvn', probs=myprobs)
	assert(math.isclose(m.pwm[0]['A'], 0.97))
	assert(math.isclose(m.pwm[1]['A'], 0.01))
	assert(math.isclose(m.pwm[2]['A'], 0.01))
	assert(math.isclose(m.pwm[3]['A'], 0.01))
	assert(math.isclose(m.pwm[4]['A'], 0.49))
	assert(math.isclose(m.pwm[5]['A'], 0.01))
	assert(math.isclose(m.pwm[6]['A'], 0.49))
	assert(math.isclose(m.pwm[7]['A'], 0.01))
	assert(math.isclose(m.pwm[8]['A'], 0.49))
	assert(math.isclose(m.pwm[9]['A'], 0.01))
	assert(math.isclose(m.pwm[10]['A'], 0.01))
	assert(math.isclose(m.pwm[11]['A'], 0.33))
	assert(math.isclose(m.pwm[12]['A'], 0.33))
	assert(math.isclose(m.pwm[13]['A'], 0.33))
	assert(math.isclose(m.pwm[14]['A'], 0.25))
	assert(math.isclose(m.pwm[15]['A'], 0.70))
	assert(math.isclose(m.pwm[16]['A'], 0.10))
	assert(math.isclose(m.pwm[17]['A'], 0.10))
	assert(math.isclose(m.pwm[18]['A'], 0.10))
	assert(math.isclose(m.pwm[19]['A'], 0.40))
	assert(math.isclose(m.pwm[20]['A'], 0.10))
	assert(math.isclose(m.pwm[21]['A'], 0.40))
	assert(math.isclose(m.pwm[22]['A'], 0.10))
	assert(math.isclose(m.pwm[23]['A'], 0.40))
	assert(math.isclose(m.pwm[24]['A'], 0.10))
	assert(math.isclose(m.pwm[25]['A'], 0.10))
	assert(math.isclose(m.pwm[26]['A'], 0.30))
	assert(math.isclose(m.pwm[27]['A'], 0.30))
	assert(math.isclose(m.pwm[28]['A'], 0.30))
	assert(math.isclose(m.pwm[29]['A'], 0.25))

def test_pwm2string():
	pwm1 = [
		{'A': 0.97, 'C': 0.01, 'G': 0.01, 'T': 0.01},
		{'A': 0.94, 'C': 0.02, 'G': 0.02, 'T': 0.02},
		{'A': 0.91, 'C': 0.03, 'G': 0.03, 'T': 0.03},
		{'A': 0.88, 'C': 0.04, 'G': 0.04, 'T': 0.04},
		{'A': 0.85, 'C': 0.05, 'G': 0.05, 'T': 0.05},
		{'A': 0.82, 'C': 0.06, 'G': 0.06, 'T': 0.06},
		{'A': 0.79, 'C': 0.07, 'G': 0.07, 'T': 0.07},
		{'A': 0.76, 'C': 0.08, 'G': 0.08, 'T': 0.08},
		{'A': 0.73, 'C': 0.09, 'G': 0.09, 'T': 0.09},
		{'A': 0.70, 'C': 0.10, 'G': 0.10, 'T': 0.10},
	]
	m = pwm.PWM(pwm=pwm1)
	s = pwm.pwm2string(m, probs=myprobs)
	assert(s == 'AAAAAaaaaa')


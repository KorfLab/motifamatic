import pytest
import math

import motiflib

def test_checkprobs_negative():
	p = [1.0, -1.0]
	with pytest.raises(AssertionError):
		motiflib.check_probs(p)

def test_checkprobs_sum():
	p = [0.5, 0.6]
	with pytest.raises(AssertionError):
		motiflib.check_probs(p)

def test_entropy1():
	p = [0.5, 0.5]
	assert(math.isclose(1.0, motiflib.entropy(p)))

def test_entropy2():
	p = [0.25, 0.25, 0.25, 0.2500001] # does not quite sum to 1, ok
	assert(math.isclose(2.0, motiflib.entropy(p), abs_tol=1e-3))

"""
Tests to add
dl1() Jenna
dl2() Sai
dkl() Claire
"""

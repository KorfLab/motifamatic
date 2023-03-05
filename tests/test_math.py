import pytest
import math
import sys

import motiflib

def test_check_probs_ok():
	p = [0.4, 0.3, 0.2, 0.1]
	motiflib.check_probs(p)

def test_checkprobs_negative():
	p = [1.0, -1.0]
	with pytest.raises(AssertionError):
		motiflib.check_probs(p)

def test_checkprobs_sum():
	p = [0.5, 0.6]
	with pytest.raises(AssertionError):
		motiflib.check_probs(p)

def test_dl1_ok():
	ps = [0.4, 0.3, 0.2, 0.1]
	qs = [0.1, 0.2, 0.3, 0.4]
	assert(math.isclose(motiflib.dl1(ps, qs), 0.8))

def test_dl1_empty():
	ps = []
	qs = []
	with pytest.raises(AssertionError):
		motiflib.dl1(ps, qs)

def test_dl1_negative():
	ps = [1, -1, 1]
	qs = [-0.5, 0.5, 1]
	with pytest.raises(AssertionError):
		motiflib.dl1(ps, qs)

def test_dl1_sum():
	ps = [0.5, 0.2, 0.1]
	qs = [0.5, 0.6]
	with pytest.raises(AssertionError):
		motiflib.dl1(ps, qs)

def test_entropy():
	ps = [0.25, 0.25, 0.25, 0.25]
	assert(motiflib.entropy(ps), 2)

def test_dkl_ok():
	ps = [0.2, 0.4, 0.1, 0.3]
	qs = [0.2, 0.4, 0.1, 0.3] 
	assert(math.isclose(motiflib.dkl(ps, qs), 0))

def test_dkl_zeroes():
	ps = [0.4, 0.6, 0.3, 0.0]
	qs = [0.1, 0.3, 0.3, 0.3]
	with pytest.raises(AssertionError):
		motiflib.dkl(ps, qs)

def test_dkl_empty():
	ps = [0.2, 0.4, 0.4]
	qs = []
	with pytest.raises(AssertionError):
		motiflib.dkl(ps, qs)
		
def test_dkl_negative():
	ps = [0.8, 0.1, 0.1]
	qs = [0.2, -0.8, 0.8]
	with pytest.raises(AssertionError):
		motiflib.dkl(ps, qs)

def test_dkl_sum():
	ps = [0.4, 0.2]
	qs = [0.1, 0.1, 0.1, 0.1]
	with pytest.raises(AssertionError):
		motiflib.dkl(ps, qs)
	
"""
entropy() has no tests: Claire
dl2() has no tests: Meghana
dkl() has no tests: Sai
"""

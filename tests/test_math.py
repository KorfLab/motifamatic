import pytest
import math

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

"""
entropy() has no tests: Claire
dl1() has no tests: Jenna
dl2() has no tests: Meghana
dkl() has no tests: Sai
"""

import io
import pytest
import sys

import motiflib

def test_filepointer_stdin():
	fp = motiflib.get_filepointer('-')
	assert(fp is sys.stdin)

def test_filepointer_error():
	with pytest.raises(ValueError):
		fp = motiflib.get_filepointer([]) # other thing

def test_fasta_file():
	seqs = [seq for name, seq in motiflib.read_fasta('tests/example.fasta')]
	assert(len(seqs) == 8)

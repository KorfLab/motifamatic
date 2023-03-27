import io
import pytest
import sys

import tools

def test_filepointer_stdin():
	fp = tools.get_filepointer('-')
	assert(fp is sys.stdin)

def test_filepointer_error():
	with pytest.raises(ValueError):
		fp = tools.get_filepointer([]) # other thing

def test_fasta_file():
	seqs = [seq for name, seq in tools.read_fasta('tests/example.fasta')]
	assert(len(seqs) == 8)


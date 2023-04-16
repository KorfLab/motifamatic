import math
import pytest
import pwm
import tools
import mm
import motif_finder

def test_random_motif():
	m = pwm.random_motif(5, name='random', source='random')
	assert(m.name == 'random')
	assert(m.source == 'random')
	assert(m.length == 5)

def test_from_seqs():
	seqs = [seq for name, seq in tools.read_fasta('tests/example.fasta')]
	m = pwm.PWM(seqs=seqs)
	assert(math.isclose(m.pwm[0]['A'], 1.0))
	assert(math.isclose(m.pwm[1]['A'], 0.25))

def test_from_pwm():
	pwm1 = [
		{'A': 1.00, 'C': 0.00, 'G': 0.00, 'T': 0.00},
		{'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25},
	]
	m = pwm.PWM(pwm=pwm1)
	assert(math.isclose(m.pwm[0]['A'], 1.0))
	assert(math.isclose(m.pwm[1]['A'], 0.25))

def test_error_from_nothing():
	with pytest.raises(AssertionError):
		m = pwm.PWM()
	
def test_anrwdistance():
	beg_motif = ['ATAACTTTGTTCAACTTAAAATGAAGCTCCAGTAGGGAGCACGACAGAGATGTCCAGATTGAGTCAGCGT', 'CTGGTGACAACTTGGAAGGATTAAAGAGGAGAATTTAAAAAAATTAAGTTGTCAAGTGTCACGTGGCGAC','AAAAAACTTTTTATCTGCAGCCCTTGCTTAAATTTGTGGTAGACAAATGGCTTCTAAATTTAATTTTAGC']
	end_motif = ['ATTGTTCAACTTAAAATGAAGCTCCAGTAGGGTCCAGATTGAGTAACTTCAGATCTGCGT', 'CTGGTGACGGAAGGATTAAAGAGGAGAATTTACAAGTGTCACGTGGCGGATTTAAACTTC','AAAATTTATCTGCAGCCCTTGCTTAAATTTGTTCTAACTAACTTTTAACTTGATCTTAGC']
	bkg0 = mm.MM(beg_motif, order=0)
	bkg1 = mm.MM(end_motif, order=0)
	found1 = motif_finder.kmer_finder(beg_motif, bkg0, motif_finder.anrwdistance, 5)
	found2 = motif_finder.kmer_finder(end_motif, bkg1, motif_finder.anr, 5)
	motif_score = 0
	for i in found1:
		if i[0] == "GATCT":
			motif_score = i[1]
	for j in found2:
		if j[0] == "GATCT":
			assert(motif_score > i[1])

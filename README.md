motifamatic
===========

Software for finding, comparing, and displaying nucleotide motifs (position
weight matrices).

## Intent ##

While there are excellent motif-finding software packages that already exist
(e.g. MEME suite), they aren't optimized for some of our needs (e.g. IME
signals). Motifs are also a simple sequence model that make a good playground
for learning bioinformatics concepts.

## Motif-finding ##

Given a collection of related sequences, find the sub-sequence pattern with the
most surprising distribution. For example, in the following sequences, `ACGT`
is found multiple times (shown in uppercase).

```
ttgaacACGTtgtgcttcaaacggcctatgtaagtgcggtttgg
tagcttcctcagagtgccgctgcgtcacgaACGTctaggttgac
tagtccctatctactgtatagtgtccctcggaacctacgccact
acaagccACGTtgtgcaaccgaattcactcctgACGTagcaggc
gggaatcaacctaccgtgACGTccaagtgggtcaatcgccggcc
gttgcctgcatggtatctccttgtcggctgtaatacgaagaaca
```

Is it likely for `ACGT` to be present in 4/6 sequences of length 40? Is it
likely for there to be 5 copies of `ACGT` among these sequences? Depending on
your probability model, the distribution of `ACGT` might be highly unexpected
and therefore worth reporting.

## Position Weight Matrices ##

Nucleotide patterns are generally not 100% conserved. That is, they are rarely
just `ACGT` but more like `ACG` followed sometimes by `C` and sometimes by `T`.
There are several ways to represent such ambiguities in text. Biologists might
write `ACG(C/T)` while those more on the computational spectrum might follow
regex sytax and write `ACG[CT]`. The more proper way is to use IUPAC ambiguity
symbols: `ACGY`. While text representations are convenient, they show
possiblities rather than probabilities. For example, what if `C` is much more
common than `T`?

A Position Weight Matrix (PWM) specifies the probabilities of nucleotides at
each position of a motif of length n. For example, below is a PWM describing
splice donor sites in C. elegans. The site is 5 nt long and displays the
probabilties of A, C, G, and T in a column at each position. The most common
sequence generated by this PWM is GTAAG.

```
     1      2      3      4      5
A  0.000  0.000  0.586  0.694  0.092
C  0.000  0.000  0.014  0.056  0.035
G  1.000  0.000  0.254  0.102  0.766
T  0.000  1.000  0.147  0.148  0.106
```

Representing this PWM in text is a matter of preference. The following all
capture the pattern to various degrees. Note the use of lowercase letters to
represent a weak single nucleotide preference.

+ GTAAG - only canonical letters
+ GTRAG - noting the preference for both A and G at the 3rd position
+ GTaag - using lowercase letters instead of IUPAC symbols
+ GTrag - maybe the best representation
+ GTNNN - conflating possibility and probability

## Progression ##

+ kmers
+ regular expressions
+ discretized position weight matrices
+ position weight matrices

The simplest model for a motif is that it does not allow for any ambiguities.
For the splice donor above, this is simply `GTAAG`. Even with this simple
model, we can evaluate the probability of finding the pattern in a collection
of sequences. Therefore, the first motif-finding algorithm is to find the most
surprising k-mers. In a collection of splice donor sites, we would expect
`GTAAG` to be slightly more common than `GTGAG`.

Regular expressions offer a convenient method to find ambiguities that are
equally likely. For example `GT[AG]AG` is one way to describe the motif.
However, is this better than `GT[ACGT]AG`? It depends on the probability model
and collection of sequences.

Another way to represent PWMs is to discretize the probability space and assign
it to a letter. For example, `A` doesn't mean 100% `A`, and `R` doesn't mean
50% `A` and 50% `G`. Instead, we can set `A` to 97% and give the other
nucleotides 1%. We can do related transformations with ambiguity codes and
lowercase letters.

Finally, we can also represent the individual nucleotide probabilities as
floating point values of arbitrary precision.

## Internal Representations ##

A motif is an object of type `PWM`. The `PWM` class defines the following
attributes:

+ motif.name - hopefully unique
+ motif.source - could be TRANSFAC, JASPAR, etc.
+ motif.pwm - a list of dictionaries containing nucleotide probabilities
+ motif.length
+ motif.entropy - the sum of individual entropies (2 - entropy really)

The internal representation of motif.pwm is as follows:

```
[
    {
        "A": 1.0,
        "C": 0.0,
        "G": 0.0,
        "T": 0.0
    },
    {
        "A": 0.5,
        "C": 0.5,
        "G": 0.0,
        "T": 0.0
    },
    {
        "A": 0.0,
        "C": 0.0,
        "G": 0.75,
        "T": 0.25
    },
    {
        "A": 0.25,
        "C": 0.25,
        "G": 0.25,
        "T": 0.25
    }
]
```

There are several ways to construct motifs:

+ Random
+ From a case-sensitive degenerate nucleotide string
+ From a list of sequences - all equal length
+ From a list of dictionaries of probabilities
+ From reading JASPAR, TRANSFAC, or our own formatted files

The `motiflib.py` file has embedded code examples which can be run as a form of
testing via `python3 motiflib.py`.

## To Do List ##

Coding stuff

+ Class representing PWMs - done
+ Constructor from sequences - done
+ Random constructor - done
+ Read from our own PWM model file format - done
+ Read from TRANSFAC format - done
+ Read from JASPAR format - done
+ Compare ungapped motifs - done
+ Compare gapped motifs - WIP
+ Convert ambiguity strings to motifs - done
+ Convert motifs to ambiguity strings - done
+ Convert regex to motifs -
+ Convert motifs to regex -
+ Display motifs as SVG - WIP
+ Make a motif-finder based on strings - WIP
+ Make a motif-finder based on regex -
+ Make a motif-finder based on discretized PWMs -
+ Make a motif-finder based on full PWMs -
+ Unit testing - WIP

Data stuff

+ Create test sequences for motif-finding - WIP
	+ motif-embedding function - WIP
		+ background probability model rather than 25% each
		+ sequences must be exactly the right length
		+ coordinate system of motifs is off
		+ reverse-complement motifs should be an option
+ Get real sequences for motif-finding

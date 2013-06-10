# A simple example with optimizer.
# This demonstrates a new setting in version 2.1 of PHASE, which
# allows equilibrium frequencies to be estimated from the
# base frequencies in the alignment, rather than being estimated
# by maximum likelihood.

# See the manual and other examples for more advanced examples.

{DATAFILE}
Data file = sequence-data/s10-TN93dG10.dna
Interleaved data file = no
{\DATAFILE}

{MODEL}
Model = TN93
Discrete gamma distribution of rates = yes
Number of gamma categories = 10
Invariant sites = no
Optimize model parameters = yes
Empirical frequencies = yes
{\MODEL}

{TREE}
Outgroup = o1
Tree file = sequence-data/s10.tre
{\TREE}

Random seed = 29072011

Output file = results/optimizer.s10-TN93dG10


# The program distphase generates a matrix of pairwise ML distances for a
# given set of sequences (base-pair models can be used with RNA stems).
# This matrix can be used with a tree reconstruction algorithm
# (UGPMA/NJ), not available in the PHASE package.

# The drawback of this program is that the parameters of the substitution
# model must be specified, and thus already need to have been found by
# another method.

# A standard DATAFILE block (see manual).
{DATAFILE}
Data file = sequence-data/s10-TN93dG10.dna
Interleaved data file = no
Heterogeneous data models = no
{\DATAFILE}

# A standard MODEL block (see manual) which also contains the name of the "model file"
{MODEL}
Model = TN93
Discrete gamma distribution of rates = yes
Number of gamma categories = 10
Invariant sites = no

# Model parameters are required with distphase, in the form of a "model file".
# A skeleton for this file, which can then be filled out with user-defined
# parameters, can be produced with the 'simulate' program; or the file can be
# generated with optimized parameters by using mlphase or optimizer.
Model parameters file = sequence-data/s10-TN93dG10.mod
{\MODEL}

Random seed = 29072011

Output file = results/distphase.s10-TN93dG10.out
#Output format = lower-triangular


# Searching for the optimal model and branch lengths in the
# ML framework with optimizer.
# The likelihoods of a set of candidate tree topologies are compared.

# A standard DATAFILE block (see manual).
{DATAFILE}
Data file = sequence-data/hiv6.dna
Interleaved data file = yes
{\DATAFILE}

# A standard MODEL block (see manual), with an option to optimize or
# use user-defined substitution parameters.
{MODEL}
Model = TN93
Discrete gamma distribution of rates = no
Invariant sites = yes

# If this value is 'no', you must specify a parameters file as the next option.
Optimize model parameters = yes

# This optional field must be used if you do not optimize the substitution
# parameters. It can also be used if you wish to start the ML optimizations
# from a specific set of parameter values
# Starting model parameters file = notusedhere.mod
{\MODEL}

# A TREE block, specifying the set of candidate tree topologies.
{TREE}
# A file containing one or more candidate topologies.
Tree file = sequence-data/hiv6.tre

# You must specify an outgroup although it is used for representational
# purposes only, and it does not affect the results.
# This outgroup must be the name of a species in your datafile or the name
# of a monophyletic clade in your clusters file (see below).
Outgroup = outgroup

# Optional: specify a file that contains monophyletic clades.
# It only affects the presentation of the phylogeny, not the max-likelihood.
Clusters file = sequence-data/hiv6.cls
{\TREE}

Random seed = 29072011

Output file = results/optimizer.hiv6-TN93I


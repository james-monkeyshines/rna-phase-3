# Phylogenetic tree reconstruction in the ML framework with mlphase.
# The dataset in this example is small and mlphase can be used.

# A standard DATAFILE block (see manual).
{DATAFILE}
Data file = sequence-data/hiv6.dna
Interleaved data file = yes
{\DATAFILE}

# A standard MODEL block (see manual), with an option to optimize or
# use user-defined substitution parameters.
{MODEL}
Model = HKY85
Discrete gamma distribution of rates = no
Invariant sites                      = yes

# If this value is 'no', you must specify a parameters file as the next option.
Optimize model parameters = yes

# This optional field must be used if you do not optimize the substitution
# parameters. It can also be used if you wish to start the ML optimizations
# from a specific set of parameter values
# Starting model parameters file = notusedhere.mod
{\MODEL}

# A TREE block
{TREE}
# You must specify an outgroup although it is used for representational
# purposes only, and it does not affect the results.
# This outgroup must be the name of a species in your datafile or the name
# of a monophyletic clade in your clusters file (see below).
Outgroup = outgroup

# See manual for the available heuristic/exhaustive search methods.
Search algorithm = Star decomposition

# Optional: specify a file that contains monophyletic clades.
# Tree topologies that do not match these constraints are not evaluated.
Clusters file = sequence-data/hiv6.cls
{\TREE}

Random seed = 29072011

Output file = results/mlphase.hiv6-HKY85I.out


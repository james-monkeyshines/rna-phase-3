# A simple example with mcmcphase.

# A standard DATAFILE block for analysing DNA sequences
{DATAFILE}
Data file = sequence-data/whales31-cytb.dna
Interleaved data file = no
# Ignore the class line at the end of this sequence file and
# group the nucleotides in a single class.
Heterogeneous data models = no
{\DATAFILE}


# A standard 4-state DNA model, across-site rate heterogeneity is accounted
# for with the discrete gamma model (6 rate categories)
{MODEL}
    Model = TN93
    Discrete gamma distribution of rates = yes
    Number of gamma categories  = 6
    Invariant sites             = no
{\MODEL}

# The standard TREE block for mcmcphase
{TREE}
Tree = Unrooted MCMC tree
Outgroup = Chevrotain
{\TREE}

# Tuning parameters
{PERTURBATION}
Tree, proposal priority = 20
Model, proposal priority = 1
{PERTURBATION_TREE}
Topology changes, proposal priority = 1
Branch lengths, proposal priority = 10
# One has to specify the prior on branch lenghts, a uniform(0,10) could have been used
# but see the scientific Literature for possible issues.
Branch lengths, prior = exponential(10)
{\PERTURBATION_TREE}
{PERTURBATION_MODEL}
    Frequencies, proposal priority = 1
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL}

{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 100000
Sampling iterations = 100000
Sampling period = 100

Output file   = results/mcmcphase.whales31-TN93dG6.mcmc
Output format = phylip


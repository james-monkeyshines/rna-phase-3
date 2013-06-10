# Bayesian inference from standard DNA sequence with MCMC.
# A simple example to use with mcmcphase/mcmcsummarize

# A standard DATAFILE block when analysing DNA sequences
{DATAFILE}
Data file = sequence-data/hiv6.dna
Interleaved data file = yes
{\DATAFILE}

# A standard MODEL block when analysing DNA sequences.
# We define a HKY85 model here.
# An invariant category is used to model across-sites rate variation.
{MODEL}
Model = HKY85
Discrete gamma distribution of rates = no
Invariant sites = yes
# Substitution parameters are initialized from the sequences by default.
# Alternatively, you can specify the initial parameters with a "model file".
# This is important if you wish the parameters to remain constant during
# the MCMC run.
# Starting model parameters file = not_used_here
{\MODEL}

{TREE}
# Here, we assume the molecular clock and use an ultrametric tree with a root.
# Nevertheless, you are strongly encouraged to use the standard "Unrooted MCMC tree"
# in most cases.
Tree = Rooted MCMC tree with molecular clock

# The "outgroup" is optional with rooted trees. When specified, it is used to
# fix the position of the root. This field would be compulsory if an unrooted
# tree was used but it has no effect on the results in this case.
# Here a clusters file containing monophyletic clades is specified and the
# outgroup is set to be one of these clades.
Outgroup = outgroup
# Monophyletic clades can be specified to enforce some restrictions on the space
# of tree topology.
Clusters file = sequence-data/hiv6.cls
# By default, the MCMC run starts from a random tree but this can be changed
# This is used when trying to estimate Bayesian posterior probability
# distributions assuming a fixed tree (see below).
# Starting tree file = not used here
{\TREE}

# Tuning parameters for the MCMC run (see manual for more info).
{PERTURBATION}
# Relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 2
Model, proposal priority = 1

{PERTURBATION_TREE}
Topology changes, proposal priority = 1
Branch lengths, proposal priority = 5
# A uniform prior could be used here:
Tree height, prior = exponential(1)
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
    Frequencies, proposal priority = 3
    # example: a prior different than the default dirichlet(1,1,1,1) for frequencies
    Frequencies, prior = dirichlet(2,2,2,2)
    # example: a prior different than the default uniform prior for the "rate ratios", a hyperparameter is introduced
    Rate ratios, prior = exponential(exponential(1))
    Rate ratios, proposal priority = 1
    Rate ratios exponential hyperparameter, proposal priority = 1
    Invariant parameter, proposal priority = 1
{\PERTURBATION_MODEL}
{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 50000
Sampling iterations = 100000
Sampling period = 100

Output file   = results-check/mcmcphase.hiv6-HKY85I.mcmc
Output format = phylip


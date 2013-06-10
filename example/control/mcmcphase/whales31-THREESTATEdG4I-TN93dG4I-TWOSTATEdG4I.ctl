# An example to analyse protein-coding genes with different DNA models
# for the different codon positions with mcmcphase.
# This is similar to Gibson et al. but with a different dataset (and
# keeping the third codon position).

# We have to set up a MIXED model but we cannot rely on the
# default behaviour loops=MODEL1 and stems=MODEL2 here
{DATAFILE}
Data file = sequence-data/whales31-cytb.dna
Interleaved data file = no
# Read the "model line" from the data file. Each nucleotide is associated
# to one of the three classes
Heterogeneous data models = yes
{\DATAFILE}

# Use a MIXED model and specify the substitution model to use with each class.
{MODEL}
Model = MIXED
Number of models = 3
{MODEL1}
    Model = THREESTATE
    Discrete gamma distribution of rates = yes
    Number of gamma categories  = 4
    Invariant sites             = yes
{\MODEL1}
{MODEL2}
    Model = TN93
    Discrete gamma distribution of rates = yes
    Number of gamma categories  = 4
    Invariant sites             = yes
{\MODEL2}
{MODEL3}
    Model = TWOSTATE
    Discrete gamma distribution of rates = yes
    Number of gamma categories  = 4
    Invariant sites             = yes
{\MODEL3}
{\MODEL}

{TREE}
Tree = Unrooted MCMC tree
Outgroup = Chevrotain
{\TREE}

# # # # # # # # # # # # # # # # # # # # # # # #  The MCMC PERTURBATION Section # # # # # # # # # # # # 

{PERTURBATION}

# Tree/model relative proposal probabilities
Tree, proposal priority = 10
Model, proposal priority = 1

{PERTURBATION_TREE}
Topology changes, proposal priority = 1
Branch lengths, proposal priority = 8
Branch lengths, prior = exponential(10)
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
# We tune (approximately) the probabilities for each substitution model by comparing the
# number of free parameters in each
Model 1, proposal priority = 6
Model 2, proposal priority = 8
Model 3, proposal priority = 3
# Proposal probability for the average substitution rates of MODEL2 and MODEL3
# (MODEL1 is constrained equal to 1.0)
Average rates, proposal priority = 1
{PERTURBATION_MODEL1}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
    Invariant parameter, proposal priority = 1
{\PERTURBATION_MODEL1}
{PERTURBATION_MODEL2}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
    Invariant parameter, proposal priority = 1
{\PERTURBATION_MODEL2}
{PERTURBATION_MODEL3}
    Frequencies, proposal priority = 2
    Gamma parameter, proposal priority = 1
    Invariant parameter, proposal priority = 1
{\PERTURBATION_MODEL3}
{\PERTURBATION_MODEL}

{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 300000
Sampling iterations = 600000
Sampling period = 200

Output file   = results/mcmcphase.whales31-THREESTATEdG4I-TN93dG4I-TWOSTATEdG4I.mcmc
Output format = phylip


# Using different models for the LSU & SSU genes, an example with mcmcphase

{DATAFILE}
Data file = sequence-data/tol40.rna
Interleaved data file = yes
Interleaved structure = no
# DO NOT use the "automatic method" to analyse this dataset:
# we need the nucleotides to be distributed according to the
# "model line" in the data file
Heterogeneous data models = yes
{\DATAFILE}


# Set up a MIXED model with two different TN93 substitution models for the
# two blocks with unpaired nucleotides (one for each gene) and two 7D for
# the pairs (one for each gene)
{MODEL}
Model = MIXED
Number of models = 4
 {MODEL1} # LSU pairs
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
 {\MODEL1}
 {MODEL2} # LSU unpaired nuc
  Model = TN93
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
 {\MODEL2}
 {MODEL3} # SSU pairs
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
 {\MODEL3}
 {MODEL4} # SSU unpaired nuc
  Model = TN93
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites             = no
 {\MODEL4}

{\MODEL}

# Use a standard unrooted tree.
# The outgroup is compulsory but does not affect the results.
{TREE}
Tree = Unrooted MCMC tree
Outgroup = Thermotoga_maritima
{\TREE}


# Tuning parameters for the MCMC runs.
{PERTURBATION}

# Relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 5
Model, proposal priority = 1

{PERTURBATION_TREE}
# We use 10/40 for topology change vs branch length changes.
# It is not exactly equivalent to 1/4 because this is also given relative to the
# proposal priority for hyperparameters that are introduced with the
# the prior on branch lengths (Hyperpriors, proposal priority)
Topology changes, proposal priority = 10
Branch lengths, proposal priority = 40
Hyperpriors, proposal priority = 1

# We use a vague prior exp(lambda) on branch lengths rather than the default exp(10)
Branch lengths, prior = exponential(uniform(0,100))
# A lambda hyperparameter has been introduced. It needs a "proposal priority"
# but this is not used because it is the only hyperparameter
Branch lengths exponential hyperparameter, proposal priority = 1
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
# Relative probabilities for the proposals on the two models and the average substitution rate of MODEL2 
Model 1, proposal priority = 7
Model 2, proposal priority = 8
Model 3, proposal priority = 7
Model 4, proposal priority = 8
# Relative probabilities for the proposals that change the average substitution rate of MODEL2, MODEL3, MODEL4
Average rates, proposal priority = 1
{PERTURBATION_MODEL1}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL1}
{PERTURBATION_MODEL2}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL2}
{PERTURBATION_MODEL3}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL3}
{PERTURBATION_MODEL4}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL4}
{\PERTURBATION_MODEL}

{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 750000
Sampling iterations = 1500000
Sampling period = 150


Output file   = results/mcmcphase.tol40.mcmc
Output format = phylip


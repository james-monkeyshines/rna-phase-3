# An example similar to Hudelot et al. (2003)

# A standard DATAFILE block for RNA sequences having a secondary structure.
{DATAFILE}
Data file = sequence-data/mammals69.rna
Interleaved data file = no
# Use the "automatic method" to analyse this dataset:
# unpaired nucleotides ('.' in the secondary structure) are
# handled by the MODEL1 of the MIXED model (see below).
# pairs (corresponding parenthesis in the secondary structure)
# are handled by the MODEL2 of the MIXED model (see below)
Heterogeneous data models = auto
{\DATAFILE}


#Set up a MIXED model with REV for loops and 7D for stems
{MODEL}
Model = MIXED
Number of models = 2
  {MODEL1}
  Model = REV
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
  {\MODEL1}
  {MODEL2}
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
{\MODEL2}
{\MODEL}

# Use a standard unrooted tree.
# The outgroup is compulsory but does not affect the results.
{TREE}
Tree = Unrooted MCMC tree
Outgroup = ORNANAMIT
{\TREE}


# Tuning parameters for the MCMC runs.
{PERTURBATION}

# Relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 8
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
Model 1, proposal priority = 10
Model 2, proposal priority = 10
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
{\PERTURBATION_MODEL}

{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 750000
Sampling iterations = 1500000
Sampling period = 150


Output file   = results/mcmcphase.mammals69-REVdG6-7DdG6.mcmc
Output format = phylip


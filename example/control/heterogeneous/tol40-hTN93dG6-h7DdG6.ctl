# An example of mcmcphase, with different models for the LSU & SSU genes.

{DATAFILE}
Data file = sequence-data/tol40.rna
Interleaved data file = yes
Interleaved structure = no

# We use the "automatic method" and we ignore the "model line"
# in the data file. LSU and SSU are considered as a single
# gene (although loops and stems are partitioned).
Heterogeneous data models = auto
{\DATAFILE}


# Set up a "standard" MIXED model for RNA genes and make it HETEROGENEOUS
{MODEL}
 Model = HETEROGENEOUS
 # We use 3 different MIXED models, all their parameters are independent
 Number of models = 3
  {BASEMODEL}
   Model = MIXED
   Number of models = 2
    {MODEL1}
     Model = TN93
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
  {\BASEMODEL}
{\MODEL}

# Use the rooted heterogeneous tree for time-heterogeneous model.
# The outgroup is optional but it is recommended to use one and to force
# the position of the root (the method does not recover very well the
# position of the root and it does not mix properly in our experience).
{TREE}
Tree = Heterogeneous MCMC tree
# We use a clade file to define 2 monophyletic clades (eubacteria
# and eukaryotes). These 2 clades are "protected" during the run.
Clusters file = sequence-data/tol40.cls
# We can now position the root on the branch leading to eubacteria
Outgroup = eubacteria
{\TREE}


# Tuning parameters for the MCMC runs
{PERTURBATION}

# Relative proposals probabilities between the tree and the substitution model
Tree, proposal priority = 3
Model, proposal priority = 1

{PERTURBATION_TREE}
Topology changes, proposal priority = 1
Branch lengths, proposal priority = 7
Branch lengths, prior = exponential(10)
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
Model, proposal priority = 20
# Ancestral frequencies for the MODEL1 (loops)
{ANCESTRAL_FREQUENCIES1}
  Frequencies, proposal priority = 1
{\ANCESTRAL_FREQUENCIES1}
# Ancestral frequencies for the MODEL2 (base-pair)
{ANCESTRAL_FREQUENCIES2}
  Frequencies, proposal priority = 1
{\ANCESTRAL_FREQUENCIES2}
# IMPORTANT: keep the next field to 0
Average rates, proposal priority = 0
  # Now the "standard" perturbation block for a MIXED model
  {PERTURBATION_BASEMODEL}
  Model 1, proposal priority = 7
  Model 2, proposal priority = 8
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
  {\PERTURBATION_BASEMODEL}
{\PERTURBATION_MODEL}
{\PERTURBATION}

Random seed = 29072011

Burnin iterations = 2000000
Sampling iterations = 10000000
Sampling period = 1000


Output file   = results/heterogeneous.tol40-hTN93dG6-h7DdG6.mcmc
Output format = phylip


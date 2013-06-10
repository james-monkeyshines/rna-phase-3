# An example modelling compositional heterogeneity across-sites
# when composition changes smoothly with the rate of evolution

# A simulation similar to Gowri-Shankar and Rattray (2006)

# The standard DATAFILE block when RNA sequences with a secondary structure are used.
{DATAFILE}
Data file = sequence-data/primates16.rna
Interleaved data file = no
Heterogeneous data models = auto
{\DATAFILE}

# We use a TN93 model for loops and a 7D for stems. Both substitution models
# account for rate heterogeneity across sites with a discrete gamma distribution
# (8 categories). Different compositional parameters are used in each rate
# category with the base-pair model.
{MODEL}
Model = MIXED
Number of models = 2
 {MODEL1}
  Model = TN93
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 8
  Invariant sites            = no
 {\MODEL1}
 {MODEL2}
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 8
  Number of frequencies sets = 8
  Invariant sites            = no
 {\MODEL2}
{\MODEL}

# The standard TREE block with mcmcphase
{TREE}
Tree = Unrooted MCMC tree
Outgroup = BOSTAUMIT
{\TREE}


######################## The MCMC PERTURBATION Section ############
# Tune the proposal probabilities
{PERTURBATION}
Tree, proposal priority = 1
Model, proposal priority = 6

{PERTURBATION_TREE}
Topology changes, proposal priority = 2
Branch lengths, proposal priority = 15
Branch lengths, prior = exponential(10)
{\PERTURBATION_TREE}

{PERTURBATION_MODEL}
Model 1, proposal priority = 4
Model 2, proposal priority = 30
Average rates, proposal priority = 1
{PERTURBATION_MODEL1}
    Frequencies, proposal priority = 2
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL1}
{PERTURBATION_MODEL2}
    Frequencies, proposal priority = 30
    Frequencies, initial Dirichlet tuning parameter = 70000
    Frequencies variation prior = gp(exponential(20),exponential(1),uniform(0,10),uniform(0,10),0.005,0)
    Frequencies variation prior theta1 hyperparameter, proposal priority = 3
    Frequencies variation prior scale hyperparameter, proposal priority = 3
    Frequencies variation prior theta2 hyperparameter, proposal priority = 3
    Frequencies variation prior linear trend hyperparameter, proposal priority = 3
    Rate ratios, proposal priority = 1
    Gamma parameter, proposal priority = 1
{\PERTURBATION_MODEL2}
{\PERTURBATION_MODEL}
{\PERTURBATION}

Random seed = 29072011

# This is just for illustation purpose (approx 1 to 2 days of computation).
# Multiply burn-in and sampling by 10 and check for convergence if reliable
# results are required.
Burnin iterations = 1000000
Sampling iterations = 4000000
Sampling period = 1000

Output file   = results/heterogeneous.primates16-TN93dG8-7DdF8.mcmc
Output format = phylip


# An example with the program likelihood.
# likelihood computes the likelihood of a phylogeny when the
# substitution parameters are known

# This program is useful because:
# 1)it can perform ancestral sequence reconstruction
# 2)it can compute the Bayesian Posterior Probabilities
# for the different categories of a mixture model (it can be used
# to compute site-specific substitution rate, Yang, 95)
# 3)it can also output the site-specific loglikelihoods used in
# some statistical tests.
# 4)it can also be used for the method proposed in Gowri-Shankar
# et al.(2006) to approximate the equilibrium frequencies in
# different rate categories

# A standard DATAFILE block for RNA sequences
{DATAFILE}
Data file = sequence-data/s70-TN93dG6-7DdG6.rna
Interleaved data file = no
Heterogeneous data models = auto
{\DATAFILE}

# A standard MODEL block for RNA sequences
{MODEL}
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

# User-specified model parameters are compulsory
# (see the manual/examples for more info on "model files").
Model parameters file = results-check/optimizer.s70-TN93dG6-7DdG6.mod
{\MODEL}

{TREE}
Tree file = results-check/optimizer.s70-TN93dG6-7DdG6.tre
{\TREE}

# The five following fields are optional

# Filename to output results for the ancestral sequences marginal reconstruction.
# BEWARE1: 7-state models reconstruct MM pairs, R/Y models reconstruct 0/1 sequences.
# BEWARE2: When used with base-pair models (for instance), the program REPEATS the
#          BPP for the pair at each nucleotide position.
Ancestral sequences = results/likelihood.s70-TN93dG6-7DdG6.anc

# Filename to output the results of the site-specific "BPP" rate categories estimation.
# BEWARE: When used with base pair models (for instance), the program REPEATS the MAP
#         category estimate and the BPP at each nucleotide position
Site-specific substitution rates = results/likelihood.s70-TN93dG6-7DdG6.rat

# Filename to output the site-specific loglikelihood.
# BEWARE1: These site-specific likelihoods are not given in the order you might expect.
#          Invariant sites are at the end for technical reasons. Moreover, sites are
#          classed according to their data type (loops/stems,..) with MIXED models.
# BEWARE2: With base-pair models (for instance), the site-specific loglik appears
#          ONLY ONCE (it would probably not serve any purpose to repeat it and it
#          could be error-prone).
Site-specific loglikelihood = results/likelihood.s70-TN93dG6-7DdG6.ssl

# Compute rate-specific composition.
# Composition for a given rate category is computed from the frequencies
# observed at each site weighted by the BPP of the site-specific rate category.
Rate vs. composition = yes

# You can limit the rate-specific composition estimation to a set of species:
# For a single species simply use its name;
# For all the species use "all" or simply omit the field (default behaviour);
# For a specific subset specify the name of a file (and use the name of one
# species per line in this file).
Selected species = all


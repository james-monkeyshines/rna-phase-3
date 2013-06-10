# A control file to be used with simulate, with some advanced settings:
#  Generate a RNA sequence with secondary structure information.
#  Generate a random tree.
# Please see the other example control file for the basic settings.

# Set up a MIXED model
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
{\MODEL}

# Change the following field to yes to generate the skeleton file
# for this MIXED model.
Retrieve the name of the model's parameters = no

Model parameters file = sequence-data/s70-TN93dG6-7DdG6.mod

Random seed = 29072011

# We generate a random tree (but the number of species is specified by the user)
# see Aldous, 1996 and Yang, 1997 for information on these topics.
{TREE}
# The generated phylogeny (tree topology + branch lengths)
# is stored in the following file
Generate tree = yes
Tree file = results/simulate.s70.tre

# Total number of species in the generated tree
# (including the outgroup if there is one).
Number of species = 70
# Optional parameter to create an outgroup cluster.
Number of species in the outgroup = 7

# The process used to generate the tree; options are:
# Yule process, Birth-death process, Uniform process, Beta-splitting process.
Generating model = Birth-death process

# The Beta-splitting process requires the specification of beta in [-2;+inf]
# Beta parameter=-2 equivalent to comb
# Beta parameter=-1.5 equivalent to the uniform process
# Beta parameter=0 equivalent to a Yule or birth-death process
# Beta parameter=+inf or beta=inf or beta=+infinity is a symmetric binary trie

# Yule process and Birth-death process incorporate automatically
# a probability distribution for branch lengths which cannot be changed.
# When using the uniform and/or beta-splitting process one has to specify
# a prior distribution for the branch lengths.
# Choose among Uniform, Exponential, Pure-birth process, Birth-death process
# If using Uniform then you also have to specify two extra parameters:
# "Branch prior, lower bound" and "Branch prior, upper bound"
# If using Exponential then you have to specify "Branch prior, exponential parameter" (=1/mean)

# When using the Yule or the Birth-death process as a generating model or when
# using the Pure-birth process or the Birth-death process for the prior on branch
# lengths, you have to specify a "Birth rate", and a "Death rate" if appropriate
# NB: distance from root to tips is supposed to be 1.0 for this procedure
# See Yang, Rannala 1997 for more information.
Birth rate = 10
Death rate = 5

# When using a Yule or birth-death process, the height of the tree is
# rescaled from 1.0 to the user's choice.
Tree height = .8

# When using the Yule or the Birth-Death process, you can specify a species
# sampling (see Yang,Rannala (1997)) (it is 1.0 by default).
# Accounting for the fact that your dataset does not contain all existing
# species affects the prior on branch lenghts.
Species sampling = .05
# You can choose to use a different value for the outgroup (if you previously
# used the field "Number of species in the outgroup").
Outgroup sampling = .001
{\TREE}


# The number of symbols to generate (a symbol can be a pair or a complete codon
# depending on your model). Since a MIXED model is used here, you have to specify
# a length for each model/class/type of data.
Number of symbols from class 1 = 1500
# For the stems, there are 2000 nucleotides but only 1000 symbols.
Number of symbols from class 2 = 1000

# The output of the program
# IMPORTANT 1: PHASE is generating the nucleotides for the first class, then the nucleotides for the second class
#              before concatenating the results. Heterogeneous data types are not intertwined.
#              In this case it means that the first 1500 nucleotides are not paired, the 2000 remaining nucleotides
#              are generated with the RNA model.
# IMPORTANT 2: Note that the two nucleotides of a pair are at position i/i+1 in the generated sequences
# IMPORTANT 3: 7-state RNA models are using AA for the mismatch state, the THREESTATE model is using C for the Y state
#              the TWOSTATE model is using A and C for its two states.
Output file = results/simulate.s70-TN93dG6-7DdG6.rna


# Optional parameters follow. They are used to produce a sequence file fully
# compatible with other programs in the PHASE package and they can be left out.

# The data file type to use at the first line of the sequence file
Data file type = STRUCT

# To produce a correct structure line in this case, PHASE needs the following:
Structure for the elements of class 1 = .
Structure for the elements of class 2 = ()
# with a codon model we would have used "123" (without quotes)

# This is the number of nucleotides for each data type which is used to
# compute the total length of the alignment properly.
Number of nucleotides from class 1 = 1500
Number of nucleotides from class 2 = 2000

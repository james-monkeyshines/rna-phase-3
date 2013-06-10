# Use this control file with simulate to evolve DNA sequences along
# the branches of a tree.

# The user has to specify the substitution model used.
# Here we use the TN93 model and rate heterogeneity across sites
# is modelled with the discrete gamma model (10 categories).
{MODEL}
    Model = TN93
    Discrete gamma distribution of rates = yes
    Number of gamma categories = 10
    Invariant sites = no
{\MODEL}

# The user also has to specify the parameters for this substitution model.
# Parameters are contained in a "model file"
Model parameters file = sequence-data/s10-TN93dG10.mod

# The format of this "model file" is not straightforward and depends on the
# substitution model. It is recommended to let PHASE generate a skeleton file
# for you first. Change the following field to "yes" to do so, simulate will
# not generate any sequences and will simply create/overwrite the model file.
# Replace the default parameters with your own values and do not forget to
# change the following field back to "no" 
Retrieve the name of the model's parameters = no

# Change the random seed to generate different set of sequences
Random seed = 29072011

# Simulate on a fixed tree; see the manual (or the other example control
# file) for details on random tree generation.
{TREE}
Tree file = sequence-data/s10.tre
Generate tree = no
{\TREE}

# In this simple case we have to specify the lenght of the sequences, ie, the
# number of nucleotides, but things can become more complicated with complex
# substitution models.
# When using a base-pair or a codon model, you do not specify a number of nucleotides
# but a number of pairs or a triplets.
# When a MIXED model is used (ie more that one type of data), you have to specify a
# number of symbols for each substitution model.
Number of symbols from class 1 = 4000

# The file where the sequences are written
Output file = results/simulate.s10-TN93dG10.dna

# Optional parameters follow. They are used to produce a sequence file fully
# compatible with PHASE (ie, that can be used directly with other programs
# of the package) but they are not compulsory. Manual editing of these files
# might be easier so you should not bother with these unless it is necessary.

# The third token to use in the first line of your sequence file (can also be
# CODON or STRUCT; see manual). For simple nucleotide sequences generated
# with a DNA substitution model use 'DNA'.
Data file type = DNA

# If you use STRUCT in the previous case, you have to tell simulate how the
# structure line should be produced. For instance:
# Structure for the elements of class 1 = .
# There should be more than one field when a MIXED model is used.
# If class x is a nucleotide or amino-acid model use: Structure for the elements of class x = "."
# If class x is a doublet model use: Structure for the elements of class x =  "()"
# If class x is with a codon model use: Structure for the elements of class x = "123"

# The last optional parameter is the total number of nucleotides in the
# alignment (the second field in the first line).
Number of nucleotides from class 1 = 4000


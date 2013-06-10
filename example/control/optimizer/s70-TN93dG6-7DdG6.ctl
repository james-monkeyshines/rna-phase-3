# An example using a MIXED model with optimizer.
# Loops and stems of (simulated) RNA sequences are
# handled with different substitution processes.

# A standard DATAFILE block for a MIXED model with RNA loops and stems.
{DATAFILE}
Data file = sequence-data/s70-TN93dG6-7DdG6.rna
Interleaved data file = no
# Use "auto" for the following field when analysing standard RNA sequences.
Heterogeneous data models = auto
{\DATAFILE}

# A standard MODEL block for a MIXED model:
# MODEL1 is used for unpaired nucleotides;
# MODEL2 is used with pairs.
# (See the manual and other examples for more complex cases.)
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

# Use fixed model parameters (ie only optimize branch lengths).
Optimize model parameters = no
Starting model parameters file = sequence-data/s70-TN93dG6-7DdG6.mod
{\MODEL}

{TREE}
Tree file = sequence-data/s70.tre
# An outgroup is mandatory, but is only used for displaying the tree.
# It can be a single species, but here we use a monophyletic clade
# defined in the clusters file.
Outgroup = outgroup
Clusters file = sequence-data/s70.cls
{\TREE}

Random seed = 29072011

Output file = results/optimizer.s70-TN93dG6-7DdG6


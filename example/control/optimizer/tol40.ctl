# Using different models for the LSU & SSU genes, an example with optmiizer.

{DATAFILE}
Data file = sequence-data/tol40.rna
Interleaved data file = yes
# The structure and the model line are not interleaved.
Interleaved structure = no
# DO NOT use the "automatic method" to analyse this dataset;
# we need the nucleotides to be distributed according to the
# "model line" in the data file.
Heterogeneous data models = yes
{\DATAFILE}

# A MIXED model with two different TN93 substitution models for the
# two blocks with unpaired nucleotides (one for each gene) and
# two 7D models for the pairs (one for each gene).
{MODEL}
Model = MIXED
Number of models = 4
 {MODEL1} #LSU pairs
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
 {\MODEL1}
 {MODEL2} #LSU unpaired nuc
  Model = TN93
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
 {\MODEL2}
 {MODEL3} #SSU pairs
  Model = RNA7D
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
 {\MODEL3}
 {MODEL4} #SSU unpaired nuc
  Model = TN93
  Discrete gamma distribution of rates = yes
  Number of gamma categories = 6
  Invariant sites = no
 {\MODEL4}

{\MODEL}

# A TREE block, specifying a topology whose branch lengths will be optimized.
{TREE}
Tree file = sequence-data/tol40.tre

# You must specify an outgroup although it is used for representational
# purposes only, and it does not affect the results.
Outgroup = Thermotoga_maritima
{\TREE}

Random seed = 29072011

Output file   = results-check/optimizer.tol40


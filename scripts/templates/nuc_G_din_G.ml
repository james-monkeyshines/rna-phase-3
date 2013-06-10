# Template of a control file for use with mlphase or optimizer.

# The name of the data file needs to be added programmatically.
# For simplicity we assume that the data isn't interleaved.
{DATAFILE}
Data file = 
Interleaved data file = no
Heterogeneous data models = auto
{\DATAFILE}

{MODEL}
Model = MIXED
Number of models = 2
Optimize model parameters = yes
Empirical frequencies = yes

{MODEL1}
Model = REV
Discrete gamma distribution of rates = yes
Number of gamma categories = 4
Invariant sites = no
{\MODEL1}

{MODEL2}
Model = RNA7A
Discrete gamma distribution of rates = yes
Number of gamma categories = 4
Invariant sites = no
{\MODEL2}

{\MODEL}

{TREE}
# mlphase needs to know how to search tree space
Search algorithm = Star decomposition
#Search algorithm = Exhaustive

# optimizer needs a tree file, to be added programmatically.
Tree file = 

# The outgroup is mandatory, and must match a species in the datafile;
# the name of this species needs to be added programmatically.
Outgroup = 
{\TREE}

# Add a random seed programmatically, rather than using the system time (the default).
Random seed = 

# The name of the output file needs to be added programmatically.
Output file = 


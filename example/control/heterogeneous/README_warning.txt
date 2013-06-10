Heterogeneity across sites and over time
----------------------------------------
The methods can be tried but this is still work in progress.
Monitor for convergence carefully and check the outputs
of the software. You are welcome to ask for help.

We are still working on the Gaussian Process method proposed
in Gowri-Shankar et al. (2006) to account for compositional
heterogeneity across sites. We are currently revising the
MCMC proposal mechanisms to improve the mixing. The method
proposed in this version of PHASE cannot handle more than a
dozen species in a reasonable amount of time. We also have
to test the behaviour of this GP prior on a broad range of
datasets.

Current work on time-heterogeneous methods is more advanced.
We completed a reversible-jump MCMC algorithm that allows for
a variable number of substitution models on the branches of
the tree. The number of models is determined by the amount
of heterogeneity evidenced by the data.
This version of PHASE requires the user to specify
a fixed number of models and it is not possible to choose
which substitution parameters are heterogeneous/constant
over time. Nevertheless, you can ask for the
development version we are currently working with.


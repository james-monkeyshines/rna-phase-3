# The following example will create subdirectories called 'control' and
# 'results' in the 'example' directory, to store intermediate files and
# the output files from PHASE. The model selection results will be saved
# in a file example/RF00403_results.txt.

# The MCMC tree search can be time-consuming, so that option has been
# commented out below; simply uncomment if you want to use that option.

perl model_selection.pl \
  --alignment example/RF00403.fa \
  --structure example/RF00403.structure.txt \
  --tree_file example/RF00403.bionj.nh \
  --out_dir example \
#  --mcmc_tree_search


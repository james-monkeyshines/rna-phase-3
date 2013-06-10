#!/usr/bin/perl
use strict;
use warnings;

use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Getopt::Long;

use FindBin;
use lib $FindBin::Bin;

use PhaseUtils
	qw (read_from_file save_to_file
		from_fasta to_phylip_s
		parse_structure sample_size pair_gaps
		phase_seq_file phase_ctl_file
		phase_ml phase_mcmc phase_parse phase_mcmc_tree
		aic aicc delta_aic
);

=pod

=head1 Name

model_selection.pl

=head1 Description

Calculate likelihoods and AIC values for a range of RNA models.
(Execute 'model_selection.pl --help' for further details.)

=head1 Author

James Allen
(james@monkeyshines.co.uk)

=head1 Copyright

Copyright 2012 James Allen; University of Manchester

=head1 Licence

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

=cut

################################################################################
# Arguments - Start
my(	$help, $verbose, $overwrite, $tree_search, $exhaustive_tree_search,
	$alignment_file, $structure_file, $tree_file,
	$out_dir, $phase_optimizer, $phase_mcmc, $template_dir, $sample_type,
);

GetOptions(
	'help' => \$help,
	'verbose' => \$verbose,
	'replace' => \$overwrite,
	'mcmc_tree_search' => \$tree_search,
	'exhaustive_mcmc' => \$exhaustive_tree_search,
	'alignment_file=s' => \$alignment_file,
	'structure_file=s' => \$structure_file,
	'tree_file=s' => \$tree_file,
	'out_dir:s' => \$out_dir,
	'phase_optimizer:s' => \$phase_optimizer,
	'phase_mcmc:s' => \$phase_mcmc,
	'phase_template_dir:s' => \$template_dir,
	'calculate_sample_size:s' => \$sample_type,
);

sub usage {
  print qq~

model_selection.pl: Model Selection with RNA Models

This script performs maximum likelihood inference of model parameters and
branch lengths for a fixed topology, using a range of DNA and RNA models.
AIC values are calculated in order to choose the model which best fits the
data. The script can then perform MCMC tree search using this best model.

* Mandatory parameters:
  --a[lignment_file] <file>
  --s[tructure_file] <file>
  --t[ree_file] <file>
  Three input files are required, an alignment in Fasta format, a structure
  in dot-bracket notation, and a tree in Newick format. Readseq is included
  in this distribution of PHASE (../readseq), so please use that to convert
  to Fasta format, if necessary.

* Optional parameters:
  --o[ut_dir] <dir>
  If an output directory is not specified, results are saved in the current
  directory. The PHASE input and output files are saved in subdirectories,
  'control' and 'results', respectively. The output is summarised in a text
  file in the output directory.

  --phase_o[ptimizer] <file>
  --phase_m[cmc] <file>
  This script will first look for the required PHASE executable files in the
  directory '../bin'. If not found, the executables are assumed to be in your
  path. If they are in neither of these locations, paths must be provided as
  parameters.

  --phase_t[emplate_dir] <dir>
  The script uses templates to generate the control files that are required
  by PHASE; these are stored in the './templates' directory. If this
  directory cannot be found, a paths must be provided as a parameter.

  --c[alculate_sample_size] ('characters'|'sites')
  Calculation of corrected AIC (AICc) values requires a sample size.
  It is not clear what this means in the context of sequence alignments,
  but two commonly used approximations are the number of characters
  and the number of sites.

* Optional flags:
  --v[erbose]
  Display informational messages about what the script is doing.

  --r[eplace]
  Existing results will not be overwritten unless this option is used.

  --m[cmc_tree_search]
  After model selection, do MCMC tree search with the best model. The
  ML estimates are used as a starting point for the model parameters.
  The MCMC search has 150,000 burn-in iterations, 300,000 sampling
  iterations, and a sampling period of 100. See the *.mcmc files in the
  './templates' directory for further details. The majority-rule consensus
  tree (nodes labelled with posterior probabilities) is saved in the
  output directory.
  For comparison, the MCMC tree search is repeated for the simplest
  nucleotide+G model (HKY+G, unless you edit the script).

  --e[xhaustive_mcmc]
  After model selection, do MCMC tree search with all models (parameters
  the same as the -m option). Note that this could be computationally
  intensive, so is not generally recommended.

~;
  print 'Author: James Allen (james@monkeyshines.co.uk)'."\n\n";
  print "Usage: $0 --align <file> -str <file> -tree <file> [options]\n";
  print qq~
  Please provide values for:
    --a[lignment_file]  Fasta-format alignment
    --s[tructure_file]  Secondary structure in dot-bracket notation
    --t[ree_file]       Newick-format tree
  
  Optional parameters:
    --o[ut_dir]                default '.'
    --phase_o[ptimizer]        default '../bin/optimizer'
    --phase_m[cmc]             default '../bin/mcmcphase'
    --phase_t[emplate_dir]     default './templates'
    --c[alculate_sample_size]  default 'characters'

  Optional flags:
    --h[elp]              View this help page
    --v[erbose]           Display progress on screen
    --r[eplace]           Replace any existing results files without warning
    --m[cmc_tree_search]  After model selection, do tree search with best model
    --e[xhaustive_mcmc]   After model selection, do tree search with all models
  
  ~;
  print "\n";
  exit;
}

usage() if (defined $help);

unless ( $alignment_file && $structure_file && $tree_file ) {
	usage;
}

$sample_type = 'characters' unless $sample_type;
$verbose = 0 unless $verbose;
$overwrite = 0 unless $overwrite;
$tree_search = 0 unless $tree_search;
$exhaustive_tree_search = 0 unless $exhaustive_tree_search;
# Arguments - End
################################################################################

################################################################################
# Model Set-Up - Start
# The models to be executed have one of four possible patterns:
#   1) a single nucleotide+G model for stems and loops
#   2) different nucleotide+G models for stems and loops
#   3) a nucleotide+G model for loops and a dinucleotide model for stems
#   4) a nucleotide+G model for loops and a dinucleotide+G model for stems
# where +G indicates gamma-distributed rates-across-sites. Each pattern
# has a corresponding template file, which contains most of the parameters
# necessary to execute PHASE (the model name and filenames are inserted
# by this script).
# Two lists of nucleotide and dinucleotide models are specified, and are
# combined with the four patterns to generate RNA models for all combinations
# of nucleotide and dinucleotide model. By default, 2 nucleotide models are
# used, and 16 dinucleotide models. This represents all of the dinucleotide
# models that PHASE knows about; additional nucleotide models (JC69, K80, F81,
# and TN93) are available, but are not used by default to keep the total number
# of RNA models manageable - add these models to the @nuc_models list if
# you're especially interested in modelling loop regions.

my %models = (
	'nuc+G' => 'nuc_G.ml',
	'nuc+G_nuc+G' => 'nuc_G_nuc_G.ml',
	'nuc+G_din' => 'nuc_G_din.ml',
	'nuc+G_din+G' => 'nuc_G_din_G.ml',
);
my @nuc_models = qw (HKY85 REV);
my @dinuc_models = qw (RNA7A RNA7B RNA7C RNA7D RNA7E RNA7F RNA7G
	RNA16A RNA16B RNA16C RNA16D RNA16E RNA16F RNA16I RNA16J RNA16K);
# Model Set-Up - End
################################################################################

################################################################################
# Sort Out Files and Directories - Start
if (!-e $alignment_file) {
	error_msg("Sequence file '$alignment_file' does not exist.");
}
if (!-e $structure_file) {
	error_msg("Structure file '$structure_file' does not exist.");
}
if (!-e $tree_file) {
	error_msg("Tree file '$tree_file' does not exist.");
}

$out_dir = '.' unless $out_dir;
if (!-e $out_dir) {
	make_path($out_dir);
	print "Writing output to new directory, '$out_dir'\n" if $verbose;
} else {
	print "Writing output to existing directory, '$out_dir'\n" if $verbose;
}
my $control_dir = "$out_dir/control";
my $results_dir = "$out_dir/results";
make_path($control_dir) unless -e $control_dir;
make_path($results_dir) unless -e $results_dir;

# If values are not given for the optional parameters, first look for
# them where they should be locally; then check the user's path.
if (!defined($phase_optimizer)) {
	$phase_optimizer = dirname($0)."/../bin/optimizer";
	if (!-e $phase_optimizer) {
		$phase_optimizer = "optimizer";
	}
}
my $optimizer_out = `$phase_optimizer 2>&1`;
if (!defined($optimizer_out) || $optimizer_out !~ /usage\s+:/m) {
	error_msg("Could not find the PHASE 'optimizer' program ($phase_optimizer).\n".
			"Use the --phase_optimizer option to specify its location.\n");
}

my $phase_mcmcsummarize;
if ($tree_search || $exhaustive_tree_search) {
	if (!defined($phase_mcmc)) {
		$phase_mcmc = dirname($0)."/../bin/mcmcphase";
		if (!-e $phase_mcmc) {
			$phase_mcmc = "mcmcphase";
		}
	}
	my $mcmc_out = `$phase_mcmc 2>&1`;
	if ($mcmc_out !~ /usage\s+:/m) {
		error_msg("Could not find the PHASE 'mcmcphase' program ($phase_mcmc).\n".
				"Use the --phase_mcmc option to specify its location.\n");
	}
	
	($phase_mcmcsummarize = $phase_mcmc) =~ s/phase$/summarize/;
}

if (!defined($template_dir)) {
	$template_dir = dirname($0)."/templates";
}
if (!-e $template_dir) {
	error_msg("Could not find the template directory '$template_dir'.\n".
			"Use the --phase_template_dir option to specify its location.\n");
}
# Sort Out Files and Directories - End
################################################################################

################################################################################
# Generate Sequence Files - Start
# PHASE requires the sequence file to be in different formats for the
# different mixed models. Also, any pairs in which one nucleotide is a gap
# character are replaced with a pair of gaps; see the PHASE 3.0 changelog
# for further details.
my $seq_data = read_from_file($alignment_file);
if ($seq_data !~ /^\s*>/m) {
	error_msg("Sequence file '$alignment_file' does not seem to be Fasta-format.\n".
			"Other formats can be converted to Fasta with Readseq, which is ".
			"included with PHASE.");
}
my ($seqs, $order) = from_fasta($seq_data);
my $seq_base = basename($alignment_file);
$seq_base =~ s/\.\w+$//;
my $seq_count = scalar(keys(%$seqs));

my $structure = read_from_file($structure_file);
(my $stem_pattern = $structure) =~ s/\./1 /g;
$stem_pattern =~ s/[^1 ]/2 /g;

pair_gaps($seqs, $structure);
my $sample_size = sample_size($seqs, $structure, $sample_type);
my %sequences;
while (my ($id, $seq) = each(%$seqs)) {
	if (exists $sequences{$seq}) {
		error_msg("The sequences for '".$sequences{$seq}.
			"' and '$id' are identical.", 1);
	} else {
		$sequences{$seq} = $id;
	}
	if ($seq =~ /^[\-\.\s]+$/m) {
		error_msg("All-gap sequence for '$id'.", 1);
	}
}

my $phy_file = "$control_dir/$seq_base.phy";
if ($overwrite || !-e $phy_file) {
	save_to_file(to_phylip_s($seqs, $order), $phy_file);
} else {
	error_msg("The sequence file '$phy_file' already exists from a previous ".
		"analysis, and will not be overwritten. Either use the --replace ".
		"flag, or provide a different --out_dir.", 1);
}

my $nh_file = "$control_dir/$seq_base.nh";
if ($overwrite || !-e $nh_file) {
	save_to_file(read_from_file($tree_file), $nh_file);
} else {
	error_msg("The tree file '$nh_file' already exists from a previous ".
		"analysis, and will not be overwritten. Either use the --replace ".
		"flag, or provide a different --out_dir.", 1);
}

foreach my $model_type (sort keys %models) {
	my $phase_seq_file = "$control_dir/$seq_base\_$model_type.phy";
	if ($overwrite || !-e $phase_seq_file) {
		my $phase_class = undef;
		my $phase_structure = undef;
		if ($model_type =~ /din/) {
			$phase_structure = $structure;
		} else {
			$phase_class = $stem_pattern;
		}
		
		phase_seq_file(
			$phy_file,
			$phase_class,
			$phase_structure,
			$phase_seq_file
		);
	}
}
# Generate Sequence Files - End
################################################################################

################################################################################
# Set up Mixed Models - Start
# (The single nucleotide model is included in these
# 'mixed' models, for the sake of easier programming.)
my %mixed_models;
foreach my $model_type (sort keys %models) {
	my $nuc_models = () = $model_type =~ /nuc/g;
	my $dinuc_models = () = $model_type =~ /din/g;
	my $combinations =
		scalar(@nuc_models)**$nuc_models * scalar(@dinuc_models)**$dinuc_models;
	
	my @phase_models;
	# It'll always be a nucleotide model first.
	for (my $i = 0; $i < $combinations/scalar(@nuc_models); $i++) {
		for my $nuc_model (@nuc_models) {
			(my $model_name = $model_type) =~ s/nuc/$nuc_model/;
			push @phase_models, $model_name;
		}
	}
	
	if ($nuc_models > 1) {
		my $cursor = 0;
		for my $nuc_model (@nuc_models) {
			for (my $i = 0; $i < $combinations/scalar(@nuc_models); $i++) {
				$phase_models[$i+$cursor] =~ s/nuc/$nuc_model/;
			}
			$cursor += scalar(@nuc_models);
		}
	}
	
	if ($dinuc_models > 0) {
		my $cursor = 0;
		for my $dinuc_model (@dinuc_models) {
			for (my $i = 0; $i < $combinations/scalar(@dinuc_models); $i++) {
				$phase_models[$i+$cursor] =~ s/din/$dinuc_model/;
			}
			$cursor += scalar(@nuc_models);
		}
	}
	
	$mixed_models{$model_type} = \@phase_models;
}
# Set up Mixed Models - End
################################################################################

################################################################################
# Calculate Likelihoods - Start
my @columns = (
	'Model Name', 'Free Parameters', 'lnL', 'lnL Adjustment',
	'AIC', 'AICc', 'Delta AIC', 'Delta AICc');
my $results = join("\t", @columns);

my (%results, @aic_list, @aicc_list, $aic_counter);

foreach my $model_type (sort keys %models) {
	my $phase_seq_file = "$control_dir/$seq_base\_$model_type.phy";
	my $template_file = "$template_dir/".$models{$model_type};
	
	foreach my $phase_model (@{$mixed_models{$model_type}}) {
		my $ctl_file = "$control_dir/$seq_base\_$phase_model.ctl";
		my $out_base = "$results_dir/$seq_base\_$phase_model";
		my $out_file = "$out_base.out";
		
		if ($overwrite || !-e $out_file) {
			print "Executing model '$phase_model'.\n" if $verbose;
			(my $base_models = $phase_model) =~ s/\+G//g;
			my @models = split(/_/, $base_models);
			
			phase_ctl_file(
				$template_file,
				$phase_seq_file,
				$out_base,
				undef,
				\@models,
				$nh_file,
				undef,
				$ctl_file);
			phase_ml($ctl_file, 'fatal', $phase_optimizer);
		} else {
			print "Existing results for model '$phase_model'.\n" if $verbose;
		}
		
		if (-e $out_file) {
			my %phase_results = phase_parse($out_file);
			unless ($phase_results{'ln_likelihood'}) {
				print STDERR "PHASE failed for model '$phase_model'.\n";
				unlink $out_file;
				next;
			}
			if ($phase_model =~ /RNA7/) {
				my @freq_types = ('equal', 'empirical');
				foreach my $freq_type (@freq_types) {
					my $id = "$phase_model.$freq_type";
					my $lnL = $phase_results{'ln_likelihood'} + $phase_results{"adjustment_$freq_type"};
					my $parameters = $phase_results{'parameters'};
					$parameters += 9 if $freq_type eq 'empirical';
					$results{$id}{'lnL'} = $lnL;
					$results{$id}{'lnL Adjustment'} = $phase_results{"adjustment_$freq_type"};
					$results{$id}{'Free Parameters'} = $parameters;
					$results{$id}{'AIC Index'} = $aic_counter++;
					$results{$id}{'Model Type'} = $model_type;
					push @aic_list, aic($lnL, $parameters);
					push @aicc_list, aicc($lnL, $parameters, $sample_size);
				}
			} else {
				my $lnL = $phase_results{'ln_likelihood'};
				my $parameters = $phase_results{'parameters'};
				$results{$phase_model}{'lnL'} = $lnL;
				$results{$phase_model}{'Free Parameters'} = $parameters;
				$results{$phase_model}{'AIC Index'} = $aic_counter++;
				$results{$phase_model}{'Model Type'} = $model_type;
				push @aic_list, aic($lnL, $parameters);
				push @aicc_list, aicc($lnL, $parameters, $sample_size);
			}
		}
	}
}

my @delta_aic = delta_aic(\@aic_list);
my @delta_aicc = delta_aic(\@aicc_list);

foreach my $model (sort keys %results) {
	my $aic_index = $results{$model}{'AIC Index'};
	$results .= "\n$model".
				"\t".$results{$model}{'Free Parameters'}.
				"\t".$results{$model}{'lnL'}.
				"\t";
	$results .= $results{$model}{'lnL Adjustment'} || "";
	$results .= "\t".$aic_list[$aic_index].
				"\t".$aicc_list[$aic_index].
				"\t".$delta_aic[$aic_index].
				"\t".$delta_aicc[$aic_index];
}
save_to_file($results, "$out_dir/$seq_base\_results.txt");
# Calculate Likelihoods - End
################################################################################

################################################################################
# Report Best Model - Start
my ($index, $best_model);
foreach my $i (0..$#delta_aicc) {
	if ($delta_aicc[$i] == 0) {
		$index = $i;
		last;
	}
}
foreach my $model (sort keys %results) {
	if ($results{$model}{'AIC Index'} == $index) {
		$best_model = $model;
		last;
	}
}
print "\nModel $best_model is the best fit to the data (using AICc).\n".
		"Full details are saved in '$out_dir/$seq_base\_results.txt'.\n";

my $mod_file = "$results_dir/$seq_base\_$best_model.mod";
$mod_file =~ s/(\.equal|\.empirical)//;
my $mod = read_from_file($mod_file);
if ($mod =~ /substitution rate ratio = 0\.0+\s/m) {
	error_msg(
		"Zero-valued rate of change between partitions in the best ".
		"model; check your alignments and other results, to verify ".
		"that performing these analyses is sensible.\n", 1);
}
# Report Best Model - End
################################################################################

################################################################################
# Tree Search - Start
if ($tree_search || $exhaustive_tree_search) {
	if ($seq_count < 4) {
		error_msg(
			"There are fewer than 4 sequences in the alignment, ".
			"so tree search will not be conducted.\n", 1);
	} else {
		my @tree_search_models;
		if (!$exhaustive_tree_search) {
			push @tree_search_models, $nuc_models[0].'+G';
			if ($best_model ne $nuc_models[0].'+G') {
				push @tree_search_models, $best_model;
			}
		} else {
			@tree_search_models = sort keys %results;
		}
		
		foreach my $ts_model (@tree_search_models) {
			my $model_type = $results{$ts_model}{'Model Type'};
			(my $trimmed = $ts_model) =~ s/(\.equal|\.empirical)//;
			
			my $phase_seq_file = "$control_dir/$seq_base\_$model_type.phy";
			my $template_file = "$template_dir/".$models{$model_type};
			$template_file =~ s/ml$/tree.mcmc/;
			
			my $start_mod_file = "$results_dir/$seq_base\_$trimmed.mod";
			my $ctl_file = "$control_dir/$seq_base\_$trimmed.mcmc.ctl";
			my $out_base = "$results_dir/$seq_base\_$trimmed.mcmc";
			my $out_file = "$out_base.out";
			my $tree_con_file = "$out_dir/$seq_base.consensus_$ts_model.nh";
			
			if ($overwrite || !-e $out_file) {
				print "\nExecuting tree search with model '$ts_model'.\n" if $verbose;
				
				# Since we've got estimates of model parameters, might as well
				# use them as a starting point.
				(my $base_models = $trimmed) =~ s/\+G//g;
				my @models = split(/_/, $base_models);
				
				phase_ctl_file(
					$template_file,
					$phase_seq_file,
					$out_base,
					$start_mod_file,
					\@models,
					undef,
					undef,
					$ctl_file);
				
				phase_mcmc($ctl_file, 'fatal', $phase_mcmc);
			} else {
				print "\nExisting tree search with model '$ts_model'.\n" if $verbose;
			}
			
			if ($overwrite || !-e $tree_con_file) {
				my ($consensus, undef) = phase_mcmc_tree($ctl_file, 'fatal', $phase_mcmcsummarize);
				save_to_file($consensus, $tree_con_file);
			}
			
			print "\nMajority-rule consensus topology calculated with model ".
					"$ts_model, saved in '$tree_con_file'.\n";
		}
	}
}
# Tree Search - End
################################################################################

################################################################################
sub error_msg {
	my ($msg, $warn) = @_;
	print STDERR $warn ? "\nWARNING: " : "\nERROR: ";
	print STDERR "$msg\n";
	exit unless $warn;
}
################################################################################


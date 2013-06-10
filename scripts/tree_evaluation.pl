#!/usr/bin/perl
use strict;
use warnings;

use File::Path qw (make_path);
use Getopt::Long;

use PhaseUtils
	qw (read_from_file read_delim save_to_file from_fasta to_phylip_s);

=pod

=head1 Name

tree_evaluation.pl

=head1 Description

For a given alignment and set of trees, run Phyml (fixed topology),
to generate sitewise likelihoods. Then feed the sitewise likelihoods
to Consel to do AU/SH-tests.

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
my(	$help, $verbose,
	$alignment_file, $trees_file, $ids_file,
	$out_dir, $phyml, $consel, $makermt, $catpv,
);

GetOptions(
	'help' => \$help,
	'verbose' => \$verbose,
	'alignment_file=s' => \$alignment_file,
	'trees_file=s' => \$trees_file,
	'ids_file=s' => \$ids_file,
	'out_dir:s' => \$out_dir,
	'phyml:s' => \$phyml,
	'consel_consel:s' => \$consel,
	'consel_makermt:s' => \$makermt,
	'consel_catpv:s' => \$catpv,
);

sub usage {
  print qq~

tree_evaluation.pl: SH/AU-Tests with Consel, via Phyml

This script performs maximum likelihood inference of model parameters and
branch lengths for a set of fixed topologies with the GTR+G model. Phyml
provides sitewise likelihoods, and these are  passed to Consel to calculate
the significance of any differences between the trees.

* Mandatory parameters:
  --a[lignment_file] <file>
  --t[rees_file] <file>
  Two input files are required, an alignment in Fasta format, and a set of
  trees in Newick format in a single file, separated by newlines. Readseq
  is included in this distribution of PHASE (../readseq), so please use
  that to convert to Fasta format, if necessary.

* Optional parameters:
  --i[ds_file] <file>
  By default, the trees in the 'trees_file' will be numbered from 1 to n in
  the output file, according to their order in the file. Mapping between trees
  and results can be awkward when there are many trees, so a file can be given
  that associates each tree with a name, and this name will then be used in the
  output. Each row in the 'ids_file' correponds to the same row in the
  'trees_file', and the number of rows must match.

  --o[ut_dir] <dir>
  If an output directory is not specified, results are saved in the current
  directory. The Consel and Phyml output is not saved, as these files can
  be very large.

  --p[hyml] <file>
  This script will look for a Phyml executable file called 'phyml' in your
  path. If not found then a path must be provided as a parameter.

  --consel_co[nsel] <file>
  --consel_m[akermt] <file>
  --consel_ca[tpv] <file>
  This script will look for Consel executable files called 'consel',
  'makermt' and 'catpv' in your  path. If not found then a path must
  be provided as a parameter.

* Optional flags:
  --v[erbose]
  Display informational messages about what the script is doing.

~;
  print 'Author: James Allen (james@monkeyshines.co.uk)'."\n\n";
  print "Usage: $0 --align <file> -str <file> -tree <file> [options]\n";
  print qq~
  Please provide values for:
    --a[lignment_file]  Fasta-format alignment
    --t[rees_file]      Newick-format trees
  
  Optional parameters:
    --i[ds_file]        default <ordinal numbers>
    --o[ut_dir]         default '.'
    --p[hyml]           default 'phyml'
    --consel_co[nsel]   default 'consel'
    --consel_m[akermt]  default 'makermt'
    --consel_ca[tpv]    default 'catpv'

  Optional flags:
    --h[elp]     View this help page
    --v[erbose]  Display progress on screen
  
  ~;
  print "\n";
  exit;
}

usage() if (defined $help);

unless ( $alignment_file && $trees_file ) {
	usage;
}

$verbose = 0 unless $verbose;
# Arguments - End
################################################################################

################################################################################
# Sort Out Files and Directories - Start
if (!-e $alignment_file) {
	error_msg("Sequence file '$alignment_file' does not exist.");
}
if (!-e $trees_file) {
	error_msg("Tree file '$trees_file' does not exist.");
}
my @ids;
if ($ids_file) {
	if(-e $ids_file) {
		@ids = split(/\n/, read_from_file($ids_file));
		my @trees = split(/\n/, read_from_file($trees_file));
		if (scalar(@ids) != scalar(@trees)) {
			error_msg("Mismatch between the number of trees in '$trees_file' and the number of ids in '$ids_file'");
		}
	} else {
		error_msg("IDs file '$ids_file' does not exist.");
	}
}

$out_dir = '.' unless $out_dir;
if (!-e $out_dir) {
	make_path($out_dir);
	print "Writing output to new directory, '$out_dir'\n" if $verbose;
} else {
	print "Writing output to existing directory, '$out_dir'\n" if $verbose;
}

if (!defined($phyml)) {
	$phyml = "phyml";
}
my $phyml_out = `$phyml --quiet 2>&1`;
if (!defined($phyml_out) || $phyml_out !~ /PhyML/m) {
	error_msg("Could not find the Phyml program ($phyml).\n".
			"Use the --phyml option to specify its location.\n");
}

if (!defined($consel)) {
	$consel = "consel";
}
my $consel_out = `$consel --help 2>&1`;
if (!defined($consel_out) || $consel_out !~ /consel\.c/m) {
	error_msg("Could not find the Consel program ($consel).\n".
			"Use the --consel_consel option to specify its location.\n");
}
if (!defined($makermt)) {
	$makermt = "makermt";
}
my $makermt_out = `$makermt --help 2>&1`;
if (!defined($makermt_out) || $makermt_out !~ /makermt\.c/m) {
	error_msg("Could not find the Consel program ($makermt).\n".
			"Use the --consel_makermt option to specify its location.\n");
}
if (!defined($catpv)) {
	$catpv = "catpv";
}
my $catpv_out = `$catpv --help 2>&1`;
if (!defined($catpv_out) || $catpv_out !~ /Terminated/m) {
	error_msg("Could not find the Consel program ($catpv).\n".
			"Use the --consel_catpv option to specify its location.\n");
}
# Sort Out Files and Directories - End
################################################################################

################################################################################
# Execute CONSEL - Start
sub consel {
	my ($sitewise_file) = @_;
	my %consel;
	my @cols = qw(rank item obs au np bp pp kh sh wkh wsh);
	
	(my $base_name = $sitewise_file) =~ s/\.\w+$//;
	my $makermt = `makermt --phyml $sitewise_file`;
	if ($makermt =~ /exit normally/) {
		my $consel = `consel $base_name.rmt`;
		if ($consel =~ /exit normally/) {
			my $catpv = `catpv $base_name.pv`;
			if ($catpv !~ /Terminated/) {
				# The output table generated by catpv is great for human readability,
				# but less so for computers, so reformat as tab-delimited.
				$catpv =~ s/^\n//gm;
				$catpv =~ s/^\# reading.*\n//gm;
				$catpv =~ s/[\#\|]//gm;
				$catpv =~ s/^ +//gm;
				$catpv =~ s/ +/\t/gm;
				
				save_to_file($catpv, "$base_name.consel");
				%consel = read_delim("$base_name.consel", \@cols, 'item');
			} else {
				print "CONSEL: catpv execution failed ($sitewise_file)\n";
			}
		} else {
			print "CONSEL: consel execution failed ($sitewise_file)\n";
		}
	} else {
		print "CONSEL: makermt execution failed ($sitewise_file)\n";
	}
	
	unlink "$base_name.rmt";
	return wantarray ? %consel : \%consel;
}
# Execute CONSEL - End
################################################################################

################################################################################
# Calculate Tree Significance - Start
sub main {
	my $tree_eval_file = "$out_dir/tree_evaluation.txt";
	save_to_file("Tree\tRank\tSH\tAU", $tree_eval_file);
	
	my $phy_file = "$alignment_file.phy";
	unless (-e $phy_file) {
		my ($seqs, $order) = from_fasta(read_from_file($alignment_file));
		save_to_file(to_phylip_s($seqs, $order), $phy_file);
	}
	`$phyml -i $phy_file -m GTR -u $trees_file -b 0 -o lr --print_site_lnl`;
	my $sitewise_file = "$phy_file\_phyml_lk.txt";
	if (-e $sitewise_file) {
		my %consel = consel("$phy_file\_phyml_lk.txt");
		
		my $tree_eval;
		foreach my $item (sort keys %consel) {
			if ($ids_file) {
				$tree_eval .= "\n".$ids[$item-1];
			} else {
				$tree_eval .= "\n".$item;
			}
			$tree_eval .= "\t".$consel{$item}{'rank'};
			$tree_eval .= "\t".$consel{$item}{'sh'};
			$tree_eval .= "\t".$consel{$item}{'au'};
		}
		save_to_file($tree_eval, $tree_eval_file, 'append');
		
		# Tidy up results, to make the directory contents manageable.
		unlink "$phy_file\_phyml_lk.ci";
		unlink "$phy_file\_phyml_lk.consel";
		unlink "$phy_file\_phyml_lk.pv";
		unlink "$phy_file\_phyml_lk.vt";
		unlink "$phy_file\_phyml_lk.txt";
		unlink "$phy_file\_phyml_stats.txt";
		unlink "$phy_file\_phyml_tree.txt";
		unlink "$phy_file";
	} else {
		print "PhyML failed to generate sitewise likelihoods\n";
	}
}
# Calculate Tree Significance - End
################################################################################

################################################################################
sub error_msg {
	my ($msg, $warn) = @_;
	print STDERR $warn ? "\nWARNING: " : "\nERROR: ";
	print STDERR "$msg\n";
	exit unless $warn;
}
################################################################################

main();


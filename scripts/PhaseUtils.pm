package PhaseUtils;

use warnings;
use strict;
use Carp;

use Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
	read_from_file read_delim save_to_file
	from_fasta to_phylip_s
	parse_structure sample_size pair_gaps
	phase_seq_file phase_ctl_file
	phase_ml phase_mcmc phase_parse phase_mcmc_tree
	aic aicc delta_aic
);

=pod

=head1 Name

PhaseUtils

=head1 Description

A module containing functions that are useful for scripting
execution of the PHASE programs.

=head1 Functions

=head2 Sequence File Parsing

The world hardly needs more code for parsing sequence data, but for the sake
of having everything in one place, here are functions to go to and from Fasta
and Phylip format. Ideally I'd just use Readseq, but it truncates taxa names,
which can cause problems (PHASE is happy with the 'relaxed' format that PAML
uses, where names can be longer than 10 characters as long as there are two
trailing spaces). Also, it's useful to have basic functionality if Readseq
doesn't work/isn't available etc.

=over

=item B<($seqs, $order) = from_fasta($data)>

Convert Fasta-formatted $data into a hashref of IDs and sequences ($seqs),
and a hashref that tracks the order of the sequences in the file ($order).

=item B<$data = to_fasta($seqs[, $order])>

Convert a hashref of IDs and sequences ($seqs) into Fasta-formatted $data,
optionally using a hashref defining the order of the sequences ($order).

=item B<($seqs, $order) = from_phylip_s($data)>

Convert Phylip-sequential-formatted $data into a hashref of IDs and
sequences ($seqs), and a hashref that tracks the order of the sequences
in the file ($order). The function is not strict about the length of
taxa names, and considers the first instance of whitespace in a line
to delineate taxon and sequence.

=item B<$data = to_phylip_s($seqs[, $order])>

Convert a hashref of IDs and sequences ($seqs) into Phylip-sequential-
formatted $data, optionally using a hashref defining the order of the
sequences ($order). The function is not strict about the length of
taxa names, and uses two spaces to delineate taxon and sequence.

=back

=head2 Structure Parsing

=over

=item B<(%|$)stems = parse_structure($structure)>

From a dot-bracket $structure, generate a hash or hashref with
keys that represent the position of opening brackets, and values
that represent the position of the appropriate closing bracket.

=item B<$size = sample_size($seqs, $structure[, $sample_type])>

Calculation of corrected AIC (AICc) values requires a 'sample size'.
It is not clear what this means in the context of sequence alignments,
but two commonly used approximations are the number of characters
and the number of sites. Use $sample_type to specify either 'characters'
(the default) or 'sites'. $seqs is a hashref of IDs and sequences,
and $structure is a dot-bracket structure; the latter is required so
that a base-pair counts as 1, rather than as 2 separate bases.

=item B<C<undef> = pair_gaps($seqs, $structure)>

Sometimes a structure specifies a base-pair that, in one or more taxa
consists of a gap and a base, eg A-. The most practical and sensible
course of action is to replace this with a pair of gaps (see the PHASE
manual for further discussion). The hashref $seqs, of IDs and sequences,
is edited to ensure that all sequences have either two bases or two gap
characters at base-pairs defined by the dot-bracket $structure.

=back

=head2 RNA Model Analysis with PHASE

=over

=item B<$new_seq_file = phase_seq_file($seq_file[, $class][, $structure][, $new_seq_file])>

PHASE requires sequence data in a modified version of Phylip format (see
the PHASE manual for details). This script modifies a standard Phylip format
file ($seq_file) appropriately. The $class parameter, if given, defines a
partition of the data in a mixed model, and should be a space-separated string
of numbers indicating the partition of sites, e.g. '1 1 2 2 2 1 1'. For RNA
sequences a secondary structure is required, and $structure should contain
this in standard dot-bracket notation. The modified sequence file is saved as
$new_seq_file, if given, or saved as '$seq_file' with a unique numbered suffix.

=item B<$ctl_file = phase_ctl_file($template_file, $seq_file, $out_file[, $mod_file][, $models][, $tree_file][, $outgroup][, $ctl_file])>

$template_file is a path to a PHASE control file, into which are inserted
the $seq_file and $out_file parameters. Starting model parameters can be
specified with $mod_file. By default REV is used for DNA models, and 7A
for RNA models; to use a different model (or models, in the case of a
mixed model), supply an arrayref, $models. For ML analysis of a fixed tree,
use $tree_file. PHASE always requires an outgroup for formatting trees;
the choice doesn't affect the results, and by default the last taxa in
$seq_file is used. The modified sequence file is saved as $ctl_file, if
given, or saved as '$seq_file' with a unique numbered suffix.

Examples of template files, for a range of common evolutionary models,
are provided in ./templates/.

=item B<$out_file = phase_ml($ctl_file[, $fatal][, $program_path])>

$ctl_file is a path to a PHASE control file, which will be supplied to either
the 'mlphase' or the 'optimizer' program, depending on whether or not a tree is
specified in the control file. By default, the Perl code will abort execution
if an error is detected when running PHASE; set $fatal to 0 (zero) to display
a warning and continue execution. The location of the executable file can be
specified ($program_path); by default the programs are assumed to be in the
user's execution path. The output is stored as specified by the control file;
the output file path is returned.

=item B<[$|%]results = phase_ml_parse($out_file)>

$out_file is a path to a PHASE results file from an ML analysis (ie either
'mlphase' or 'optimizer').. The %results hash (or hash reference, if called
in scalar context) extracts the key information from the PHASE output,
namely the following keys: 'ln_likelihood', 'tree_length', 'tree', and
'parameters' (the number of free parameters in the model). For 7-state
models, two more values, 'adjustment_equal' and 'adjustment_empirical',
are available, representing the adjustment required to make the
log-likelihood comparable with a 16-state log-likelihood, for frequencies
of individual mismatch states calculated either as equal or empirically.

=item B<$out_file = phase_mcmc($ctl_file[, $fatal][, $program_path])>

$ctl_file is a path to a PHASE control file, which will be supplied to
the 'mcmcphase' program. By default, the Perl code will abort execution
if an error is detected when running PHASE; set $fatal to 0 (zero) to display
a warning and continue execution. The location of the executable file can be
specified ($program_path); by default the program is assumed to be in the
user's execution path. The output is stored as specified by the control file;
the output file path is returned.

=item B<($support_tree, $blength_tree) = phase_mcmc_tree($ctl_file[, $fatal][, $program_path])>

Once an MCMC analysis has been conducted, a consensus tree can be created.
$ctl_file is a path to the PHASE control file that was used to perform the
MCMC, and it will be supplied to the 'mcmcsummarize' program. By default,
the Perl code will abort execution if an error is detected when running
PHASE; set $fatal to 0 (zero) to display a warning and continue execution.
The location of the executable file can be specified ($program_path); by
default the program is assumed to be in the user's execution path. Two trees
are returned: $support_tree is annotated with posterior probabilities;
$blength_tree has branch lengths, and should be used with caution, if at all
(see PHASE manual for details, section headed 'Using mcmcsummarize').

=back

=back

=head2 Calculate AIC

=over

=item B<$aic = aic($lnL, $parameters)>

Calculate AIC from a log-likelihood $lnL and the number of free $parameters.

=item B<$aicc = aicc($lnL, $k, $n)>

Calculate the corrected AIC from a log-likelihood $lnL, the number of free
parameters $k, and the sample size $n.

=item B<(@|$)aic_delta = delta_aic($aic_list)>

For an arrayref of AIC values $aic_list, calculate a corresponding array(ref)
of delta AIC values, defined as AIC-min(AIC).

=back

=head2 Random Numbers

=over

=item B<$seed = random_seed()>

Generate an 8-digit pseudo-random number.

=item B<$timestamp_id = timestamp_id([$max_random])>

Return a string consisting of a timestamp and a random number between 0 and
1000 (by default - specify $max_random to use a different upper limit). This
function is useful for generating names for temporary files.

=back

=head2 File IO

=over

=item B<$data = read_from_file($file)>

Read data from the given file location into a string.

=item B<C<undef> = save_to_file($data, $file[, $append])>

Save $data in the location given by $file; if $append is true, then
append to, rather than overwrite, any existing data in that file.

=back

=head2 Numeric Calculations

=over

=item B<$min = min($value1, $value2[, $value3[, ...]])>

=item B<$min = min(@values)>

Find the minimum of two or more values, or the minimum value in an array.

=item B<$max = max($value1, $value2[, $value3[, ...]])>

=item B<$max = max(@values)>

Find the maximum of two or more values, or the maximum value in an array.

=item B<$sum = sum($value1, $value2[, $value3[, ...]])>

=item B<$sum = sum(@values)>

Find the sum of two or more values, or the sum of all values in an array.

=back

=head2 Script Execution

=over

=item B<$out = execute_check($execute, $check[, $check_success][, $program][, $fatal])>

Execute shell scripts/programs, and check for errors. $execute should
consist of a program name and any parameters, but should not redirect either
the standard error or output streams, as the function will do this and place
both in the $out return value. To test whether a program executed without
errors, the $check parameter should contain text that appears in the program's
error message; if this text is found in the output, the error is passed on to
the user. Alternatively, text that appears in the valid output of the program
can be provided as the $check value, if $check_success is true. The $program
parameter is used to report which program failed to the user. If the $fatal
parameter is true then execution ceases after the error is reported, otherwise
it will continue after displaying a warning (default is the latter).

=back

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
# SEQUENCE FILE PARSING
################################################################################
sub from_fasta {
	my ($in) = @_;
	my (%seqs, %order);
	%seqs = $in =~ /^>(\S+)[^\n]*\n([^>]*)/gms;
	foreach my $id (keys %seqs) {
		$seqs{$id} =~ s/[\n\r ]//gms;
	}
	my @ids = $in =~ /^>(\S+)/gms;
	for (my $i=0; $i<scalar(@ids); $i++) {
		$order{$ids[$i]} = $i+1;
	}
	return (\%seqs, \%order);
}

sub to_fasta {
	my ($seqs, $order) = @_;
	my $out;
	my @ids = keys(%$seqs);
	unless ($order) {
		for (my $i=0; $i<scalar(@ids); $i++) {
			$$order{$ids[$i]} = $i+1;
		}
	}
	foreach my $id (sort {$$order{$a} <=> $$order{$b}} keys %$order) {
		if (exists($$seqs{$id})) {
			(my $seq = $$seqs{$id}) =~ s/([\w\-\.]{1,80})/$1\n/gms;
			$out .= ">$id\n$seq\n";
		}
	}
	return $out;
}

sub from_phylip_s {
	my ($in) = @_;
	my (%seqs, %order);
	my $counter = 0;
	my @lines = split(/\n/, $in);
	foreach my $line (@lines) {
		if ($line =~ /^(\S+)\s+(\S+)$/) {
			$seqs{$1} .= $2;
			unless (exists($order{$1})) {
				$order{$1} = ++$counter;
			}
		}
	}
	return (\%seqs, \%order);
}

sub to_phylip_s {
	my ($seqs, $order) = @_;
	my $out;
	my @ids = keys(%$seqs);
	my $seq_count = scalar(keys(%$seqs));
	my $seq_length = length($$seqs{$ids[0]});
	$out = "  $seq_count $seq_length\n";
	unless ($order) {
		for (my $i=0; $i<scalar(@ids); $i++) {
			$$order{$ids[$i]} = $i+1;
		}
	}
	foreach my $id (sort {$$order{$a} <=> $$order{$b}} keys %$order) {
		if (exists($$seqs{$id})) {
			my $seq = $$seqs{$id};
			$out .= "$id  $seq\n";
		}
	}
	return $out;
}

sub from_phylip_i {
	my ($data, $relaxed) = @_;
	$relaxed = 1 unless defined $relaxed;
	my (%seqs, %order, @ids);
	my $counter = 0;
	
	my @lines = split(/\n[\n \t]*/, $data);
	my $info = shift @lines;
	my ($count, $length) = $info =~ /(\d+)\s+(\d+)/;
	
	foreach my $i (1..$count) {
		my $line = shift @lines;
		my ($id, $seq);
		if ($relaxed) {
			($id, $seq) = $line =~ /^(\S+)\s+(.+)/;
		} else {
			($id, $seq) = $line =~ /^(.{10})(.+)/;
		}
		$id =~ s/\s+$//g;
		$seqs{$id} = $seq;
		$order{$id} = ++$counter;
		push @ids, $id;
	}
	
	while (@lines) {
		foreach my $i (1..$count) {
			$seqs{$ids[$i-1]} .= shift @lines;
		}
	}
	
	foreach my $id (keys %seqs) {
		$seqs{$id} =~ s/\s+//g;
	}
	
	return (\%seqs, \%order);
}

sub to_phylip_i {
	my ($seqs, $order, $relaxed) = @_;
	$relaxed = 1 unless defined $relaxed;
	my $out;
	my @ids = keys(%$seqs);
	my $seq_count = scalar(keys(%$seqs));
	my $seq_length = length($$seqs{$ids[0]});
	$out = "  $seq_count $seq_length\n";
	unless ($order) {
		for (my $i=0; $i<scalar(@ids); $i++) {
			$$order{$ids[$i]} = $i+1;
		}
	}
	foreach my $id (sort {$$order{$a} <=> $$order{$b}} keys %$order) {
		(my $seq = $$seqs{$id}) =~ s/(.{1,10})/$1 /gm;
		$seq =~ s/(.{1,54}) /$1\n/gm;
		my @seq = split(/\n/, $seq);
		
		if ($relaxed) {
			$out .= "$id ";
		} else {
			$out .= sprintf "%-10.10s", $id;
		}
		$out .= shift @seq;
		$out .= "\n";
		
		$$seqs{$id} = \@seq;
	}
	
	my $chunks = scalar(@{$$seqs{$ids[0]}});
	foreach my $i (1..$chunks) {
		$out .= "\n";
		foreach my $id (sort {$$order{$a} <=> $$order{$b}} keys %$order) {
			$out .= "          ";
			$out .= shift @{$$seqs{$id}};
			$out .= "\n";
		}
	}
	
	return $out;
}
################################################################################

################################################################################
# STRUCTURE PARSING
################################################################################
sub parse_structure {
	my ($structure) = @_;
	my %pairs = ('(' => ')', '<' => '>', '[' => ']', 'A' => 'a', 'B' => 'b');
	my %pairs_rev = (')' => '(', '>' => '<', ']' => '[', 'a' => 'A', 'b' => 'B');
	
	my @structure = split(//, $structure);
	my $position = 0;
	my (%stems, %stack);
	foreach my $base (@structure) {
		if ($base !~ /[\-\.]/) {
			if (exists($pairs{$base})) {
				push @{$stack{$base}}, $position;
			} elsif (exists($pairs_rev{$base})) {
				my $partner = $pairs_rev{$base};
				my $paired_position = pop @{$stack{$partner}};
				$stems{$paired_position} = $position;
			} else {
				die "Unrecognised character '$base'";
			}
		}
		$position++;
	}
	
	return wantarray ? %stems : \%stems;
}

sub sample_size {
	my ($seqs, $structure, $sample_type) = @_;
	$sample_type = 'characters' unless $sample_type;
	my $size = 0;
	
	if ($sample_type =~ /sites/i) {
		my $unpaired = () = $structure =~ /\./g;
		my $paired = length($structure) - $unpaired;
		my $size = $unpaired + ($paired/2);
	} elsif ($sample_type =~ /characters/i) {
		my %stems = parse_structure($structure);
		foreach my $id (keys %$seqs) {
			my @seq = split(//, $$seqs{$id});
			my $pairs = 0;
			while (my ($pos1, $pos2) = each(%stems)) {
				if ($seq[$pos1] ne '-' && $seq[$pos2] ne '-') {
					$pairs++;
				}
			}
			$size += () = $$seqs{$id} =~ /[^\-]/g;
			$size -= $pairs;
		}
	} else {
		die "Unrecognised option '$sample_type' for sample type.";
	}
	return $size;
}

sub pair_gaps {
	my ($seqs, $structure) = @_;
	my %stems = parse_structure($structure);
	foreach my $id (keys %$seqs) {
		my @seq = split(//, $$seqs{$id});
		while (my ($pos1, $pos2) = each(%stems)) {
			if ($seq[$pos1] eq '-' && $seq[$pos2] ne '-') {
				$seq[$pos2] = '-';
			} elsif ($seq[$pos2] eq '-' && $seq[$pos1] ne '-') {
				$seq[$pos1] = '-';
			}
		}
		$$seqs{$id} = join('', @seq);
	}
	return;
}
################################################################################

################################################################################
# RNA MODEL ANALYSIS WITH PHASE
################################################################################
sub phase_seq_file {
	my ($seq_file, $class, $structure, $new_seq_file) = @_;
	
	unless ($new_seq_file) {
		my $timestamp_id = timestamp_id();
		($new_seq_file = $seq_file) =~ s/(\.\w+)$/_$timestamp_id$1/;
	}
	
	# The data file needs a flag for the data type,
	# and, optionally, a structure and/or class definitions.
	my $seq = read_from_file($seq_file);
	if ($structure) {
		$seq =~ s/^(\s*\d+\s+\d+)\n/$1 STRUCT\n$structure\n\n/m;
	} else {
		$seq =~ s/^(\s*\d+\s+\d+)\n/$1 DNA\n/m;
	}
	if ($class) {
		$seq .= "\n$class\n";
	}
	save_to_file($seq, $new_seq_file);
	
	return $new_seq_file;
}

sub phase_ctl_file {
	my ($template_file, $seq_file, $out_file, $mod_file, $models, $tree_file, $outgroup, $ctl_file) = @_;
	
	my $ctl = read_from_file($template_file);
	$ctl =~ s/(Data file *= *)/$1$seq_file/;
	$ctl =~ s/(Output file *= *)/$1$out_file/;
	
	if ($mod_file) {
		$ctl =~ s/(Starting model parameters file *= *)/$1$mod_file/;
	} else {
		$ctl =~ s/(Starting model parameters file *= *)/#$1/;
	}
	
	# The templates default to REV for DNA models, and 7A for RNA models.
	if ($models) {
		if (scalar(@$models) == 1) {
			$ctl =~ s/(Model *= *)\S*/$1$$models[0]/;
		} else {
			for (my $i=1; $i<= scalar(@$models); $i++) {
				$ctl =~ s/(\{MODEL$i\}\nModel *= *)\S*/$1$$models[$i-1]/;
			}
		}
	}
	
	if ($tree_file) {
		$ctl =~ s/(Tree file *= *)/$1$tree_file/;
	} else {
		$ctl =~ s/(Tree file *= *)/#$1/;
	}
	
	unless ($outgroup) {
		# Need to get outgroup from seq_file - it's all a
		# bit arbitrary, so just pick the last species.
		my $seq = read_from_file($seq_file);
		($outgroup) = $seq =~ /^(\S+)\s+\S+[\s\d]*\Z/m;
	}
	$ctl =~ s/(Outgroup *= *)/$1$outgroup/;
	
	my $random_seed = random_seed();
	# Make the 10-digit random seed sufficiently short to fit in a C++ long int variable.
	$random_seed =~ s/\d$//;
	$ctl =~ s/(Random seed *= *)/$1$random_seed/;
	
	unless ($ctl_file) {
		my $timestamp_id = timestamp_id();
		($ctl_file = $seq_file) =~ s/(\.\w+)$/_$timestamp_id.ctl/;
	}
	
	save_to_file($ctl, $ctl_file);
	return $ctl_file;
}

sub phase_ml {
	my ($ctl_file, $fatal, $program) = @_;
	die "A control file is required" unless $ctl_file;
	$fatal = "fatal" unless defined $fatal;
	
	my $ctl = read_from_file($ctl_file);
	
	# Extract output file name from the control file.
	my ($out_file) = $ctl =~ /Output file\s*=\s*(\S+)\n/ms;
	# PHASE prompts to overwrite an existing output file,
	# which is no good for batch-mode execution...
	unlink($out_file, "$out_file.out", "$out_file.mod", "$out_file.tre");
	
	# Use optimizer if given tree topology, mlphase otherwise, and delete
	# either the 'Search algorithm' or the 'Tree file' option accordingly.
	my ($tree_file) = $ctl =~ /Tree file\s*=\s*(\S+)\n/;
	
	if ($tree_file) {
		$program = "optimizer" unless $program;
		$ctl =~ s/Search algorithm\s*=\s*[^\n]*\n//;
		$out_file .= ".out";
	} else {
		$program = "mlphase" unless $program;
		$ctl =~ s/Tree file\s*=\s*[^\n]*\n//;
	}
	save_to_file($ctl, $ctl_file);
	
	execute_check("$program \"$ctl_file\"", 'maxLikelihood', 1, "PHASE", $fatal);
	
	return $out_file;
}

sub phase_parse {
	my ($out_file) = @_;
	die "A results file is required" unless $out_file;
	
	my $tree;
	my %results;
	my $out_data = read_from_file($out_file);
	
	($tree, $results{'ln_likelihood'}) = $out_data =~ /^([^\n]+)\nmaxLikelihood.+(\-\d+\.*\d+)/m;
	if ($tree) {
		my @branch_length = $tree =~ /(\d+\.*\d+)/;
		foreach my $length (@branch_length) {
			$results{'tree_length'} += $length;
		}
		$results{'tree'} = $tree;
	}
	($results{'parameters'}) = $out_data =~ /^free parameters.* (\d+)$/m;
	my ($adjustment_equal) = $out_data =~ /^likelihood adjustment \(equal frequencies\): (.+)$/m;
	my ($empirical_extra, $adjustment_empirical) = $out_data =~ /^likelihood adjustment \(empirical frequencies => \+(\d+) free parameters\): (.+)$/m;
	if (defined $adjustment_equal) {
		$results{'adjustment_equal'} = $adjustment_equal;
	}
	if (defined $adjustment_empirical) {
		$results{'adjustment_empirical'} = $adjustment_empirical;
		$results{'empirical_extra'} = $empirical_extra;
	}
	
	return wantarray ? %results : \%results;
}

sub phase_mcmc {
	my ($ctl_file, $fatal, $program) = @_;
	die "A control file is required" unless $ctl_file;
	$fatal = "fatal" unless defined $fatal;
	
	my $ctl = read_from_file($ctl_file);
	
	# Extract output file name from the control file.
	my ($out_file) = $ctl =~ /Output file\s*=\s*(\S+)\n/ms;
	# PHASE prompts to overwrite an existing output file,
	# which is no good for batch-mode execution...
	unlink($out_file, "$out_file.out", "$out_file.mod", "$out_file.tre");
	
	$program = "mcmcphase" unless $program;
	execute_check("$program \"$ctl_file\"", 'Best Tree', 1, "PHASE", $fatal);
	
	return $out_file;
}

sub phase_mcmc_tree {
	my ($ctl_file, $fatal, $program) = @_;
	die "A control file is required" unless $ctl_file;
	$fatal = "fatal" unless defined $fatal;
	
	$program = "mcmcsummarize" unless $program;
	execute_check("$program \"$ctl_file\"", 'Consensus', 1, "PHASE (mcmcsummarize)", $fatal);
	
	my $ctl = read_from_file($ctl_file);
	my ($out_file) = $ctl =~ /Output file\s*=\s*(\S+)\n/ms;
	my ($blength_tree, $support_tree) = read_from_file("$out_file.cnt") =~ /^(.*)\n(.*)$/m;
	return ($support_tree, $blength_tree);
}

################################################################################

################################################################################
# CALCULATE AIC
################################################################################
sub aic {
	my ($lnL, $parameters) = @_;
	my $aic = (-2*$lnL)+(2*$parameters);
	return $aic;
}

sub aicc {
	my ($lnL, $k, $n) = @_;
	my $aicc = aic($lnL, $k) + (2*$k*($k+1))/($n - $k - 1);
	return $aicc;
}

sub delta_aic {
	my ($aic_list) = @_;
	my $min_aic = min(@$aic_list);
	my @aic_delta;
	foreach my $aic (@$aic_list) {
		push @aic_delta, $aic-$min_aic;
	}
	return wantarray ? @aic_delta : \@aic_delta;
}
################################################################################

################################################################################
# RANDOM NUMBERS
################################################################################
# # Generate an 8-digit random number.
sub random_seed {
	my $seed = int(rand(100000000));
	until (length($seed) == 8) {
		$seed = int(rand(100000000));
	}
	return $seed;
}

# # The following is a nice way to get a truly random 10-digit integer,
# # but it won't work on Windows (without Cygwin).
# sub random_seed {
	# my ($source) = @_;
	# $source = "/dev/random" unless $source;
	# my $seed = "";
	# # Sometimes the '10' digit number is shorter, because there are
	# # effectively leading zeroes; we ignore those until we get 10-digits,
	# # so that the script returns numbers of consistent length.
	# until (length($seed) == 10) {
		# # Random bytes are returned, but we want integers, so
		# # convert them from octal with od.
		# $seed = `od -An -N4 -tu4 $source`;
		# $seed =~ s/\s//g;
	# }
	# return $seed;
# }

sub timestamp_id {
	my ($max_random) = @_;
	$max_random = 1000 unless $max_random;
	
	my $timestamp = strftime("%m%d%H%M%S", localtime);
	return "$timestamp\_".int(rand($max_random));
}
################################################################################

################################################################################
# FILE IO
################################################################################
sub read_from_file {
	my ($file) = @_;
	my $data;
	open(FILE, $file) || croak "Cannot open file '$file'";
	while (my $line = <FILE>) {
		$data .= $line;
	}
	close(FILE);
	return $data;
}

sub read_delim {
	my ($file, $cols, $col_name, $delim) = @_;
	my (%cols, @data, %data);
	$cols{$$cols[$_]} = $_ for (0..(scalar(@$cols)-1));
	$col_name = $$cols[0] unless $col_name;
	$delim = "\t" unless $delim;
	
	open(FILE, $file) || croak "Cannot open file '$file'";
	my $header = <FILE>;
	while (my $line = <FILE>) {
		chomp($line);
		my @row = split(/$delim/, $line);
		next unless @row;
		push @data, \@row;
	}
	close(FILE);
	
	foreach my $row (@data) {
		my $row_name = $$row[$cols{$col_name}];
		foreach my $col (keys %cols) {
			next if $col eq $row_name;
			croak "Non-unique row name '$row_name'" if exists($data{$row_name}{$col});
			if (defined $$row[$cols{$col}]) {
				$data{$row_name}{$col} = $$row[$cols{$col}];
			} else {
				$data{$row_name}{$col} = '';
			}
		}
	}
	return wantarray ? %data : \%data;
}

sub save_to_file {
	my ($data, $file, $append) = @_;
	$append = $append ? ">" : "";
	open(FILE, ">$append$file") || croak "Cannot open file '$file'";
	print FILE $data;
	close(FILE);
	return;
}
################################################################################

################################################################################
# NUMERIC CALCULATIONS
################################################################################
sub min {
	my @values = @_;
	my @sorted = sort {$a <=> $b} @values;
	return $sorted[0];
}

sub max {
	my @values = @_;
	my @sorted = sort {$b <=> $a} @values;
	return $sorted[0];
}

sub sum {
	my @values = @_;
	my $sum;
	$sum += $_ for @values;
	return $sum;
}
################################################################################

################################################################################
# SCRIPT EXECUTION
################################################################################
sub execute_check {
	my ($execute, $check, $check_success, $program, $fatal) = @_;
	$check_success = 0 unless $check_success;
	($program) = $execute =~ /^(\S)/ unless $program;
	$fatal = 0 unless $fatal;
	
	my $out = `$execute 2>&1`;
	if (!defined($out)) {
		if ($fatal) {
			croak "Aborting due to $program errors: failed to execute '$execute'";
		} else {
			carp "Continuing execution despite $program errors: failed to execute '$execute'";
		}
	}
	my $error_condition;
	if ($check_success) {
		$error_condition = $out !~ /$check/gms;
	} else {
		$error_condition = $out =~ /$check/gms;
	}
	if ($error_condition) {
		print STDERR "\n$program Output - Start\n";
		print STDERR "$out\n";
		print STDERR "$program Output - End\n\n";
		if ($fatal) {
			croak "Aborting due to $program errors: tried to execute '$execute'";
		} else {
			carp "Continuing execution despite $program errors: tried to execute '$execute'";
		}
	}
	return $out;
}
################################################################################

1;


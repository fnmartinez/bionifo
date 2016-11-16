#! /usr/bin/perl
use strict;
use warnings;

use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Basename;

use Bio::SeqIO;
use Bio::Seq;
use Bio::Perl;

my $input_file = "";
my $output_dir = "";
my $man = 0;
my $help = 0;

GetOptions("help|?" => \$help, 
           "man" => \$man,
           "input-file=s" => \$input_file,
           "output-directory=s" => \$output_dir)
or pod2usage(2);
pod2usage(1) if $help || !$input_file;
pod2usage(-exitval => 0, -verbose => 2) if $man;

$input_file = File::Spec->rel2abs($input_file);

if (!(-e $input_file || -r $input_file)) {
	print "Input file ${input_file} does not exists or is not readable.", "\n";
	exit 1;
}

if (!$output_dir) {
	$output_dir = cwd();
} elsif (!(-e $output_dir || -r $output_dir || -d $output_dir)) {
	print "Output directory ${output_dir} is not a directory, or does not exists, or is not readable", "\n";
	exit 1;
}

$output_dir = File::Spec->rel2abs($output_dir);

if (!(-w $output_dir)) {
	print "Output director ${output_dir} is not writible.", "\n";
	exit 1;
}

my @directions = ("fwd", "rev");
my $inseq = Bio::SeqIO->new(-file => $input_file, -format => 'genbank');
my $input_basename = fileparse($input_file);
(my $input_basename_extensionless = $input_basename) =~ s/\.[^.]+$//;

my $seq = $inseq->next_seq();
my $max_seq_length = 0;
my $most_probable_rf = shift;
my $most_probable_rf_dir = "";
my $most_probable_rf_index = 0;
for my $direction (@directions) {
	for (my $i = 0; $i < 3; $i++) {
		my $read_seq_no = $i+1;
		my $output_basename = "${input_basename_extensionless}-${direction}-${read_seq_no}.fas";
		my $output_file = File::Spec->catfile($output_dir, $output_basename);
		my $outseq = Bio::SeqIO->new(-file => ">" . $output_file, -format => 'fasta');
		$outseq->write_seq($seq->translate(undef, undef, $i));
		if ($seq->length > $max_seq_length) {
			$max_seq_length = $seq->length;
			$most_probable_rf = $seq->translate(undef, undef, $i);
			$most_probable_rf_dir = $direction;
			$most_probable_rf_index = $read_seq_no;
		}
	}
	$seq = $seq->revcom();
}

my $output_basename = "${max_seq_length}_${input_basename_extensionless}-${most_probable_rf_dir}-${most_probable_rf_index}.fas";
my $output_file = File::Spec->catfile($output_dir, $output_basename);
my $outseq = Bio::SeqIO->new(-file => ">" . $output_file, -format => 'fasta');
$outseq->write_seq($most_probable_rf);

__END__

=head1 NAME

ex1.pm - Sequence processing

=head1 SYNOPSIS

ex1.pm --input-file input_file.gb [--output-directory output_dir]

  Options:
    --help               Brief help message
    --man                Full documentation
    --input-file         Path where the input file is in Genbank format
    --output-directory   Optional. Current dir by default. Directory where the sequences should be placed.
    
=head1 OPTIONS

=over 8

=item B<-help>

Prints a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<--input-file>

Path to the input file in Genbank format containing a gene sequence.

=item B<--output-directory>

I<Optional>. The path with a directory where the reading sequence will be left. 
B<By default is the current working directory.> 

=back

=head1 DESCRIPTION

B<This program> will read the provided path in I<input_file>, fetching the file in Genbank format with a sequence of a gene and output the six reading frames into the output directory I<output_dir> in FASTA format.

=cut 
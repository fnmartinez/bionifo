#! /usr/bin/perl
# As seen in http://www.bioinfopoint.com/index.php/code/3-multiple-sequence-alignment-with-bioperl-and-muscle
use strict;
use warnings;

use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

use Bio::Tools::Run::Alignment::Muscle;
use Bio::AlignIO;
use Bio::SeqIO;

my $help = 0;
my $man = 0;
my $input_file = "";
my $output_basename = "blast.out";
my $output_dir = cwd();

# The parameters to be passed to MUSCLE
# Here, around 4Gb can be used for the alignment,
# and the process will have at most 100 interactions (if it does not converge before)
my @params = (quiet => 0, maxmb => '4000', maxiters => '100');

GetOptions("help|?" => \$help,
           "man" => \$man,
           "input-file=s" => \$input_file,
           "output-file=s" => \$output_basename
) or pod2usage(2);
pod2usage(1) if $help || !$input_file;
pod2usage(-exitval => 0, -verbose => 2) if $man;

$input_file = File::Spec->rel2abs($input_file);

if (!(-e $input_file || -r $input_file)) {
	print "Input file ${input_file} does not exists or is not readable.", "\n";
	exit 1;
}

my $output_file = File::Spec->catfile($output_dir, $output_basename);

if (-e $output_file) {
	print "Otput file ${output_file} already exists.", "\n";
	exit 1;
}

# The "factory" created for the alignment with the parameters above 
my $factory = Bio::Tools::Run::Alignment::Muscle->new(@params);

# The input file with protein sequences in FASTA format
my $str = Bio::SeqIO->new(-file=> $input_file, '-format' => 'fasta');

# The output file where the alignment will be writen - in clustalw format
my $out = Bio::AlignIO->new(-file   => ">".$output_file, -format => 'clustalw');


# First creates an array with all sequences to be alignment
my @protseq_array =();

while ( my $seq = $str->next_seq() ) {
	push (@protseq_array, $seq);

	print "Added sequence ".$seq->id." to the array\n";
}

# Then align the sequences in the array using the factory created before
my $aln = $factory->align(\@protseq_array);

$out->write_aln($aln);


__END__

=head1 NAME

ex2b.pl - MSA using MUSCLE algorithm

=head1 SYNOPSIS

ex2.pl --input-file input_file.gb [--output-file blast.out]

  Options:
    --help               Brief help message
    --man                Full documentation
    --verbose            Add verbosity to the script
    --input-file         Path where the input file is in Genbank format
    --output-file        Optional. I<blast.out> by default. Filename .
    
=head1 OPTIONS

=over 8

=item B<-help>

Prints a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<--input-file>

Path to the input file in Genbank format containing a gene sequence.

=item B<--output-file>

I<Optional>. The filename where the BLAST result will be put. 
B<By default is I<blast.out>.> 

=back

=head1 DESCRIPTION

B<This program> will read the provided path in I<input_file>, fetching the file in FASTA format with a sequences. The result of applying the MUSCLE algorithm is left in I<output_file>.

=cut 
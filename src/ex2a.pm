#! /usr/bin/perl
# As seen in http://search.cpan.org/dist/BioPerl/Bio/Tools/Run/RemoteBlast.pm
use strict;
use warnings;

use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

use Bio::Tools::Run::RemoteBlast;
use Bio::SeqIO;

my $help = 0;
my $man = 0;
my $v = 0;
my $input_file = "";
my $output_basename = "blast.out";
my $output_dir = cwd();

my @params = ( '-prog' => 'blastp',
		       '-data' => 'swissprot',
		       '-expect' => '1e-10',
		       '-readmethod' => 'SearchIO' );

GetOptions("help|?" => \$help,
           "man" => \$man,
           "verbose" => \$v,
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

my $factory = Bio::Tools::Run::RemoteBlast->new(@params);

my $r = $factory->submit_blast($input_file);

print "waiting...\n" if( $v > 0 );
while ( my @rids = $factory->each_rid ) {
	foreach my $rid ( @rids ) {
		my $rc = $factory->retrieve_blast($rid);
		if( !ref($rc) ) {
			if( $rc < 0 ) {
				$factory->remove_rid($rid);
			}
			print "." if ( $v > 0 );
			sleep 5;
		} else {
			my $result = $rc->next_result();
			#save the output
			$factory->save_output($output_file);
			$factory->remove_rid($rid);
			while ( my $hit = $result->next_hit ) {
				next unless ( $v > 0);
				print "\thit name is ", $hit->name, "\n";
				while( my $hsp = $hit->next_hsp ) {
					print "\t\tscore is ", $hsp->score, "\n";
				}
			}
		}
	}
}


__END__

=head1 NAME

ex2a.pl - Remote BLAST over a FASTA sequence

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

=item B<--verbose>

Prints messages about the execution that is being done. Useful since the remote execution of BLAST can take a while.

=item B<--input-file>

Path to the input file in Genbank format containing a gene sequence.

=item B<--output-file>

I<Optional>. The filename where the BLAST result will be put. 
B<By default is I<blast.out>.> 

=back

=head1 DESCRIPTION

B<This program> will read the provided path in I<input_file>, fetching the file in FASTA format with a sequence of aminoacids, result of the translation of a certain gene. The result of BLAST is left into I<output_file>.

=cut 
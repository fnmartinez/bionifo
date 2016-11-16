#! /usr/bin/perl
# As seen in http://search.cpan.org/dist/BioPerl/Bio/Tools/Run/RemoteBlast.pm
use strict;
use warnings;

use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;

use Bio::Tools::Run::RemoteBlast;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::SeqIO;

my $help = 0;
my $man = 0;
my $verbose = 0;
my $remote = 0;
my $input_file = "";
my $blast_db_name = "swissprot";
my $blast_bin_dir = "/usr/bin";
my $blast_evalue = 10;
my $output_basename = "blast.out";
my $output_dir = cwd();

my @params = ( '-prog' => 'blastp',
		       '-data' => $blast_db_name,
		       '-expect' => "1e-${blast_evalue}",
		       '-readmethod' => 'SearchIO' );

GetOptions("help|?" => \$help,
           "man" => \$man,
           "verbose" => \$verbose,
           "remote" => \$remote,
           "db_name" => \$blast_db_name,
           "prog_dir" => \$blast_bin_dir,
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
	print "Output file ${output_file} already exists.", "\n";
	exit 1;
}

if ($remote) {
	my $factory = Bio::Tools::Run::RemoteBlast->new(@params);
	
	my $r = $factory->submit_blast($input_file);
	
	print "waiting...\n" if( $verbose > 0 );
	while ( my @rids = $factory->each_rid ) {
		foreach my $rid ( @rids ) {
			my $rc = $factory->retrieve_blast($rid);
			if( !ref($rc) ) {
				if( $rc < 0 ) {
					$factory->remove_rid($rid);
				}
				print "." if ( $verbose > 0 );
				sleep 5;
			} else {
				my $result = $rc->next_result();
				#save the output
				$factory->save_output($output_file);
				$factory->remove_rid($rid);
				while ( my $hit = $result->next_hit ) {
					next unless ( $verbose > 0);
					print "\thit name is ", $hit->name, "\n";
					while( my $hsp = $hit->next_hsp ) {
						print "\t\tscore is ", $hsp->score, "\n";
					}
				}
			}
		}
	}	
} else {
	$blast_db_name = File::Spec->rel2abs($blast_db_name);

	if (!(-e $blast_db_name || -r $blast_db_name)) {
		print "Blast DB file ${blast_db_name} does not exists or is not readable.", "\n";
		exit 1;
	}
	
	$blast_bin_dir = File::Spec->rel2abs($blast_bin_dir);
	
	if (!(-e $blast_bin_dir || -r $blast_bin_dir || -d $blast_bin_dir)) {
		print "Blast binary directory ${blast_bin_dir} is not a directory, or does not exists, or is not readable", "\n";
		exit 1;
	}
	
	$blast = StandAloneBlastPlus->new(
		-db_name => $blast_db_name,
		-prog_dir => $blast_bin_dir
	);
	$result = $blast->blastp(
		-query => $input_file,
		-outfile => $output_file
		-method_args => [ '-evalue' => $blast_evalue ]
	);
}


__END__

=head1 NAME

ex2a.pm - BLAST over a FASTA sequence of proteins

=head1 SYNOPSIS

ex2a.pm --input-file input_file.gb [--output-file blast.out]

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

=item B<-remote>

Uses the remote BLAST+ algorithm instead of the local one

=item B<-verbose>

Prints messages about the execution that is being done. Useful since the remote execution of BLAST can take a while.

=item B<-input-file>

Path to the input file in Genbank format containing a gene sequence.

=item B<-output-file>

I<Optional>. The filename where the BLAST result will be put. 
B<By default is I<blast.out>.> 

=item B<db_name>

I<Only if using local>. The database name or location to be used by BLAST+
B<By default is I<swissprot>>


=item B<prog_dir>

I<Only if using local>. The directory where the BLAST+ program resides.
B<By default is I</usr/bin>>

=back

=head1 DESCRIPTION

B<This program> will read the provided path in I<input_file>, fetching the file in FASTA format with a sequence of aminoacids in a protein. Then it will apply the BLAST algorithm and leave the output in the I<output_file>.

=cut 
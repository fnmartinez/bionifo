#! /usr/bin/perl
use strict;
use warnings;

use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Basename;

use Bio::Perl;
use Bio::Seq; 
use Bio::SeqIO; 
use Bio::SearchIO;
use Bio::DB::GenBank;
use Data::Dumper; 

my $input_file = "";
my $pattern = "";
my $help = 0;
my $man = 0;
my $output_file = "output.fas"; 
my $output_format = 'fasta';

GetOptions("help|?" => \$help, 
           "man" => \$man,
           "input-file=s" => \$input_file,
           "pattern=s" => \$pattern,
           "output-file" => \$output_file)
or pod2usage(2);
pod2usage(1) if $help || !$input_file;
pod2usage(-exitval => 0, -verbose => 2) if $man;

$input_file = File::Spec->rel2abs($input_file);

if (!(-e $input_file || -r $input_file)) {
	print "Input file ${input_file} does not exists or is not readable.", "\n";
	exit 1;
}

if (!$pattern) {
	print "Pattern empty", "\n";
	exit 1;
}

$output_file = File::Spec->rel2abs($output_file);

if (-e $output_file) {
	print "Output file ${output_file} already exists.", "\n";
	exit 1;
}

my $searchio = new Bio::SearchIO (-format => 'blast',
                                  -file   => $input_file);
my $seqOut = Bio::SeqIO->new(-file => '>' . $output_file,
			                 -format => $output_format);

while (my $result = $searchio->next_result()) {
	$result->database_name;
	
	while (my $hit = $result->next_hit) {
		
		my $acc = $hit->accession();
		
		my $str = $hit->description();

		if(index($str, $pattern) != -1){
			
			print "HIT \n";
			while ( (my $khit, my $vhit) = each %{$hit}) { 
		 		print "$khit => $vhit \n"; 
			}
			print "HSP \n";
			while(my $hsp = $hit->next_hsp()) {
				while ( (my $khsp,my $vhsp) = each %{$hsp}) {
				    if ($khsp eq "EVALUE") {   
     					print "$khsp => $vhsp \n";
     				}  
     			}
			}
			print "\n";

			my $gb = Bio::DB::GenBank->new(-retrievaltype => 'tempfile' , -format => $output_format);
			my $seqIn = $gb->get_Stream_by_acc($acc);
			# write each entry in the input file to the output file 
			while (my $inseq = $seqIn->next_seq) {
        			$seqOut->write_seq($inseq);
			}

		}
	}
}

__END__
=head1 NAME

ex3.pm - Analyze a BLAST output by searching a pattern

=head1 SYNOPSIS

ex3.pm --input-file blast.out --pattern Arabidopsis [--output-file output.fas]

=head1 OPTIONS

=over 8

=item B<-help>

Prints a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=item B<-input-file>

Path to the input file containing the BLAST output.

=item B<-pattern>

Pattern to search in the blast output.

=item B<-output-file>

I<Optional>. The filename where the BLAST result will be put. 
B<By default is I<blast.out>.> 

==head1 DESCRIPTION

B<This program> will read the provided path in the I<input_file> with the output of a blast algorithm and a I<pattern>. It will try to find matches of the I<pattern> in the I<input_file> by parsing it. If there is a match, the sequence of this match will be put in the I<output_file> in FASTA format.  

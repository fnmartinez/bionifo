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


my $help = 0;
my $man = 0;
my $input_file = "";
my $output_basename = "blast_parsed.fas";
my $output_dir = cwd();
my $pattern = "";

GetOptions("help|?" => \$help,
           "man" => \$man,
           "input-file=s" => \$input_file,
           "output-file=s" => \$output_basename,
           "pattern=s" => \$pattern
) or pod2usage(2);
pod2usage(1) if $help || !$input_file || !$pattern;
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

# Get report
my $searchio = new Bio::SearchIO (-format => 'blast',
                                  -file   => $input_file);
my $output_file_format = 'fasta';
my $seq_out = Bio::SeqIO->new(-file => '>' . $output_file,
			      -format => $output_file_format);

while (my $result = $searchio->next_result()) {
	# Get report info
	$result->database_name;
	
	#Bio::Search::HitI object
	while (my $hit = $result->next_hit) {
		
		my $acc = $hit->accession();
		
		my $str = $hit->description();

		#The index function searches for one string within another
		#If the substring is not found, index returns -1.
		if(index($str, $pattern) != -1){
			
			print "HIT", "\n";
			while ( (my $khit, my $vhit) = each %{$hit}) { 
		 		print "${khit} => ${vhit}", "\n"; 
			}
			print "HSP", "\n";
			while(my $hsp = $hit->next_hsp()) {
				while ( (my $khsp,my $vhsp) = each %{$hsp}) {
				    if ($khsp eq "EVALUE") {   
     					print "${khsp} => ${vhsp}", "\n";
     				}  
     			}
			}
			print "\n";

			my $gb = Bio::DB::GenBank->new(-retrievaltype => 'tempfile' , -format => 'fasta');
			my $seqIn = $gb->get_Stream_by_acc($acc);
			# write each entry in the input file to the output file 
			while (my $inseq = $seqIn->next_seq) {
        			$seq_out->write_seq($inseq);
			}

		}
	}
}

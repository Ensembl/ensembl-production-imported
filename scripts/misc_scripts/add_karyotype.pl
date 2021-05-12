
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use Getopt::Long;
use IO::File;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::KaryotypeBand;
use Carp;

# populate core karyotype table with cytogenetic band data

my $opts ={}; 

GetOptions( $opts, '-registry=s', '-species=s', '-in=s', '-h' , '-v' );

&check_opts($opts);

my $reg = 'Bio::EnsEMBL::Registry';

#set up the registry and adaptors

$reg->load_all( $opts->{'registry'} );

my $ka = $reg->get_adaptor( $opts->{'species'} ,"core", "KaryotypeBand");  
my $sa = $reg->get_adaptor( $opts->{'species'} ,"core", "Slice");

my $fh = IO::File->new( $opts->{'in'} ) or croak "unable to open file $opts->{'in'}\n";

my $header = $fh->getline();

unless( $header =~ /chromosome	band	start	end	stain/i ){
    croak "header line for file $opts->{'in'} is not valid\n";
}


my $chroms = {};
my $existing = 0;
my $new = 0;


while( my $l = $fh->getline ){

    chomp( $l );

    my $d = [ split(/\t/,$l)];

    my $chrom = $$d[0];
    my $band  = $$d[1];
    my $start = $$d[2];
    my $end   = $$d[3];
    my $stain = $$d[4];

    my $slice = $sa->fetch_by_region( 'toplevel' , $chrom, $start, $end ) or croak "unable to obtain slice for $chrom $start $end\n";

    my $b = $ka->fetch_all_by_chr_band( $chrom , $band );

    if ( scalar @$b > 0 ){ ++$existing }
    else{

#N.B. the start end coordinates are relative to the slice you supply!
	++$new;
	my $obj = Bio::EnsEMBL::KaryotypeBand->new( -START => 1, -END => (( $end - $start ) + 1), -NAME => $band, -STAIN => $stain, -SLICE => $slice ) or 
	    croak "unable to generate karyotype object for band $band\n";


	$ka->store( $obj );

    }

    ++$chroms->{$chrom};

}

print "ignored $existing existing bands\n";
print "processed $new new bands\n";

for my $k ( sort keys %$chroms){
    print "$k\t$chroms->{$k}\n";
}

print "\n*** Remember to check that you have populated the seq_region_attrib table ***\n";
print "*** with the ordered chromosomes using the attrib_type = karyotype_rank value ***\n\n";

exit 0;

######################################################################################

sub check_opts{

   my $options = shift @_;

   if ( $options->{'h'} or $opts->{'help'} ){ &usage }

   for my $f (  'registry', 'species', 'in' ){

       unless ( $options->{$f} ) { &usage("$f parameter not supplied" ) }
   }

   unless( -f $options->{'registry'} ){ &usage( "-registry $options->{'registry'} must be a file\n" ) }
}

sub usage{

   my $message = shift @_;

   if ( $message ){ print">>>> $message\n\n" }

   print<<"EOF";

usage: uniprot2seq.pl -registry <file> -species <string> -in <file> -v

-registry      Ensmebl registry file
-species       valid species name from registry file
-in            file containing karyotye data
-v             verbose reporting

Populate an Ensembl core database karotype + seq_region_attrib table with karyotype data.
Karyotype data is supplied with the the -in file of the format

chromosome	band	start	end	stain
X	1A	1	3063322	TEL
X	1B	3063323	4968023
X	2A	4968024	5992912

All chromosome values must exist in the corresponding core database, and stain values are used
to indicate the centromere/telomere regions as follows

TEL  = telomere
ACEN = centromer

Data is written directly to the karyotpe table and the seq_region_attrib table is populated for
the corressponding chromosome sequences.

EOF

   if ( $message ){ exit 1 }
   else{ exit 0 }
}


=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

=head1 NAME


=head1 DESCRIPTION


=head1 MAINTAINER

$Author: ckong $

=cut
use warnings;
use strict;
use Bio::EnsEMBL::Utils::CliHelper;
use Log::Log4perl qw/:easy/;
use Data::Dumper;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Analysis;

Log::Log4perl->easy_init($INFO);
my $logger     = Log::Log4perl->get_logger();

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();
my $optsd      = $cli_helper->get_dba_opts();
push(@{$optsd},"file:s");
my $opts       = $cli_helper->process_args($optsd, \&usage);

$logger->info('Creating db adaptors...\n');
# use the command line options to get an array of database details
my $dba;

for my $db_args (@{$cli_helper->get_dba_args_for_opts($opts)}) {
   $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%{$db_args});
}

#for my $db_args ( @{ $cli_helper->get_dba_args_for_opts( $opts, 0 ) } ) {

#        my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$db_args);
#}

#my $db_args = $cli_helper->get_dba_args_for_opts( $srcopts, 0 );

#die "GAF file is required !\n" unless @ARGV;

my $analysis_obj = Bio::EnsEMBL::Analysis->
    new( -logic_name      => 'gaf_loader',
         -program         => 'gaf_loader.pl',
         -description     => 'GAF loader',
         -display_label   => 'GAF loader',
    );


$logger->info('Getting adaptors...\n');
my $ga    = $dba->get_GeneAdaptor() || die "problem getting gene adaptor\n";

my %count;
my $file = $opts->{file};
#my $file = $ARGV[0];
open(FILE, $file) || die "problem opening GAF file";

$logger->info('Start parsing GAF file...\n');
while (<FILE>) {
       my $line = $_;
       
       next if $line =~/^!/;
       chomp ($line);
       my @cols      = split(/\t/, $_, 17);

       die "Is the format GAF2.0?\n" unless @cols == 17;

       my ($db, $db_objID, $db_objSym, $qual, $go_id,
           $db_ref, $evidence_code, $with_or_froms, $aspect,
           $db_objName, $db_objSyn, $db_objTyp,
           $taxons, $date, $assigned_by, $annot_ext,
           $gene_product_form_id
          ) = @cols;

=pod

UniProtKB	A2P2R3	YMR084W		GO:0003674	GO_REF:0000015	ND		F	Putative glutamine--fructose-6-phosphate aminotransfer
ase [isomerizing]	YM084_YEAST|YMR084W	protein	taxon:559292	20030203	SGD		
UniProtKB	A2P2R3	YMR084W		GO:0004360	GO_REF:0000003	IEA	EC:2.6.1.16	F	Putative glutamine--fructose-6-phosphate amino
transferase [isomerizing]	YM084_YEAST|YMR084W	protein	taxon:559292	20140118	UniProt		
UniProtKB	A2P2R3	YMR084W		GO:0005575	GO_REF:0000015	ND		C	Putative glutamine--fructose-6-phosphate aminotransfer
ase [isomerizing]	YM084_YEAST|YMR084W	protein	taxon:559292	20030203	SGD		
UniProtKB	A2P2R3	YMR084W		GO:0006048	GO_REF:0000041	IEA	UniPathway:UPA00113	P	Putative glutamine--fructose-6-phospha
te aminotransferase [isomerizing]	YM084_YEAST|YMR084W	protein	taxon:559292	20140118	UniProt	

=cut
      # Check for mandatory columns 
      die "'",
      join("'\t'",
           $db, $db_objID, $db_objSym, $go_id,
           $db_ref, $evidence_code, $aspect, $db_objTyp,
           $taxons, $date, $assigned_by,
          ), "'\n"
      unless
           $db && $db_objID && $db_objSym && $go_id &&
           $db_ref && $evidence_code && $aspect && $db_objTyp && 
           $taxons && $date && $assigned_by;

      next if
        ( $go_id eq 'GO:0005575' ||   # Cellular component
          $go_id eq 'GO:0003674' ||   # Molecular function
          $go_id eq 'GO:0008150' ) && # Biological process
          $db_ref eq 'TAIR:Communication:1345790' && $evidence_code eq 'ND';

      my @qual               = split(/\|/, $qual);
      my @db_ref             = split(/\|/, $db_ref);
      my @with_or_froms      = split(/\|/, $with_or_froms);
      my @db_objSyn          = split(/\|/, $db_objSyn);
      my ($taxon, $pathogen) = split(/\|/, $taxons, 2);
      my @annot_ext          = split(/\|/, $annot_ext);

      my $gene_obj;

      # use the db_obj_name as a stable_id...
      if ($db_objName){
	   $gene_obj = $ga->fetch_by_stable_id($db_objName);
      }

      # else check if there is *one* xref for the db_objSym
      unless ($gene_obj){
        my $gene_obj_aref = $ga->fetch_all_by_external_name($db_objSym );
        if (@$gene_obj_aref == 1){
            $gene_obj = $gene_obj_aref->[0];
        }
      }

      # else use synonyms as xrefs
      unless ($gene_obj){
        if(@db_objSyn){
            my %best_guess;
            # Can have > 1 synonyms
            for my $db_object_synonym (@db_objSyn){
                my $gene_obj_aref = $ga->fetch_all_by_external_name( $db_objSyn );
                # each synonym match > 1 gene
                for my $gene_obj (@$gene_obj_aref){
		    $best_guess{$gene_obj->stable_id}++;
                }
            }
 	    my ($best_guess) = (sort {$best_guess{$b} <=> $best_guess{$a}} keys %best_guess);
            $gene_obj =  $ga->fetch_by_stable_id( $best_guess );
	}
     }

     unless ($gene_obj){
	$logger->info('Can\'t match this record\n$line\n\n');
        $count{'not found'}++;
        next;
     }

     my $string = join("\t",$gene_obj->stable_id, $go_id), "\n";
     $logger->info($string);
     #print join("\t",$gene_obj->stable_id, $go_id), "\n";

     # Annotate genes 'canonical translation'
     my $ensembl_obj      = $gene_obj->canonical_transcript->translation;
     my $ensembl_obj_type = 'Translation';

     # Annotate genes 'transcript' for genes (such as snoRNAs or annotated pseudogenes) that don't have a translation
     unless($ensembl_obj){
        $ensembl_obj      = $gene_obj->canonical_transcript;
        $ensembl_obj_type = 'Transcript';
    }

    unless($ensembl_obj){
	$logger->info('line $. : No ensembl objects can be match and annotate on!\n');
        next;
    }

    my $dbEntry = Bio::EnsEMBL::DBEntry->
        new( -primary_id   => $go_id,
             -display_id   => $go_id,
             -dbname       => 'GO',
             -linkage_type => $evidence_code,
             -analysis     => $analysis_obj,
             -linkage_annotation => $ARGV,
        );

    $ensembl_obj->add_DBEntry( $dbEntry );

    $dba->get_DBEntryAdaptor->store( $dbEntry, $ensembl_obj->dbID, $ensembl_obj_type, 1 );
}




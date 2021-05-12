#!/usr/bin/env perl
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Utils::CliHelper;
use Carp;
use Data::Dumper;
use Log::Log4perl qw(:easy);
use JSON;
use HTTP::Tiny;

my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

# get the basic options for connecting to a database server
my $optsd = [ @{ $cli_helper->get_dba_opts() } ];
push( @{$optsd}, "verbose" );
push( @{$optsd}, "url:s" );

# process the command line with the supplied options plus a help subroutine
my $opts = $cli_helper->process_args( $optsd, \&pod2usage );
if ( $opts->{verbose} ) {
  Log::Log4perl->easy_init($DEBUG);
}
else {
  Log::Log4perl->easy_init($INFO);
}

my $logger = get_logger();

$logger->info( "Loading " . $opts->{dbname} );
my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -USER   => $opts->{user},
                                     -PASS   => $opts->{pass},
                                     -HOST   => $opts->{host},
                                     -PORT   => $opts->{port},
                                     -DBNAME => $opts->{dbname} );


my $http = HTTP::Tiny->new();
 
$opts->{url} ||= 'https://plasmogem.sanger.ac.uk/genes/list';
$logger->info("Retrieving genes from ".$opts->{url});
my $response = $http->get($opts->{url}, {
  headers => { 'Content-type' => 'application/json' }
});

if(!length $response->{content}) {
    croak "No content from ".$opts->{url};
}

my %plasmogen_hash = map {$_->{gene_id}=>$_} grep {$_->{transfection_resource_count_ko} gt 0|| $_->{transfection_resource_count_ko} gt 0} @{decode_json($response->{content})->{json_data}};

my $dbentry_adaptor = $dba->get_DBEntryAdaptor();
$logger->info("Removing existing PlasmoGem cross-references");
$dba->dbc()->sql_helper()->execute_update(
    -SQL=>q/delete ox.*,x.* from object_xref ox 
    join xref x using (xref_id)
    join external_db e using (external_db_id)
    where e.db_name='plasmogem'/
    );
for my $gene (@{$dba->get_GeneAdaptor()->fetch_all()}) {
    my $plasmo = $plasmogen_hash{$gene->stable_id()};
    if(defined $plasmo) {
        my @res = ();
        if($plasmo->{library_clone_count} gt 0) {
            push @res, "library clone";
        }

        if($plasmo->{transfection_resource_count_ko} gt 0) {
            push @res, "knock-out vector";
        }

        if($plasmo->{transfection_resource_count_tag} gt 0) {
            push @res, "tagging vector";
        }

        my $desc = "PlasmoGEM resources available: ".join(', ',@res);
        $logger->info("Adding PlasmoGem cross-reference to ".$gene->stable_id().": ".$desc);
        $dbentry_adaptor->store( Bio::EnsEMBL::DBEntry->new(
                                     -PRIMARY_ID=>$gene->stable_id(),
                                     -DISPLAY_ID=>$gene->stable_id(),
                                     -DBNAME=>'plasmogem',
                                     -DESCRIPTION=>$desc,
                                     -INFO_TYPE  => 'DIRECT'
                                 ), $gene->dbID(),
                                 'Gene' )
    }
}

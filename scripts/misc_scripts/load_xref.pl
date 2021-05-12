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

=pod

=head1 NAME

  load_xrefs.pl

=head1 SYNOPSIS

  Load custom Xrefs into coredb based on object stable_ids.

=head1 DESCRIPTION

  Inserts custom Xrefs for objects with costructable adaptors (Gene, Transctript, etc.).
  Specify object type with `-object` parameter and xref tag (aka 'external_db') with `-xref_name`.
  (For a list of external_dbs see: 'select * from external_db' )
  Provide list of `stable_id \t external_value` (tab-separated) as input (/dev/stdin). 

=head1 ARGUMENTS

  perl load_xrefs.pl
         -dbname
         -host
         -port
         -user
         -pass
         -object
         -xref_name
         -help

=head1 EXAMPLE

  echo -e 'DTIU000000\tB4U79_03975' |
    perl ./load_xrefs.pl \
      -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> \
      -dbname <core_db> \
      -object 'Gene' -xref_name 'RefSeq_gene_name'

=cut

use warnings;
use strict;

use Getopt::Long;
use Pod::Usage qw(pod2usage);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($host, $port, $user, $pass, $dbname);
my ($object, $xref_name);
my $help = 0;

&GetOptions(
  'host=s'      => \$host,
  'port=s'      => \$port,
  'user=s'      => \$user,
  'pass=s'      => \$pass,
  'dbname=s'    => \$dbname,
  'object=s'    => \$object,
  'xref_name=s' => \$xref_name,
  'help|?'      => \$help,
) or pod2usage(-message => "use -help", -verbose => 1);

pod2usage(-verbose => 2) if $help;

my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $a = $core_db->get_adaptor($object);
my $dbea = $core_db->get_DBEntryAdaptor;

while (<STDIN>) {
  chomp;

  my ($stable_id, $ext_id, $ext_disp) = split/\t/;
  $stable_id =~ s/\..*$//;

  my $obj = $a->fetch_by_stable_id($stable_id);
  
  $ext_disp = $ext_id if (!$display_ne_primary_id || !defined $ext_id);

  my $xref_entry = store_ensembl_xref($dbea, $object, $obj->dbID, $xref_name, $ext_id, $ext_disp);

  if ($xref_entry && $update_display_xref) {
    eval { $obj->display_xref($xref_entry) };
    if ($@) {
      warn "Failed to update display_xref_id for $object $stable_id ($ext_id:$ext_disp)";
    } else {
      $a->update($obj);
    }
  }
}

sub store_ensembl_xref {
  my ( $dbea, $object_type, $id, $xref_name, $external_id, $external_display ) = @_;

  # make an xref
  my $entry = new Bio::EnsEMBL::DBEntry(
     -adaptor    => $dbea,
     -primary_id => $external_id,
     -display_id => $external_display,
     -version    => 0,
     -dbname     => $xref_name,
     -info_type   => $info_type,
     -info_text   => $info_text,
  );

  # store xref
  my $ignore_release = 1;
  return $entry if $dbea->store( $entry, $id, $object_type, $ignore_release);

  return;
}

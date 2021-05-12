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

  remove_entities.pl

=head1 SYNOPSIS

  Removes objects from core db based on objects types and stable_ids.

=head1 DESCRIPTION

  Removes objects from core db.
  Specify object type with `-object` parameter.
  (For a list of external_dbs see: 'select * from external_db' )
  Provide list of `stable_id` as input (/dev/stdin). 

=head1 ARGUMENTS

  perl remove_entities.pl
         -dbname
         -host
         -port
         -user
         -pass
         -object
         -help

=head1 EXAMPLE

  echo -e 'HMIM002325-PA' |
    perl ./remove_entities.pl \
      -host <db_host> -port <db_port> -user <db_user> -pass <db_pass> \
      -dbname <core_db> \
      -object 'Traslation'

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
  'help|?'      => \$help,
) or pod2usage(-message => "use -help", -verbose => 1);

pod2usage(-verbose => 2) if $help;

my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );
my $a = $core_db->get_adaptor($object);

while (<STDIN>) {
  chomp;

  my ($stable_id) = $_;
  $stable_id =~ s/\..*$//;

  my $obj = $a->fetch_by_stable_id($stable_id);
  next if !$obj;

  $a->remove($obj);
}

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

  get_trans.pl

=head1 SYNOPSIS

  get all posssible translation from core db

=head1 DESCRIPTION

  get all posssible translation from core db

=head1 ARGUMENTS

  perl get_trans.pl
         -host
         -port
         -user
         -pass
         -help
         -dbname
         -type transcript|translation|cds
         -ignore_biotypes transposable_element,...

=head1 EXAMPLE

  perl get_trans.pl $($CMD details script) -type translation -dbname anopheles_funestus_core_1906_95_3 > pep.fa

=cut

use warnings;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage qw(pod2usage);

my ($host, $port, $user, $pass, $dbname);
my ($type, $ignore_biotypes_str);


my $help = 0;

&GetOptions(
  'host=s'      => \$host,
  'port=s'      => \$port,
  'user=s'      => \$user,
  'pass=s'      => \$pass,
  'dbname=s'    => \$dbname,
  'type=s'    => \$type,
  'ignore_biotypes:s'    => \$ignore_biotypes_str,

  'help|?'      => \$help,
) or pod2usage(-message => "use -help", -verbose => 1);

if (!$type || $type !~ /^(transcript|translation|cds)$/i ) {
  pod2usage(-message => "no known type specified", -verbose => 1);
}

pod2usage(-verbose => 2) if $help;

my %ignore_biotypes = %{{ map {$_=>1} grep {!!$_} split /,/, ($ignore_biotypes_str or '') }};

my $core_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $want_translation = ($type =~ /translation/i);
my $want_cds = ($type =~ /cds/i);

my $ta = $core_db->get_adaptor("Transcript");
if ($ta) {
  my $pctrs = $ta->fetch_all(); 
  if ($pctrs) {
    while (my $tr = shift @{$pctrs}) {
      next if not ($tr);
      next if exists $ignore_biotypes{$tr->biotype};

      my $seq;
      if ($want_translation) {
        my $prot = $tr->translate(); 
        next if not ($prot);
        $seq =  $prot->seq();
      } elsif ($want_cds) {
        $seq = $tr->translateable_seq();
      } else {
        $seq = $tr->spliced_seq();
      }
      next if (!$seq);
      printf ">%s %s\n%s\n", $tr->stable_id(), $tr->dbID(), $seq;
    }
  }
}



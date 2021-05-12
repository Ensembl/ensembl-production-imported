#!/usr/bin/env/perl
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

# populate core karyotype table with data extracted from gff (for FlyBase/Dmel)

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use DBI;

my ($type, $regex_for_id, $host, $port, $user, $pass, $dbname);

GetOptions(
  "type=s", \$type,
  "regex_for_id=s", \$regex_for_id,
  "host=s", \$host,
  "P|port=i", \$port,
  "user=s", \$user,
  "p|pass=s", \$pass,
  "dbname=s", \$dbname,
);

$type = "chromosome_band" unless $type;
$regex_for_id = '([^;]+)' unless $regex_for_id;

undef $/;
my $gff = <>;
my @bands = $gff =~ /^(.+\t(?:$type)\t.+)$/gm;
my %bands;
foreach my $band (@bands) {
    if ($band =~ /ID=$regex_for_id;/) {
        my $id = $1;
        my ($chr, $start, $end) = $band =~ /^(\S+)\t[^\t]+\t[^\t]+\t(\d+)\t(\d+)/;
        $bands{$chr}{$id}{'start'} = $start;
        $bands{$chr}{$id}{'end'} = $end;
    }
}

# I'd expected to be able to do this by creating then storing a Bio::EnsEMBL
# KaryotypeBand object; but there is no method to store such an object...
my $db = DBI->connect("DBI:mysql:host=$host;port=$port;dbname=$dbname", $user, $pass);
my $chk_sth = $db->prepare("SELECT COUNT(*) FROM karyotype");
$chk_sth->execute();
my $count = ($chk_sth->fetchrow_array())[0];
if ($count > 0) {
  die "Karyotype table in $dbname has data ($count rows)";
}

my $sr_sql =
    'SELECT sr.seq_region_id '.
    'FROM seq_region sr '.
    'INNER JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id '.
    'WHERE sr.name=? AND cs.name = "chromosome";';
my $sr_sth = $db->prepare($sr_sql);

my $insert_sql =
    'INSERT INTO karyotype '.
        '(seq_region_id, seq_region_start, seq_region_end, band, stain) '.
    'VALUES (?,?,?,?,?)';
my $insert_sth = $db->prepare($insert_sql);

my $counter = 0;
my @stain = ('gpos25', 'gpos75');
foreach my $chr (sort keys %bands) {
    my $sr_id;
    $sr_sth->execute($chr);
    $sr_sth->bind_columns(\$sr_id);
    if ($sr_sth->fetch()) {
        foreach my $id (sort keys %{$bands{$chr}}) {
            $insert_sth->execute(
                $sr_id,
                $bands{$chr}{$id}{'start'},
                $bands{$chr}{$id}{'end'},
                $id,
                $stain[$counter++ % 2]
            ) || die "Error inserting into new karyotype table";
        }
    } else {
        die "Unrecognised chromosome $chr";
    }
}

print "Inserted ".($counter++)." bands.\n";

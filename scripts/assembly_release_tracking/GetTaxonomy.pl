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

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;


# Parameters:
my ($sp_taxon_tsv);
my ($ens_release);
my ($ncbi_taxon_host);
my ($ncbi_taxon_port);
my $group_ranks = 0;

my $usage="perl GetTaxonomy.pl -sp_taxon_tsv input.tsv -ens_release 113 -ncbi_taxon_host mysql-ens-xx-x -ncbi_taxon_port 1234 (OPTIONAL: -group_ranks)";

# Handle input arguments
GetOptions(
    "sp_taxon_tsv=s" => \$sp_taxon_tsv,
    "group_ranks" => \$group_ranks,
    "ens_release=s" => \$ens_release,
    "ncbi_taxon_host=s" => \$ncbi_taxon_host,
    "ncbi_taxon_port=s" => \$ncbi_taxon_port,
) or die("Error in command line arguments\n");

if ((!$sp_taxon_tsv) || (!$ncbi_taxon_host) || (!$ncbi_taxon_port) || (!$ens_release)) {
    print "$usage\n";
    die;
}

my $tax_dba =  Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
    -user   => "ensro",
    -dbname => "ncbi_taxonomy_${ens_release}",
    -host   => "${ncbi_taxon_host}",
    -port   => "${ncbi_taxon_port}");

my $node_adaptor = $tax_dba->get_TaxonomyNodeAdaptor();

# Declare the taxonomic ranks of interest (sorted)
my @target_ranks = qw(phylum class superorder order suborder family subfamily species);
# Number of ranks (plus the default addition of "GenomeDB name" as last rank)
my $num_ranks = $#target_ranks + 1;
# Create the taxonomy template: if a rank is missing, assign "N/A"
my $null_value = 'N/A';
my %rank_template = map {$_ => $null_value} @target_ranks;
# Get the taxon id for each species
my %taxon_hash;
open(my $tsv_handle, '<', $sp_taxon_tsv) or die "Cannot open $sp_taxon_tsv";
my $this_header = <$tsv_handle>;
my @head_cols = split(/[\f\n\r\t]+/, $this_header);
while ( my $line = <$tsv_handle> ) {
    my $row = map_row_to_header($line, \@head_cols);
    #print Dumper $row;
    $taxon_hash{$row->{"Name"}} = $node_adaptor->fetch_by_taxon_id($row->{"Taxon ID"});
}
close($tsv_handle);

# @taxonomy_table will have one string per species' taxonomy in tab-separated
# values (TSV) format
my @taxonomy_table;
foreach my $species (keys %taxon_hash) {
    my $node = $taxon_hash{$species};
    if (!$node) {
        print STDERR "WARNING: '$species' was not found in the taxonony tree\n";
        next;
	}
    # Copy taxonomy template and fill it in
    my %taxonomy = %rank_template;
    while ($node->name() ne 'root') {
        my $rank = $node->rank();
        if (exists($taxonomy{$rank})) {
            $taxonomy{$rank} = $node->name();
        }
        $node = $node->parent();
    }
    # Get the taxonomy in a TSV-like string, add GenomeDB name as the last
    # level/rank and append it to the taxonomy table
    my $taxonomy_str = join("\t", (map {$taxonomy{$_}} @target_ranks), $species);
    push(@taxonomy_table, $taxonomy_str);
}
# Sort the rows in alphabetical order
@taxonomy_table = sort(@taxonomy_table);
if ($group_ranks) {
    # Group ranks in taxonomy table, replacing duplicates by "''"
    # Get the first row as reference to see how many ranks are the same
    my @prev_taxonomy = split(/\t/, $taxonomy_table[0]);
    foreach my $i (1 .. $#taxonomy_table) {
        my @taxonomy = split(/\t/, $taxonomy_table[$i]);
        # Grouped version of the taxonomy (string)
        my $grouped_taxonomy = "";
        # Flag when two consecutive species have a not-null matching rank
        my $not_null_match_found = 0;
        # Loop over each rank to compare it with the value in the previous row
        foreach my $j (0 .. $num_ranks) {
            if ($prev_taxonomy[$j] eq $taxonomy[$j]) {
                # Get the flag to True if the ranks are not null and equal
                if ($taxonomy[$j] ne $null_value) {
                    $not_null_match_found = 1;
                }
            } else {
                # If there is a not-null match, all the ranks matched can be
                # replaced by "''". If not, it means all the ranks were missing
                # until this point, so leave the null value.
                my @values;
                if ($not_null_match_found) {
                    @values = ("''") x $j;
                } else {
                    @values = ($null_value) x $j;
                }
                # As soon as one rank is different, the remaining will be too
                $grouped_taxonomy .= join("\t", @values, @taxonomy[$j .. $num_ranks]);
                last;
            }
        }
        # Replace the string by its grouped version
        $taxonomy_table[$i] = $grouped_taxonomy;
        # Get the current taxonomy as reference for the next step
        @prev_taxonomy = @taxonomy;
    }
}
# Write the array into STDOUT with @target_ranks as headers
print join("\t", map(ucfirst, @target_ranks), 'Requested name'), "\n",
      join("\n", @taxonomy_table), "\n";


sub map_row_to_header {
    my ($line, $header) = @_;
    
    chomp $line;
    chomp $header;
    my @cols = split(/[\f\n\r\t]+/, $line);
    my @head_cols;
    if ( ref $header eq 'ARRAY' ) {
        @head_cols = @$header;
    } else {
        @head_cols = split(/[\f\n\r\t]+/, $header);
    }
    
    die "Number of columns in header do not match row" unless scalar @cols == scalar @head_cols;
    
    my $row;
    for ( my $i = 0; $i < scalar @cols; $i++ ) {
        $row->{$head_cols[$i]} = $cols[$i];
    }
    return $row;
}

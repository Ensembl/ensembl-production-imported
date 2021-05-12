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

package Bio::EnsEMBL::EGPipeline::RNAFeatures::DeleteUnattachedXref;

use strict;
use warnings;
use base ('Bio::EnsEMBL::EGPipeline::Xref::DeleteUnattachedXref');

sub delete_xrefs {
  my ($self, $dbh) = @_;

  my $sql = '
    DELETE xref.* FROM
      xref INNER JOIN
      external_db USING (external_db_id)
      LEFT OUTER JOIN xref_used USING (xref_id)
    WHERE
      xref_used.xref_id IS NULL AND
      db_name IN ("miRBase", "RFAM", "TRNASCAN_SE");
  ';
  my $sth = $dbh->prepare($sql);
  $sth->execute();
}

1;

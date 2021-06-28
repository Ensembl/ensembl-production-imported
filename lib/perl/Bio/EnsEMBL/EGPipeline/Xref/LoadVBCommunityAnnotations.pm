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

=pod

=head1 NAME

Bio::EnsEMBL::EGPipeline::Xref::LoadVBCommunityAnnotations

=head1 DESCRIPTION

Load VectorBase community-submitted xrefs from flat files.

=head1 Author

James Allen

=cut

package Bio::EnsEMBL::EGPipeline::Xref::LoadVBCommunityAnnotations;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::EGPipeline::Xref::LoadXref');

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use File::Spec::Functions qw(catdir);
use Path::Tiny qw(path);

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    'logic_name'           => 'brc4_community_annotation',
    'vb_external_db'       => 'BRC4_Community_Annotation',
    'citation_external_db' => 'PUBMED',
    'exclude_name_source'  => ['mirbase_gene', 'rfam_12.1_gene', 'trnascan_gene'],
    'exclude_desc_source'  => ['mirbase_gene', 'rfam_12.1_gene', 'trnascan_gene'],
  };
}

sub run {
  my ($self) = @_;
  my $db_type              = $self->param_required('db_type');
  my $logic_name           = $self->param_required('logic_name');
  my $vb_external_db       = $self->param_required('vb_external_db');
  my $citation_external_db = $self->param_required('citation_external_db');
  my $exclude_name_source  = $self->param_required('exclude_name_source');
  my $exclude_desc_source  = $self->param_required('exclude_desc_source');

  my $dba  = $self->get_DBAdaptor($db_type);
  my $aa   = $dba->get_adaptor('Analysis');
  my $ga   = $dba->get_adaptor('Gene');
  my $dbea = $dba->get_adaptor('DBEntry');

  my $analysis = $aa->fetch_by_logic_name($logic_name);

  foreach my $external_db ($vb_external_db, $citation_external_db) {
    $self->external_db_reset($dba, $external_db);
  }
  
  my %exclude_symbol   = map {$_ => 1} @{$exclude_name_source};
  my %exclude_desc     = map {$_ => 1} @{$exclude_desc_source};
  my $exclude_citation = {};
  
  my $xref_file = $self->fetch_xref_file();
  my ($annotations, $citations) = $self->parse_xref_file($xref_file);

  $self->add_symbols($ga, $dbea, $analysis, $vb_external_db, \%exclude_symbol, $annotations);

  $self->add_descriptions($ga, $analysis, $vb_external_db, \%exclude_desc, $annotations);
  
  $self->add_citations($ga, $dbea, $analysis, $citation_external_db, $exclude_citation, $citations);

  foreach my $external_db ($vb_external_db, $citation_external_db) {
    $self->external_db_update($dba, $external_db);
  }
}

sub fetch_xref_file {
  my ($self) = @_;
  my $pipeline_dir = $self->param_required('pipeline_dir');
  my $species      = $self->param_required('species');
  my $annot_dir    = $self->param_required('annotation_dir');
  
  my $file = catdir($annot_dir, "$species.tsv");

  return $file;
}

sub parse_xref_file {
  my ($self, $file) = @_;
  my %annotations;
  my %citations;
  
  my $path = path($file);
  my @lines = $path->lines({chomp => 1});
  
  # Remove header
  shift @lines;
  
  foreach my $line (@lines) {
    my (
      $curation_id,
      $stable_id,
      $status,
      $biotype,
      $model,
      $from,
      $to,
      $symbol,
      $synonym,
      $gene_family,
      $description,
      $submitter,
      $value,
      $species_name,
      $submitter_comments,
      $date,
    ) = split(/\t/, $line);
    
    $symbol = '' if $stable_id eq $symbol;
    
    if (not $status) {
      warn("No status for line '$line'");
      next;
    }
    
    $description = latinize($description);
    $symbol = latinize($symbol);
    
    if ($status eq 'ACTIVE') {
      if ($symbol ne '') {
        # Symbol and possibly empty description
        $annotations{$stable_id}{$symbol}{desc} = $description;
        if ($synonym ne '') {
          # Synonym exists for symbol
          $annotations{$stable_id}{$symbol}{synonyms}{$synonym}++;
        }
      } elsif ($synonym ne '') {
        # Synonym without a symbol => treat as a symbol
        unless (exists $annotations{$stable_id}{$synonym}) {
          $annotations{$stable_id}{$synonym}{desc} = $description;
        }
      } elsif ($description ne '') {
        # Description without a symbol
        $annotations{$stable_id}{$symbol}{desc} = $description;
      }
    } elsif ($status eq 'SYNONYM') {
      if ($symbol ne '') {
        if ($synonym ne '') {
          # Synonym exists for symbol
          $annotations{$stable_id}{$symbol}{synonyms}{$synonym}++;
        }
      } elsif ($synonym ne '') {
        # In lieu of explicit symbol, add synonym to any existing symbols
        my $added_synonym = 0;
        foreach my $symbol (keys %{$annotations{$stable_id}}) {
          if ($symbol ne $synonym) {
            $annotations{$stable_id}{$symbol}{synonyms}{$synonym}++;
            $added_synonym = 1;
          }
        }
        
        # Synonym without a symbol => treat as a symbol
        if (! $added_synonym) {
          unless (exists $annotations{$stable_id}{$synonym}) {
            $annotations{$stable_id}{$synonym}{desc} = $description;
          }
        }
      }
    } elsif ($status eq 'CITATION') {
      if ($value ne '') {
        # Citation exists
        $citations{$stable_id}{$value} = $submitter;
        if ($symbol ne '') {
          # Citation has a (possibly duplicated) symbol
          unless (exists $annotations{$stable_id}{$symbol}) {
            $annotations{$stable_id}{$symbol}{desc} = $description;
          }
        }
      }
    }
  }
  
  die "No annotations were loaded from $file" if not %annotations;
  
  $self->delete_unmatched_symbols(\%annotations);
              
  return (\%annotations, \%citations);
}

sub latinize {
  my ($str) = @_;
  
  my %letter_name = (
    'α' => 'alpha',
    'β' => 'beta',
    'γ' => 'gamma',
    'δ' => 'delta',
    'ε' => 'epsilon',
    'ζ' => 'zeta',
    'η' => 'eta',
    'θ' => 'theta',
    'ι' => 'iota',
    'κ' => 'kappa',
    'λ' => 'lamda',
    'μ' => 'mu',
    'ν' => 'nu',
    'ξ' => 'xi',
    'ο' => 'omicron',
    'π' => 'pi',
    'ρ' => 'rho',
    'ς' => 'sigma',
    'σ' => 'sigma',
    'τ' => 'tau',
    'υ' => 'upsilon',
    'φ' => 'phi',
    'χ' => 'chi',
    'ψ' => 'psi',
    'ω' => 'omega',
  );
  
  my $new_str = $str;
  for my $char (keys %letter_name) {
    $new_str =~ s/$char/$letter_name{$char}/gi;
  }
  
  if ($str ne $new_str) {
    warn("Renamed '$str' to '$new_str'\n");
  }
  
  return $new_str;
}

sub delete_unmatched_symbols {
  my ($self, $annotations) = @_;
  
  # Sometimes 'synonym' rows have a symbol with different case to
  # the corresponding 'active' symbol. Or 'synonym' rows exist
  # for deprecated data, or are assigned to the wrong stable ID.
  my %delete;
  foreach my $stable_id (keys %$annotations) {
    foreach my $symbol (keys %{$$annotations{$stable_id}}) {
      if (! defined $$annotations{$stable_id}{$symbol}{'desc'}) {
        $delete{$stable_id}{$symbol} = $$annotations{$stable_id}{$symbol}{synonyms};
        $self->warning("Unmatched symbol found with synonym: $symbol (".join(",", keys %{$delete{$stable_id}{$symbol}}).")");
      }
    }
  }
  
  foreach my $stable_id (keys %delete) {
    foreach my $bad_symbol (keys %{$delete{$stable_id}}) {
      foreach my $symbol (keys %{$$annotations{$stable_id}}) {
        if (lc($symbol) eq lc($bad_symbol) && $symbol ne $bad_symbol) {
          foreach my $synonym (keys %{$delete{$stable_id}{$bad_symbol}}) {
            $$annotations{$stable_id}{$symbol}{synonyms}{$synonym}++;
            $self->warning("Synonym for unmatched symbol moved: $bad_symbol => $symbol");
          }
        }
      }
    }
  }
  
  foreach my $stable_id (keys %delete) {
    foreach my $bad_symbol (keys %{$delete{$stable_id}}) {
      if (exists $$annotations{$stable_id}{$bad_symbol}) {
        delete $$annotations{$stable_id}{$bad_symbol};
      }
    }
  }      
}

sub add_symbols {
  my ($self, $ga, $dbea, $analysis, $external_db, $exclude, $annotations) = @_;
  
  my $synonym_sql = 'INSERT IGNORE INTO external_synonym VALUES (?, ?)';
  my $sth = $dbea->dbc->db_handle->prepare($synonym_sql);
  
  $self->remove_display_xrefs($ga, $external_db);

  foreach my $stable_id (keys %$annotations) {
    my $gene = $ga->fetch_by_stable_id($stable_id);
    if (defined $gene && ! exists($$exclude{$gene->analysis->logic_name})) {
      foreach my $symbol (keys %{$$annotations{$stable_id}}) {
        if ($symbol ne '') {
          my $desc = $$annotations{$stable_id}{$symbol}{desc};
          
          my $xref = $self->add_xref($symbol, $desc, $analysis, $external_db);
          
          $dbea->store($xref, $gene->dbID, 'Gene');

          if ($stable_id ne $xref->display_id) {
            $gene->display_xref($xref);
            $ga->update($gene);
          }
          
          # Synonyms are only added by the '$dbea->store' method if
          # the xref doesn't already exist.
          if (exists $$annotations{$stable_id}{$symbol}{synonyms}) {
            for my $synonym (keys %{$$annotations{$stable_id}{$symbol}{synonyms}}) {
              $sth->execute($xref->dbID, $synonym) or $self->throw("Failed to add synonym '$synonym'");
            }
          }
        }
      }
    }
  }
}

sub add_descriptions {
  my ($self, $ga, $analysis, $external_db, $exclude, $annotations) = @_;
  
  $self->remove_descriptions($ga, $external_db);
  
  my $db_display_name = $self->fetch_external_db_display($external_db);

  foreach my $stable_id (keys %$annotations) {
    my %descs;
    
    foreach my $symbol (keys %{$$annotations{$stable_id}}) {
      my $desc = $$annotations{$stable_id}{$symbol}{desc};
      if ($desc ne '') {
        $descs{$desc}++;
      }
    }
    
    my @descs = keys %descs;
    if (scalar(@descs) == 1) {
      my $gene = $ga->fetch_by_stable_id($stable_id);
      if (defined $gene && ! exists($$exclude{$gene->analysis->logic_name})) {
        my $desc = $descs[0];
        $gene->description($desc);
        $ga->update($gene);
      }
    } elsif (scalar(@descs) > 1) {
      $self->throw("Multiple descriptions for $stable_id: ".join(';', @descs));
    }
  }
}

sub add_citations {
  my ($self, $ga, $dbea, $analysis, $external_db, $exclude, $citations) = @_;

  foreach my $stable_id (keys %$citations) {
    my $gene = $ga->fetch_by_stable_id($stable_id);
    if (defined $gene && ! exists($$exclude{$gene->analysis->logic_name})) {
      foreach my $pubmed_id (keys %{$$citations{$stable_id}}) {
        my $submitter = $$citations{$stable_id}{$pubmed_id};
        my $desc = "PubMed ID $pubmed_id";
        $desc   .= " - Supplied by $submitter" if $submitter;
        
        my $xref = $self->add_xref($pubmed_id, $desc, $analysis, $external_db);
        $dbea->store($xref, $gene->dbID(), 'Gene');
      }
    }
  }
}

sub add_xref {
  my ($self, $acc, $desc, $analysis, $external_db) = @_;
  
  $desc = undef if $desc eq '';

  my $xref = Bio::EnsEMBL::DBEntry->new(
    -PRIMARY_ID  => $acc,
  	-DISPLAY_ID  => $acc,
    -DESCRIPTION => $desc,
    -VERSION     => undef,
  	-DBNAME      => $external_db,
  	-INFO_TYPE   => 'DIRECT',
  );
  $xref->analysis($analysis);

  return $xref;
}

sub remove_display_xrefs {
  my ($self, $dba, $external_db) = @_;

  my $sql = q/
    UPDATE
      gene g INNER JOIN
      xref x ON g.display_xref_id = x.xref_id INNER JOIN
      external_db edb USING (external_db_id)
    SET
      g.display_xref_id = NULL
    WHERE
      edb.db_name = ?
  /;
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($external_db);
}

sub remove_descriptions {
  my ($self, $dba, $external_db) = @_;

  my $sql = q/
    UPDATE
      gene g,
      external_db edb
    SET
      g.description = NULL
    WHERE
      edb.db_name = ? AND
      g.description LIKE CONCAT('%', ' [Source:', edb.db_display_name, '%')
  /;
  my $sth = $dba->dbc->db_handle->prepare($sql);
  $sth->execute($external_db);
}

1;

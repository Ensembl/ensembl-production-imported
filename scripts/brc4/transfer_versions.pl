#!/usr/env perl
use v5.14.00;
use strict;
use warnings;
use Carp;
use autodie qw(:all);
use Readonly;
use Getopt::Long qw(:config no_ignore_case);
use Log::Log4perl qw( :easy ); 
Log::Log4perl->easy_init($WARN); 
my $logger = get_logger(); 

use Bio::EnsEMBL::Registry;
use Try::Tiny;
use File::Path qw(make_path);

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

my $species = $opt{species};

# Features to transfer versions from
my @features = qw(gene);

# Get genes and transcripts from old database
my $registry = 'Bio::EnsEMBL::Registry';
my $reg_path = $opt{old_registry};
$registry->load_all($reg_path);

my %old_data;
for my $feat (@features) {
  $old_data{$feat} = get_feat_data($registry, $species, $feat);
  print STDERR scalar(keys %{$old_data{$feat}}) . " ${feat}s from old database\n";
}

# Reload registry, this time for new database
$reg_path = $opt{new_registry};
$registry->load_all($reg_path);

my %new_data;
for my $feat (@features) {
  $new_data{$feat} = get_feat_data($registry, $species, $feat);
  print STDERR scalar(keys %{$new_data{$feat}}) . " ${feat}s from new database\n";
}

# Transfer descriptions
if ($opt{descriptions}) {
  say STDERR "Gene descriptions transfer:";
  update_descriptions($registry, $species, $old_data{gene}, $opt{write});
}

if ($opt{events}) {
  # Get events from features differences
  my ($old_ids, $new_ids) = diff_events($old_data{gene}, $new_data{gene});

  # Load all events from the file
  my $del_events = load_deletes($opt{deletes}, $old_ids);
  my $file_events = load_events($opt{events}, $old_ids, $new_ids);

  my %events = (%$del_events, %$file_events);

  for my $event_type (sort keys %events) {
    my @feat_events = @{$events{$event_type}};
    say("Event: $event_type = " . scalar(@feat_events));
  }

  my $ga = $registry->get_adaptor($species, "core", 'gene');
  my $dbc = $ga->dbc;
  add_events(\%events, $old_data{gene}, $dbc, $opt{write});
}

# Transfer versions
if ($opt{versions}) {
  say STDERR "Gene versions transfer:";
  # Update the version for the features we want to transfer and update
  for my $feat (@features) {
    update_versions($registry, $species, $feat, $old_data{$feat}, $opt{write});
  }

  # Blanket version replacement for transcripts, translations, exons
  if ($opt{write}) {
    say STDERR "Transcripts, translations and exons version update";
    my $ga = $registry->get_adaptor($species, "core", 'gene');
    my $dbc = $ga->dbc;

    # Transfer the versions from genes to transcripts
    my $transcript_query = "UPDATE transcript LEFT JOIN gene USING(gene_id) SET transcript.version = gene.version";
    $dbc->do($transcript_query);

    # Also transfer the versions from transcripts to translations
    my $translation_query = "UPDATE translation LEFT JOIN transcript USING(transcript_id) SET translation.version = transcript.version";
    $dbc->do($translation_query);

    # And set the exons version to 1
    my $exon_query = "UPDATE exon SET version = 1";
    $dbc->do($exon_query);
  } else {
    say STDERR "Transcripts, translations and exons init versions were not updated (use --write to do so)";
  }
}

###############################################################################
sub load_deletes {
  my ($deletes_file, $old_ids) = @_;

  if (not $deletes_file) {
    return {};
  }

  my %events = (
    deleted => [],
  );
  open(my $deletes_fh, "<", $deletes_file);
  while (my $line = readline($deletes_fh)) {
    chomp $line;
    my $id = $line;

    if ($old_ids->{$id}) {
      push @{$events{deleted}}, {from => [$id], to => []};
      delete $old_ids->{$id};
    } else {
      say("Warning: '$id' in the deletes file, but not deleted in the new core");
    }
  }

  return \%events;
}

sub load_events {
  my ($events_file, $old_ids, $new_ids) = @_;

  if (not $events_file) {
    return {};
  }

  my %events = (
    change => [],
    new => [],
    split => [],
    merge => [],
  );
  my %merge_to = ();
  my %split_from = ();
  open(my $events_fh, "<", $events_file);
  while (my $line = readline($events_fh)) {
    chomp $line;
    my ($id1, $event_name, $id2) = split("\t", $line);

    # New gene
    if ($id1 and not $id2) {
      push @{$events{new}}, {from => [], to => [$id1]};
      if ($new_ids->{$id1}) {
        delete $new_ids->{$id1};
      }
    # Changed gene
    } elsif ($id1 eq $id2) {
      push @{$events{change}}, {from => [$id1], to => [$id1]};
    # Merge
    } elsif ($event_name eq 'merge_gene') {
      if (not $merge_to{$id1}) {
        $merge_to{$id1} = [];
      }
      push @{$merge_to{$id1}}, $id2;
    # Split
    } elsif ($event_name eq 'split_gene') {
      if (not $split_from{$id2}) {
        $split_from{$id2} = [];
      }
      push @{$split_from{$id2}}, $id1;
    }
    else {
      say "Unsupported event '$event_name'? $line";
    }
  }

  while (my ($merge_to_id, $merge_from) = each %merge_to) {
    push @{$events{merge}}, {from => $merge_from, to => [$merge_to_id]};

    if ($new_ids->{$merge_to_id}) {
      delete $new_ids->{$merge_to_id};
    }
    for my $merge_from_id (@$merge_from) {
      if ($old_ids->{$merge_from_id}) {
        delete $old_ids->{$merge_from_id};
      }
    }
  }
  while (my ($split_from_id, $split_to) = each %split_from) {
    push @{$events{split}}, {from => [$split_from_id], to => $split_to};

    if ($old_ids->{$split_from_id}) {
      delete $old_ids->{$split_from_id};
    }
    for my $split_to_id (@$split_to) {
      if ($new_ids->{$split_to_id}) {
        delete $new_ids->{$split_to_id};
      }
    }
  }

  if (%$new_ids) {
    say(scalar(%$new_ids) . " new ids not in the event file: " . join("; ", sort keys %$new_ids));
  }

  if (%$old_ids) {
    say(scalar(%$old_ids) . " old ids not in the event file: " . join("; ", sort keys %$old_ids));
  }

  close($events_fh);

  return \%events;
}

sub diff_events {
  my ($old, $new) = @_;

  my %events = (
    new => [],
    deleted => [],
  );

  my %old_ids = map { $_ => 1 } keys %$old;
  my %new_ids = map { $_ => 1 } keys %$new;

  for my $old_id (keys %old_ids) {
    if ($new_ids{$old_id}) {
      delete $old_ids{$old_id};
      delete $new_ids{$old_id};
    }
  }

  return \%old_ids, \%new_ids;
}

sub get_feat_data {
  my ($registry, $species, $feature) = @_;
  
  my %feats;
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;
    my $description;
    if ($feature eq 'gene') {
      $description = $feat->description;
    }
    $feats{$id} = { version => $version, description => $description };
  }
  
  return \%feats;
}

sub update_versions {
  my ($registry, $species, $feature, $old_feats, $update) = @_;
  
  my $update_count = 0;
  my $new_count = 0;
  
  my $fa = $registry->get_adaptor($species, "core", $feature);
  for my $feat (@{$fa->fetch_all}) {
    my $id = $feat->stable_id;
    my $version = $feat->version;

    if (not defined $version) {
      my $old_feat = $old_feats->{$id};
      
      if (not $old_feat) {
        $new_count++;
        next;
      }
      
      my $old_version = $old_feat->{version};
      my $new_version = 1;

      if (defined $old_version) {
        $new_version = $old_version + 1;
        say STDERR "Updated $feature $id must be upgrade from $old_version to $new_version" if $opt{verbose};
        $update_count++;
      } else {
        say STDERR "New $feature $id must be initialized to 1" if $opt{verbose};
        $new_count++;
      }

      if ($update) {
        $feat->version($new_version);
        $fa->update($feat);
      }
    }
  }
  
  print STDERR "$update_count $feature version updated\n";
  print STDERR "$new_count $feature version to be initialized\n";
  print STDERR "(Use --write to make the changes to the database)\n" if $update_count + $new_count > 0 and not $update;
}

sub update_descriptions {
  my ($registry, $species, $old_genes, $update) = @_;
  
  my $update_count = 0;
  my $empty_count = 0;
  my $new_count = 0;
  
  my $ga = $registry->get_adaptor($species, "core", 'gene');
  for my $gene (@{$ga->fetch_all}) {
    my $id = $gene->stable_id;
    my $description = $gene->description;

    if (not defined $description) {
      my $old_gene = $old_genes->{$id};
      if (not $old_gene) {
        $new_count++;
        next;
      }
      my $old_description = $old_gene->{description};

      if (defined $old_description) {
        my $new_description = $old_description;
        say STDERR "Transfer gene $id description: $new_description" if $opt{verbose};
        $update_count++;

        if ($update) {
          $gene->description($new_description);
          $ga->update($gene);
        }
      } else {
        $empty_count++;
      }
    }
  }
  
  print STDERR "$update_count gene descriptions transferred\n";
  print STDERR "$empty_count genes without description remain\n";
  print STDERR "$empty_count new genes, without description\n";
  print STDERR "(Use --write to update the descriptions in the database)\n" if $update_count > 0 and not $update;
}

sub add_events {
  my ($events, $data, $dbc, $write) = @_;

  my $feat = 'gene';
  my $metadata = {};
  my $mapping_id = add_mapping_session($dbc, $metadata);

  for my $event_name (keys %$events) {
    say "Storing events $event_name...";
    my @events = @{$events->{$event_name}};
    for my $event (@events) {
      my @from = @{$event->{from}};
      my @to = @{$event->{to}};

      if ($event_name eq 'new') {
        insert_event($dbc, $write, [$mapping_id, $feat, undef, undef, $to[0], 1]);
      }
      elsif ($event_name eq 'deleted') {
        my $id = $from[0];
        my $version = $data->{$id}->{version};
        insert_event($dbc, $write, [$mapping_id, $feat, $id, $version, undef, undef]);
      }
      elsif ($event_name eq 'change') {
        my $id = $from[0];
        my $old_version = $data->{$id}->{version};
        my $new_version = $old_version + 1;
        insert_event($dbc, $write, [$mapping_id, $feat, $id, $old_version, $id, $new_version]);
      }
      elsif ($event_name eq 'split') {
        my $id = $from[0];
        my $old_version = $data->{$id}->{version};
        {
          my @values = ($mapping_id, $feat, $id, $old_version, undef, undef);
          insert_event($dbc, $write, [$mapping_id, $feat, $id, $old_version, undef, undef]);
        }
        for my $split_id (@to) {
          my @values = ($mapping_id, $feat, $id, $old_version, $split_id, 1);
          insert_event($dbc, $write, [$mapping_id, $feat, $id, $old_version, $split_id, 1]);
        }
      }
      elsif ($event_name eq 'merge') {
        my $id = $to[0];
        for my $merge_id (@from) {
          my $old_version = $data->{$merge_id}->{version};
          insert_event($dbc, $write, [$mapping_id, $feat, $merge_id, $old_version, $id, 1]);
          insert_event($dbc, $write, [$mapping_id, $feat, $merge_id, $old_version, undef, undef]);
        }
      }
    }
  }
}

sub insert_event {
  my ($dbc, $write, $values) = @_;
  my $sql = "INSERT INTO stable_id_event(mapping_session_id, old_stable_id, old_version, new_stable_id, new_name, type) VALUES(?,?,?,?,?,?)";
  my $sth = $dbc->prepare($sql);
  say("Insert values: " . join(", ", map { $_ // 'undef' } @$values));

}

sub add_mapping_session {
  my ($dbc, $metadata) = @_;
  # TODO
  return 1;
}

###############################################################################
# Parameters and usage
sub usage {
  my $error = shift;
  my $help = '';
  if ($error) {
    $help = "[ $error ]\n";
  }
  $help .= <<'EOF';
    Transfer the gene versions and descriptions to a patched database

    --old_registry <path> : Ensembl registry
    --new_registry <path> : Ensembl registry
    --species <str>   : production_name of one species

    --descriptions    : Transfer the gene descriptions

    Use either of those:
    --versions        : Transfer the gene versions, and init the others
    --events <path>   : Path to an events file, to update the history and versions
    --deletes <path>  : Path to a list of deleted genes (to use with the events file)
    
    --write           : Do the actual changes (default is no changes to the database)
    
    --help            : show this help message
    --verbose         : show detailed progress
    --debug           : show even more information (for debugging purposes)
EOF
  print STDERR "$help\n";
  exit(1);
}

sub opt_check {
  my %opt = ();
  GetOptions(\%opt,
    "old_registry=s",
    "new_registry=s",
    "species=s",
    "descriptions",
    "versions",
    "events=s",
    "deletes=s",
    "write",
    "help",
    "verbose",
    "debug",
  );

  usage("Old registry needed") if not $opt{old_registry};
  usage("New registry needed") if not $opt{new_registry};
  usage("Species needed") if not $opt{species};

  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__


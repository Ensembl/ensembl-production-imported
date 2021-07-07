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
use Capture::Tiny ':all';

use Bio::EnsEMBL::ApiVersion;
use Try::Tiny;
use File::Spec::Functions;
use DBI;

###############################################################################
# MAIN
# Get command line args
our %opt = %{ opt_check() };

my $update = $opt{update};
my $server = $opt{server};

my $release = software_version();
die("Can't find the release of the current Perl API!") if not $release;
my $ens_dir = $ENV{ENSEMBL_ROOT_DIR};
die("Can't find the ENSEMBL_ROOT_DIR of the current Perl API!") if not $ens_dir;
say STDERR "Current release version is $release";
say STDERR "Ensembl root dir is $ens_dir";

my $core_sql_dir = catdir($ens_dir, 'ensembl', 'sql');
my $var_sql_dir  = catdir($ens_dir, 'ensembl-variation', 'sql');
die("Can't find core sql directory in $core_sql_dir") if not -e $core_sql_dir;
die("Can't find variation sql directory in $var_sql_dir") if not -e $var_sql_dir;

# Open server
my @connect_list = ("DBI:mysql:", "host=$opt{host}", "port=$opt{port}");
my $connect_string = join(";", @connect_list);
my $dbh = DBI->connect($connect_string, $opt{user}, $opt{pass}) or die("Can't connect to the server");

my $databases = get_databases($dbh);
say scalar(@$databases) . " databases in $opt{host}";
my @core_dbs = grep { $_ =~ /_core_\d+_\d+_\d+/ } @$databases;
my @var_dbs  = grep { $_ =~ /_variation_\d+_\d+_\d+/ } @$databases;
say scalar(@core_dbs) . " core databases";
say scalar(@var_dbs) . " variation databases";

update_dbs($server, \@core_dbs, $release, $core_sql_dir, $update);
update_dbs($server, \@var_dbs,  $release, $var_sql_dir,  $update);

###############################################################################
sub get_databases {
  my ($dbh) = @_;
  
  my $db_array = $dbh->selectcol_arrayref("SHOW DATABASES;");
  
  return $db_array;
}

sub update_dbs {
  my ($server, $dbs, $cur_release, $sql_dir, $update) = @_;
  
  for my $db (@$dbs) {
    my $db_release;
    try {
      $db_release = get_db_release($db);
    } catch {
      warn("Can't get release from db $db: $_");
      next;
    };

    my ($db_name_release) = ($db =~ /_(\d+)_\d+$/);
    
    # Check that release versions are ok
    if (not defined $db_release) {
      warn("$db\tDB release from db is not defined");
      next;
    }
    if (not defined $db_name_release) {
      warn("$db\tCan't get db release from db name");
      next;
    }
    if ($db_release != $db_name_release) {
      warn("$db\tRelease version differs between db name and meta table: $db_release vs $db_name_release\n");
      next;
    }

    # Check if db needs to be updated
    if ($cur_release >= $db_release) {
      my $new_db = update_db($server, $db, $cur_release, $update);
      update_db_release($server, $new_db, $db_release, $cur_release, $sql_dir, $update);
    } else {
      say STDERR "$db\tRelease is newer than current: $db_release to $cur_release" if $opt{verbose};
    }
  }
}

sub get_db_release {
  my ($db_name) = @_;
  
  my @connect_list = ("DBI:mysql:$db_name", "host=$opt{host}", "port=$opt{port}");
  my $connect_string = join(";", @connect_list);
  my $dbh = DBI->connect($connect_string, $opt{user}, $opt{pass}) or die("Can't connect to the server");
  
  my $col = $dbh->selectcol_arrayref("SELECT meta_value FROM meta WHERE meta_key='schema_version';");
  my $db_release = $col->[0];
  
  $dbh->disconnect();
  
  return $db_release;
}

sub update_db {
  my ($server, $db, $cur_release, $update) = @_;

  my $new_db = $db;
  $new_db =~ s/_(core|variation)_(\d+)_(\d+)_(\d+)$/_$1_$2_${cur_release}_$4/;

  # One time prefix removal step, to be removed
  my $remove_prefix = 1;
  if ($remove_prefix) {
    $new_db =~ s/^[^_]+_//;
  }

  if ($db eq $new_db) {
    say "$db\tNo update to the name" if $opt{verbose};
  } else {
    rename_db($server, $db, $new_db, $update);
  }
  return $new_db;
}

sub rename_db {
  my ($server, $old_name, $new_name, $update) = @_;
  
  my $command = "rename_db $server $old_name $new_name";
  say STDERR "Use command '$command'" if $opt{debug};
  
  system($command) if $update;
}

sub update_db_release {
  my ($server, $new_db, $db_release, $cur_release, $sql_dir, $update) = @_;
  
  for (my $rel = $db_release+1; $rel <= $cur_release; ++$rel) {
    # Get sql files
    my @sql_files = get_sql_files($sql_dir, $rel);
    for my $file (@sql_files) {
      my $command = "$server $new_db < $file";
      say STDERR "Use command '$command'" if $opt{debug};
      
      system($command) if $update;
    }
  }
}

sub get_sql_files {
  my ($sql_dir, $release) = @_;
  
  opendir(my $DIR, $sql_dir);
  my @files;
  while (my $file = readdir $DIR) {
    push @files, catfile($sql_dir, $file) if $file =~ /^patch_${release}_/;
  }
  return @files;
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
    Update the version of all databases in a given server

    --host
    --port
    --user
    --pass            : server parameters to use
    --server <str>    : server short name from mysql-cmd
    
    --update          : Actually rename databases
    
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
    "host=s",
    "port=s",
    "user=s",
    "pass=s",
    "server=s",
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Server host needed") if not $opt{host};
  usage("Server port needed") if not $opt{port};
  usage("Server user needed") if not $opt{user};
  usage("Server short name needed") if not $opt{server};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__


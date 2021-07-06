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
my %opt = %{ opt_check() };

my $update = $opt{update};

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

update_dbs(\@core_dbs, $release, $core_sql_dir, $update);
update_dbs(\@var_dbs,  $release, $var_sql_dir,  $update);

###############################################################################
sub get_databases {
  my ($dbh) = @_;
  
  my $db_array = $dbh->selectcol_arrayref("SHOW DATABASES;");
  
  return $db_array;
}

sub update_dbs {
  my ($dbs, $cur_release, $sql_dir, $update) = @_;
  
  for my $db (@$dbs) {
    try {
      my $db_release = get_db_release($db);
      next if not $db_release;

      my ($db_name_release) = ($db =~ /_(\d+)_\d+$/);

      if ($db_release != $db_name_release) {
        warn("$db\tRelease version differs between db name and meta table: $db_release vs $db_name_release\n");
        die;
      }

      if ($cur_release > $db_release) {
        say STDERR "$db\tRelease is older than current: need updating from $db_release to $cur_release";
      } elsif($cur_release == $db_release) {
        say STDERR "$db\tUp to date";
      } else {
        say STDERR "$db\tRelease is newer than current: $db_release to $cur_release";
      }
      #lastA
    } catch {
      warn("Can't get release from db $db");
      next;
    };
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
    "update",
    "help",
    "verbose",
    "debug",
  );

  usage("Server host needed") if not $opt{host};
  usage("Server port needed") if not $opt{port};
  usage("Server user needed") if not $opt{user};
  usage("Server pass needed") if not $opt{pass};
  usage()                if $opt{help};
  Log::Log4perl->easy_init($INFO) if $opt{verbose};
  Log::Log4perl->easy_init($DEBUG) if $opt{debug};
  return \%opt;
}

__END__


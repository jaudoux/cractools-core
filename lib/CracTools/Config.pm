package CracTools::Config;
# ABSTRACT: Manage and access CracTools configuration file

use strict;
use warnings;
use POSIX;

use Config::FileManager 1.6;
use Config::Simple;
use File::Basename;
use File::HomeDir;
use Carp;

use CracTools;

require Exporter;
our @ISA = qw(Exporter);

our $CONFIG_NAME = $CracTools::PACKAGE_NAME.".cfg";

our @EXPORT = qw(LoadConfig PrintVersion getConfVar);

our %config;

#load config
our $cfg = new Config::FileManager(
    "toolname" => $CracTools::PACKAGE_NAME, # Mandatory
    "version" => "$CracTools::VERSION", # Not mandatory
    "filename" => "$CONFIG_NAME", # Not mandatory
    "paths" => [
		".",
    File::HomeDir->my_home,
		".".$CracTools::PACKAGE_NAME,
		"__APPDIR__",
		"/usr/local/etc/".$CracTools::PACKAGE_NAME,
		"/etc/".$CracTools::PACKAGE_NAME
	       ], # Not mandatory
    "interactive" => 1, # Not mandatory
    );

my $default_content = "# Default configuration file __VERSION__\n#\n\n";

$cfg->defaultContent($default_content);

=head2 PrintVersion

Print (in an uniformized way) the version information of the CracUtil script.

Usage:

  PrintVersion();

=cut

sub PrintVersion() {
  printf( "Script '%s' from %s v. %s (%s v. %s)\n",
	  basename($0),
	  $CracTools::PACKAGE_NAME, $CracTools::VERSION,
	  $CracTools::PACKAGE_NAME, $CracTools::VERSION);
}

=head2 LoadConfig

Usage:

  my $default_cfg = LoadConfig(); # Use (and get) default config file (see L<Config::FileManager>)
  LoadConfig("your_file.conf");

=cut

sub LoadConfig(;$) {
    my ($config_file) = @_;
    if (!defined $config_file) {
	$cfg->update();
	$config_file = $cfg->getPath();
    }
    Config::Simple->import_from($config_file, \%config);
    return $config_file;
}

=head2 getConfVar

Usage:

  my $var = getConfVar("GENOME"); 

=cut

sub getConfVar(;$) {
  my $var_name = shift;
  my $die = shift;
  if(defined $config{$var_name}) {
    return $config{$var_name};
  } else {
    if(defined $die && $die eq 1) {
      croak("Config variable \"$var_name\" not found.");
    } else {
      return undef;
    }
  }
}

1; 

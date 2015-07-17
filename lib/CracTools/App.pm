use strict;
use warnings;
package CracTools::App;
# ABSTRACT: CracTools App::Cmd

use App::Cmd::Setup -app;

sub global_opt_spec {
  return (
    [ "help", "log additional output" ],
  );
}

1;

package CracTools::Const;
# ABSTRACT: Constants for the CracTools-core

use strict;
use warnings;
use Exporter qw(import);

=head1 SYNOPSIS

  # get a constant variable
  my $NA = $CracTools::Const::NOT_AVAILABLE;

=head1 DESCRIPTION

This module contains some constants that are defined for all the
CracTools pipelines.

=head1 CONSTANTS

=over

=item NOT_AVAILABLE

=cut

our $NOT_AVAILABLE = 'NA';

=item NUCLEOTIDES => [ A, C, G, T ]

=cut

our $NUCLEOTIDES = [A, C, G, T ];

=back

=cut

1;

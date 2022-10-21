#!/usr/bin/perl -w

my $species = shift;

while(<>) {
  if (/$species/i) {
    print $_;
    $_ = <>;
    print $_;
  }
}

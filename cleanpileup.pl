#!/usr/bin/perl
my $file = $ARGV[0];
open MYFILE, "<", $file;
while (<MYFILE>){
  s/\^.//g; #this removes references to read start qualities
  while (m/[+-]([0-9]+)/){
    s/[+-]$1[ACGTNacgtn]{$1}/X/g
    }; #this masks adjacent indels with a X
  print "$_";
};

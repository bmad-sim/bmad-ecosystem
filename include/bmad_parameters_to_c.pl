#!/usr/bin/perl

$file_out = "bmad_parameters.h";
open (F_OUT, ">$file_out") || die ("Cannot Open File: F_OUT\n");
print F_OUT "#ifndef BMAD_PARAMETERS\n\n";
print F_OUT "namespace Bmad {\n";

&searchit ("../modules/bmad_struct.f90");
&searchit ("../modules/bmad_basic_mod.f90");

print F_OUT "}\n\n";
print F_OUT "#define BMAD_PARAMETERS\n";
print F_OUT "#endif\n";

print "Created: $file_out\n";

#------------------------------------------------------------

sub searchit {

my($file_in) = @_;
open (F_IN, $file_in) || die ("Cannot Open File: $file_in\n");

$params_here = 0;


while (<F_IN>) {
  chomp;

  if (/\!/) {s/\!\Q$'\E//;}
  if (!/integer, *parameter/ && !$params_here) {next;}
  if (m@\(/@) {next;}

  $_ = uc;  # upper case
  s/INTEGER, *PARAMETER *::/  const int/;
  s/\$//g;

  if (/\&/) {
    s/\&//;
    $params_here = 1;
    s/\s*$/\n/;
  } else {
    $params_here = 0;
    s/\s*$/\;\n/;
  }

  print F_OUT;

}

}


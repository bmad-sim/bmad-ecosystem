#!/usr/bin/perl

`../../production/bin/element_attributes tex`;
`pdflatex bmad`;
`makeindex bmad`;
`makeindex bmad.rdx -o bmad.rnd`;
`pdflatex bmad`;

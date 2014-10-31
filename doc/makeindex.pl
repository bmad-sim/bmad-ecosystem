#!/usr/bin/perl

`pdflatex bmad`;
`makeindex bmad`;
`makeindex bmad.rdx -o bmad.rnd`;
`pdflatex bmad`;

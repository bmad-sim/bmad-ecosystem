#!/usr/bin/perl

# To web page

`pdflatex ibs_ring`; 
`chmod g+w ibs_ring.pdf`;
`cp ibs_ring.pdf /nfs/classe/www/html/bmad/ibs_ring.pdf`;
`cp ibs_ring.pdf ../../../bmad-doc/other_manuals/`;

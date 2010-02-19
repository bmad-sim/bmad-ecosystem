#!/usr/bin/perl

`latex synrad3d`; 
`dvips -o synrad3d.ps -z synrad3d`; 
`ps2pdf synrad3d.ps`;
`scp synrad3d.pdf   dcs\@lnx209.lns.cornell.edu:/home/dcs/public_html/`;

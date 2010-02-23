#!/usr/bin/perl

`latex synrad3d`; 
`dvips -o synrad3d.ps -z synrad3d`; 
`ps2pdf synrad3d.ps`;
`cp synrad3d.pdf /home/dcs/public_html/`;

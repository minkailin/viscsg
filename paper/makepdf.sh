#!/bin/bash
latex paper3d.tex
dvips -o paper3d.ps paper3d.dvi
ps2pdf -sPAPERSIZE=a4 paper3d.ps
echo 'done'

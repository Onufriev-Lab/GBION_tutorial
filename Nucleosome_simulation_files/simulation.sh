#!/bin/bash

pmemd.cuda -O -i min.in -o min.out -p nucleosome.top -c nucleosome.crd -r min.ncrst -inf min.mdinfo -ref nucleosome.crd

pmemd.cuda -O -i heat.in -o heat.out -p nucleosome.top -c min.ncrst -r heat.ncrst -x heat.nc -inf heat.mdinfo -ref nucleosome.crd

pmemd.cuda -O -i equil.in -o equil.out -p nucleosome.top -c heat.ncrst -r equil.ncrst -x equil.nc -inf equil.mdinfo -ref nucleosome.crd

pmemd.cuda -O -i prod.in -o prod.out -p nucleosome.top -c equil.ncrst -r prod.ncrst -x prod.nc -inf prod.mdinfo

ambpdb -p nucleosome.top -c prod.ncrst > prod.pdb

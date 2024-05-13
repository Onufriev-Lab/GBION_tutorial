#!/bin/bash

pmemd.cuda -O -i min.in -o min.out -p dna.top -c dna.crd -r min.ncrst -inf min.mdinfo -ref dna.crd

pmemd.cuda -O -i heat.in -o heat.out -p dna.top -c min.ncrst -r heat.ncrst -x heat.trj -inf heat.mdinfo -ref dna.crd

pmemd.cuda -O -i equil.in -o equil.out -p dna.top -c heat.ncrst -r equil.ncrst -x equil.trj -inf equil.mdinfo

pmemd.cuda -O -i prod.in -o prod.out -p dna.top -c equil.ncrst -r prod.ncrst -x prod.trj -inf prod.mdinfo

ambpdb -p dna.top -c prod.ncrst > prod.pdb
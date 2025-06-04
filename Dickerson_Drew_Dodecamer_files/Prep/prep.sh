#!/bin/bash

tleap -f tleap.script

python disang.py

cp dna.* ../

cp disang_NaCl.txt ../

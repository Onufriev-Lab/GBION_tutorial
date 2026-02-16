#!/bin/bash

tleap -f tleap.script

python disang.py

cp nucleosome.* ../

cp disang_KCl.txt ../

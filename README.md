# Tutorial GBION model

This tutorial aims to demonstrate simulation of a DNA and nucleosome in implicit solvent in combination with explicit ions. We assume that user has basic skills of running MD simulations using AMBER package. It also requires installing python with NumPy and Matplotlib libraries and CHIMERAX([link for downloading](https://www.cgl.ucsf.edu/chimerax/download.html)).

The outline of this tutorial:
1. Brief introduction into GBION model
2. **Part 1**: Instruction for simulation of Dickerson-Drew dodecamer using the model
	1. Preparing the topology, initial coordinates and restraints files for simulation
	2. Energy minimization
	3. Heating the system
	4. Equilibration of the DNA in ionic solution
	5. Production run
	6. Analysis of the DNA stability
	7. Visualization of trajectory using CHIMERAX
3. **Part 2**: Instruction for simulation of the nucleosome using the model
	1. Preparing the topology, initial coordinates and restraints files for simulation
	2. Energy minimization
	3. Heating the system
	4. Histone tails equilibration
	5. Production run
	6. Analysis of the nucleosome core stability
	7. Visualization of trajectory using CHIMERAX
## 1. Brief introduction
Large systems, such as nucleosomes, require 5- or 6-digits numbers of water molecules for the proper simulation with explicit water model. Another approach to treating water (implicit water model) allows taking water into account as a continuous environment around the molecule. The solvation free energy of the solvated molecules defines their behavior in the solution. The solvation energy $\Delta G_{solv}$ can be represented as a sum of electrostatic $\Delta G_{el}$ and non-polar $\Delta G_{np}$ contributions:

$\Delta G_{solv} = \Delta G_{el} + \Delta G_{np}$ 

which are estimated independently in most cases. One of the most popular approximation to calculating $\Delta G_{np}$ is based on the assumption that $\Delta G_{np}$ is proportional to the solvent-accessible surface area (SASA).

One of the most widely used approximation for calculation $\Delta G_{el}$ is the generalized Born (GB) model. The canonical GB approximation is based on the equation originally proposed by Still *et al.*:

$\Delta G_{el} = -\frac{1}{2} \left(\frac{1}{\epsilon_{in}} -\frac{1}{\epsilon_{out}}\right)\sum_{i,j} \frac{q_{i}q_{j}}{\sqrt{d^2_{ij} + R_{i}R_{j}\exp\left(-d^2_{ij}/(4 R_i R_j)\right)}}$,

where $\epsilon_{in}$ and $\epsilon_{out}$ are the dielectric constants of the solute and the solvent, respectively, $d_{ij}$ is the distance between solute atoms $i$ and $j$, and $q_{i}$ are the atomic charges. The key parameters modulating the interaction energy are *the effective Born radii* $R_i$. $R_i^{-1}$ characterizes the average degree of solvent exposure of atom $i$.

The GB model describes behaviour of the topologically connected structures really well, but it does not take into account descrete ions around solute. The GBION model was developed as an extension of GB model to simulate ions around DNA with implicit solvent. It implies additional coefficients to the GB equation:

$\Delta G_{el} = -\frac{1}{2} \left(\frac{1}{\mathbf{\epsilon_{in}(a,b)}} -\frac{1}{\epsilon_{out}}\right)\sum_{i,j} \frac{q_{i}q_{j}}{\sqrt{d^2_{ij} + R_{i}R_{j}\exp\left(-d^2_{ij}/(\mathbf {\gamma(a,b)} R_i R_j)\right)}}$

The expression above emphasizes the main idea of the GBION model that the functional form of charge-charge interaction is different for charges that are connected through the solute or the solvent. This is achieved by variations of $\mathbf{\epsilon_{in}}$ and $\mathbf{\gamma(a,b)}$  for different pairs of interacting atoms separately: solute-solute, solute-ion and ion-ion.

## 2. SImulation of Dickerson-Drew B-DNA Dodecamer using GBION model

### 2.1 Preparation of the system for simulations

#### 2.1.1 Preparation of the DNA for simulations

1. **Download “1BNA” (Dickerson-Drew dodecamer)** from RCSB: [https://www.rcsb.org/structure/1BNA](https://www.rcsb.org/structure/1BNA) or use the bundled file in `Dickerson_Drew_Dodecamer_files/Prep/1bna.pdb`.
2. **Remove waters (if you downloaded from RCSB):**
    - Open `1bna.pdb` in Chimera.
    - Select → Residue → HOH (select all water molecules).
    - Actions → Atoms/Bonds → Delete.
    - File → Save, save as “1bna.pdb”

The version in `Dickerson_Drew_Dodecamer_files/Prep` is already water-free.
#### 2.1.2 Building topology + placing ions using tleap
In directory `Dickerson_Drew_Dodecamer_files/Prep`, see a file (or create your own) named `tleap.script` with:
```
source leaprc.DNA.OL15
loadoff atomic_ions.lib
source leaprc.water.opc
set default PBradii mbondi3
mol = loadpdb 1bna.pdb
addions mol Na+ 36
addions mol Cl- 14
saveamberparm mol dna.top dna.crd
savepdb mol dna.pdb
quit
```

`source leaprc.DNA.OL15` – load the OL15 DNA force field.

`loadoff atomic_ions.lib` & `source leaprc.water.opc`: load ion parameters.

`set default PBradii mbondi3` – use mbondi3 radii (suitable for implicit solvent we use).

`mol = loadpdb 1bna.pdb` – loading structure of DNA

`addions mol Na+ 36` & `addions mol Cl- 14`: add 36 sodium and 14 chloride ions to neutralize and mimic roughly 150 mM.

`saveamberparm mol dna.top dna.crd` – saving topology and initial coordinates of the structure for simulation

`savepdb mol dna.pdb` – saving the system to PDB file

To run the script, change working directory to `Dickerson_Drew_Dodecamer_files/Prep`and type in command line:

`tleap -f tleap.script`

Output:  
* `dna.top` (AMBER topology)
* `dna.crd` (initial coordinates)
* `dna.pdb` (PDB file suitable for visualization of the initial state of the system)
#### 2.1.3 Generating distance-restraints for ions (disang.py)

Implicit‐solvent MD has no periodic box—ions would drift away. We use distance‐restraints to keep ions near the DNA. In the same `Prep/` directory type:

`python disang.py` – this provided script finds center of mass of a specified molecule, finds specified number of atoms closest to the center of mass and produces file with distance restraints for ions around the molecule. The restraints would limit the distance from the center of mass of spcified number of center atoms to each ion. One can specify parameters below to apply the script to their own molecules:

`Molecule_file = '1bna.pdb'` – specifies PDB-file with simulated biomolecule

`namedisang = "disang_NaCl.txt"` – specifies output restraints file

`Number_of_atoms_for_com = 10` – specifies number of atoms closest to the center of mass of the molecule to restraint ions with

`number_of_ions= 50` – specifies number of ions to restrain

`starting_index_of_ion= 759` – atom index of the first ion to restrain.

**line 50:** `        disang.write("  iat=-1,-1,r1=0.0,r2=0.0,r3=40,r4=50,rk2=0.0, rk3=20.0,\n")` – r3, r4 and rk3 specify parameters explained below.

This produces `disang_NaCl.txt`, containing blocks like:
```
&rst
iresid=0,
  iat=-1,-1,r1=0.0,r2=0.0,r3=40,r4=50,rk2=0.0, rk3=20.0,
igr1=588,545,209,186,166,153,140,139,134,108,igr2=808
 /
```

- `iresid=0`: atom indices are specified directly.
- `iat=-1,-1`: instructs AMBER to read `igr1` and `igr2` groups, then restrain their centers of mass.  
    • `igr1`: (set of 10 atoms) → DNA’s geometric center.  
    • `igr2`: index of an ion.
- `r1=0.0, r2=0.0, r3=40, r4=50`: define a flat-bottom restraint (0–40 Å flat, then a parabola up to 50 Å).
- `rk3=20.0`: force constant (20 kcal · mol<sup>−1</sup> · Å<sup>−2</sup>) as atoms go from 40→50 Å.

In general the graphs for restraining force on the distance looks like this:
![restraints](Pictures/restraints.png)
From 0 to `r1` force linearly depends on the distance, from `r1` to `r2` parabolically, `r2-r3` is a flat region, from `r3` to `r4`  – parabolically, from `r4` to $\infty$  – linearly. In our case 0-`r3` is a flat region:
![restraints_flat](Pictures/restraints_flat.png)

Copy `dna.*` and `disang_NaCl.txt` up one level:

`cp dna.* ../`

`cp *.txt ../`

If you prefer, run the provided script:

`bash prep.sh`

This will auto-generate `dna.top`,`dna.crd`,`dna.pdb` and `disang_NaCl.txt` and copy up.

The simulation stages and analysis are described below. Script `simulation.sh` generates all the trajectories, but does not provide any analysis.
### 2.2 Energy minimization

In the molecule, some clashes may appear during assembling of the system. The energy minimization step is necessary to remove bad contacts. Without minimization, the energy of contacting atoms may be high enough to crash the simulation. During minimization, atoms will be moved to find the closest structure with acceptable energy.

For minimization process we use pmemd.cuda program of AMBER. The input file `min.in` for this step provided in `Dickerson_Drew_Dodecamer_files/` consists of the lines presented below:

```
Minimize
 &cntrl
  imin=1,
  ntx=1,
  igb=8,
  irest=0,
  maxcyc=2000,
  ncyc=1000,
  ntpr=100,
  ntwx=0,
  ntr=1,
  cut=9999.0,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
RESTRAIN DNA
0.1
RES 1 24
END
END
```

`imin=1`  – turn minimization regime on

`ntx=1` – read coordinates from input coordinates file

`igb=8` – specify implicit solvent model GBneck2

`irest=0` – ignore input velocities

`maxcyc=2000` – limit of minimization cycles

`ncyc=1000`  – the number of minimization cycles with the steepest descent algorithm applied. The conjugate gradient algorithm is used for another 1000 steps

`ntpr=100` – write down output every 100 cycles

`ntwx=0` – do not write coordinate trajectory file

`ntr=1` – apply restraints of the group of atoms specified below to the reference coordinates

`cut=9999.0` – Cutoff distance of nonbonded interaction calculation in angstroms. The higher the number the more interacting atoms are considered and the more accurate and computationally expensive the calculaition is. For implicit solvent simulation huge number is usually used.

`gbion=3` – turn on GBION model

`nmropt=1` – turn on distance restraints for ions

`intdiel=1` – internal dielectric of the solute molecule

`gbsa=3` – take into account the energy of the surface tension

**Parameters implemented into GB approximation of interaction energy of different atom pairs:**

`gi_coef_1_p=1,` – $K_{GB}$ for pair solute atom – cation

`gi_coef_1_n=0.05,` – $K_{GB}$ for pair solute atom – anion

`gi_coef_2_pp=1,` – $K_{GB}$ for pair cation – cation

`gi_coef_2_pn=0.05,` – $K_{GB}$ for pair cation – anion

`gi_coef_2_nn=1,` – $K_{GB}$ for pair anion – anion

`intdiel_ion_1_p=54,` – $K_{\epsilon}$ for pair solute atom – cation

`intdiel_ion_1_n=10,` – $K_{\epsilon}$ for pair solute atom – anion

`intdiel_ion_2_pp=54,` – $K_{\epsilon}$ for pair cation – cation

`intdiel_ion_2_pn=10,` – $K_{\epsilon}$ for pair anion – cation

`intdiel_ion_2_nn=10,` – $K_{\epsilon}$ for pair anion – anion

`gb_neckscale_ion_1_p=1,` – $K_{NS}$ for pair solute atom – cation

`gb_neckscale_ion_1_n=1,` – $K_{NS}$ for pair solute atom – anion

`gb_neckscale_ion_2_pp=1,` – $K_{NS}$ for pair cation – cation

`gb_neckscale_ion_2_pn=1,` – $K_{NS}$ for pair anion – cation

`gb_neckscale_ion_2_nn=1,` – $K_{NS}$ for pair anion – anion

**Parameters of restraints**

`&wt type='END'` – no conditions are varied during the simulation

`DISANG=disang_NaCl.txt` – read restraints for ions from file

`RESTRAIN DNA` – specifying restraints for DNA

`0.1` – restraint constant for DNA

`RES 1 24` – specifying residues included into nucleosome, restraints will be applied to these residues.

To run the energy minimization change working directory to `Dickerson_Drew_Dodecamer_files` and type in command line:

`pmemd.cuda -O -i min.in -o min.out -p dna.top -c dna.crd -r min.ncrst -inf min.mdinfo -ref dna.crd`

Here flag `-O` induces overwriting the output files, `-i min.in` specifies file with input parameters, `-o min.out` specifies file, where output values will be written, `-p dna.top` specifies topology file, `-c dna.crd` – file with initial coordinates, `-r min.ncrst` – file with final coordinates, `-inf min.mdinfo` – file with intermediate values of energies and performance metrics, `-ref dna.crd` – reference coordinates for restraints for DNA atoms.

The output file of the simulation (`min.out`) should look like this:
```
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Mon May 23 02:33:10 2022
| Run on 01/05/2024 at 20:19:55

|   Executable path: pmemd
| Working directory: /home/YOUR_DIRECTORY
|          Hostname: strugatsky.cbb.lan

  [-O]verwriting output
```

and so on. If the simulation goes as should, there will be a section with results of the simulation:

```
--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------
  
  
  
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -5.3174E+03     1.4958E+01     8.2467E+01     C5        204

 BOND    =       83.8034  ANGLE   =      189.4818  DIHED      =      613.8166
 VDWAALS =     -374.4397  EEL     =     3356.1615  EGB        =    -6438.5121
 1-4 VDW =      254.7749  1-4 EEL =    -3002.4574  RESTRAINT  =        0.0000
 NMR restraints: Bond =    0.000   Angle =     0.000   Torsion =     0.000
===============================================================================
```

This goes on upto NSTEP of 2000.

After the simulation files `min.out`, `min.ncrst` should appear, the last one will be used as a starting point for further simulations.
### 2.3 Heating
In this step the system will be heated from 0 K to 300 K linearly. The input file for this step `heat.in` includes the lines below:

```
Heat
 &cntrl
  imin=0,
  igb=8,
  ntx=1,
  irest=0,
  nstlim=20000,
  dt=0.002,
  ntf=2,
  ntc=2,
  tempi=0.0,
  temp0=300.0,
  ntpr=100,
  ntwx=100,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  ntr=1,
  gamma_ln=1,
  nmropt=1,
  ig=-1,
  gbion=3,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
&wt type='TEMP0', istep1=0, istep2=19999, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=19999, istep2=20000, value1=300.0, value2=300.0 /
&wt type='END' /
DISANG=disang_NaCl.txt
&end
RESTRAIN DNA
0.1
RES 1 24
END
END
```

The parameters of the simulation are:

`imin=0` – MD simulation without minimization

`nstlim=20000` – length of the simulation in time steps

`dt=0.002` – time step of simulation in ps

`ntf=2` – turning calculation of the force for SHAKE constrained bonds of

`ntc=2` – Enable SHAKE to constrain all bonds involving hydrogen

`tempi=0.0` – initial temperature of the system in K

`temp0=300.0` – final temperature of the system in K

`ntpr=100` – write values to out file every 100 steps

`ntwx=100` – add snapshot to trajectory file every 100 steps

`ntb=0` – no periodic boundary conditions

`ntp=0` – turning off barostat

`ntt=3` – turning on Langevin thermostat

`gamma_ln=1` – Langevin thermostat collision frequency. In case of implicit water also controls speed of atoms

`ig=-1`– random seed for Langevin dynamics

`&wt type='TEMP0', istep1=0, istep2=19999, value1=0.0, value2=300.0` – defining heating of the system from 0 K to 300 K

To run the heating of the system type in command line being in the directory  `Dickerson_Drew_Dodecamer_files`:

`pmemd.cuda -O -i heat.in -o heat.out -p dna.top -c min.ncrst -ref dna.crd -r heat.ncrst -x heat.nc -inf heat.mdinfo`

AMBER would produce the next files after this simulation: 
* `heat.out` – contains aggregate parameters of the system, such as temperature and energies
* `heat.nc` – contains time series of the atom coordinates of the system
* `heat.ncrst` – final coordinates of the system

### 2.4 Equilibration of the system

Now let ions relax around DNA at 300 K. The input file for this step `equil.in` consists of the lines below:

```
equilibration
 &cntrl
  imin=0,
  igb=8,
  ntx=5,
  irest=1,
  nstlim=1500000,
  dt=0.002,
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=1000,
  ntwx=1000,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  gamma_ln=0.05,
  ig=-1,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

`irest=1`, `ntx=5`: read coordinates+velocities from `heat.ncrst`.

`nstlim=1 500 000`, `dt=0.002`: 3 ns total.

No `RESTRAIN DNA` (unless you wish to restrain DNA lightly; here we allow DNA to sample freely).

Ions remain constrained by `DISANG`.

To run the simulation, type in the command line being in  `Dickerson_Drew_Dodecamer_files`: 

`pmemd.cuda -O -i equil.in -o equil.out -p dna.top -c heat.ncrst -r equil.ncrst -x equil.nc -inf equil.mdinfo -ref dna.crd`

Outputs:  
• `equil.out`  
• `equil.nc` (3 ns trajectory)  
• `equil.ncrst` (final coordinates+velocities)

### 2.5 Production run

After equilibration, run production MD to sample DNA conformation. See in `Dickerson_Drew_Dodecamer_files` or create file `prod.in`:

```
Production
 &cntrl
  imin=0,
  igb=8,
  ntx=5,
  irest=1,
  nstlim=5000000,
  dt=0.002,
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=500,
  ntwx=500,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  gamma_ln=0.05,
  ig=-1,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

* `nstlim=5 000 000`, `dt=0.002`: 10 ns production run.
* `ntpr=500`, `ntwx=500`: print/write output once per 1 ps.
* DNA is unrestrained; ions remain semi-restrained by `DISANG`.

Run production simulation:

`pmemd.cuda -O -i prod.in -o prod.out -p dna.top -c equil.ncrst -r prod.ncrst -x prod.nc -inf prod.mdinfo -ref dna.crd`

Outputs:  
* `prod.out` (energies)
* `prod.nc` (trajectory, every 1 ps)
* `prod.ncrst` (final snapshot)

### 2.6 Analysis of DNA stability

Use **CPPTRAJ** to compute RMSD of DNA heavy atoms over the production trajectory. Input file `cpptraj.in` for this analysis should consist of the lines below:

```
trajin prod.nc
rms ToFirst :1-24&!@H= out rms_dna.txt
run
quit
```

To run the analysis, type:

`cpptraj -p dna.top < cpptraj.in`

To visualize the data, use python script:

```
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('rms_dna.txt', skiprows=2)
time = data[:,0]  # ps
rms  = data[:,1]  # Å

plt.figure(figsize=(10,6))
plt.plot(time, rms, linewidth=2)
plt.xlabel('Time (ps)', fontsize=14)
plt.ylabel('RMSD (Å)', fontsize=14)
plt.title('DNA RMSD Over 10 ns Production', fontsize=16)
plt.axhline(np.mean(rms), linestyle='--', color='gray')
plt.tight_layout()
plt.savefig('rmsd_dna.png', dpi=300)
```

To run the provided script, type:

`python graph.py`

The file `rmsd_dna.png` would contain the graph of RMSD vs. time. It should look like this: ![rmsd_dna](Pictures/rmsd_dna.png)
### 2.7 Visualization of trajectory using ChimeraX

1. **Open** ChimeraX.
2. Drag and drop file **dna.pdb**
3. **Trajectory**: Drag and drop **prod.nc** (MD trajectory), choose`Amber netCDF coordinates` in pop-up window.
4. Click **Play** to observe DNA + ions’ dynamics.
You should observe something like this:![Prod_frame_DNA](Pictures/Prod_frame_DNA.png)


## 3. SImulation of a nucleosome using GBION model

### 3.1 Preparation of the system for simulations

#### 3.1.1 Building topology + placing ions using tleap
It is required to use `Nucleosome.pdb` file as a structure simulated. It can be found in directory `Nucleosome_simulation_files/Prep`. In the same directory see a file (or create your own) named `tleap.script` with:
```
source leaprc.DNA.OL15
source leaprc.protein.ff19SB
loadoff atomic_ions.lib
source leaprc.water.opc
loadAmberParams frcmod.ionsjc_tip4pew
set default PBradii mbondi3
mol = loadpdb Nucleosome.pdb
addions mol K+ 5256
addions mol Cl- 5119
saveamberparm mol nucleosome.top nucleosome.crd
savepdb mol nucleosome.pdb
quit
```

`source leaprc.DNA.OL15` – load the OL15 DNA force field.

`loadoff atomic_ions.lib`, `source leaprc.water.opc` & `loadAmberParams frcmod.ionsjc_tip4pew`: load ion parameters.

`set default PBradii mbondi3` – use mbondi3 radii (suitable for implicit solvent we use).

`mol = loadpdb 1bna.pdb` – loading structure of DNA.

`addions mol K+ 5256` & `addions mol Cl- 5119`: add 5256 potassium and 5119 chloride ions to neutralize and mimic roughly 150 mM.

`saveamberparm mol nucleosome.top nucleosome.crd` – saving topology and initial coordinates of the structure for simulation.

`savepdb mol nucleosome.pdb` – saving the system to PDB file.

To run the script, change working directory to `Nucleosome_simulation_files/Prep`and type in command line:

`tleap -f tleap.script`

Output:  
* `nucleosome.top` (AMBER topology)
* `nucleosome.crd` (initial coordinates)
* `nucleosome.pdb` (PDB file suitable for visualization of the initial state of the system)
#### 3.1.2 Generating distance-restraints for ions (disang.py)

See the general rationale under distance restraints in previous section. In the same `Nucleosome_simulation_files/Prep/` directory type:

`python disang.py`

`Molecule_file = 'Nucleosome.pdb'` – specifies PDB-file with simulated biomolecule

`namedisang = "disang_KCl.txt"` – specifies output restraints file

`number_of_ions= 10375` – specifies number of ions to restrain

`starting_index_of_ion= 24942` – atom index of the first ion to restrain.

**line 50:** `        disang.write("  iat=-1,-1,r1=0.0,r2=0.0,r3=240,r4=250,rk2=0.0, rk3=20.0,\n")` – r3, r4 and rk3 specify parameters explained earlier.

This produces `disang_KCl.txt`, containing blocks like:
```
&rst
iresid=0,
  iat=-1,-1,r1=0.0,r2=0.0,r3=240,r4=250,rk2=0.0, rk3=20.0,
igr1=5446,5445,5441,3752,3750,3748,3747,3746,3730,2169,igr2=35316
 /
```

- `iresid=0`: atom indices are specified directly.
- `iat=-1,-1`: instructs AMBER to read `igr1` and `igr2` groups, then restrain their centers of mass.  
    • `igr1`: (set of 10 atoms) → DNA’s geometric center.  
    • `igr2`: index of an ion.
- `r1=0.0, r2=0.0, r3=40, r4=50`: define a flat-bottom restraint (0–240 Å flat, then a parabola up to 250 Å).
- `rk3=20.0`: force constant (20 kcal · mol<sup>−1</sup> · Å<sup>−2</sup>) as atoms go from 240→250 Å.

Copy `nucleosome.*` and `disang_NaCl.txt` up one level:

`cp nucleosome.* ../`

`cp disang_KCl.txt ../`

If you prefer, run the provided script:

`bash prep.sh`

This will auto-generate `dna.top`,`dna.crd`,`dna.pdb` and `disang_NaCl.txt` and copy up.

### 3.2 Energy minimization

In the molecule, some clashes may appear during assembling of the system. The energy minimization step is necessary to remove bad contacts. Without minimization, the energy of contacting atoms may be high enough to crash the simulation. During minimization, atoms will be moved to find the closest structure with acceptable energy.

For minimization process we use pmemd.cuda program of AMBER. The input file `min.in` for this step provided in `Nucleosome_simulation_files/` consists of the lines presented below:

```
Minimize
 &cntrl
  imin=1,
  ntx=1,
  igb=8,
  irest=0,
  maxcyc=2000,
  ncyc=1000,
  ntpr=100,
  ntwx=0,
  ntr=1,
  cut=9999.0,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
RESTRAIN DNA
20.0
RES 1 1264
END
END
```

`imin=1`  – turn minimization regime on

`ntx=1` – read coordinates from input coordinates file

`igb=8` – specify implicit solvent model GBneck2

`irest=0` – ignore input velocities

`maxcyc=2000` – limit of minimization cycles

`ncyc=1000`  – the number of minimization cycles with the steepest descent algorithm applied. The conjugate gradient algorithm is used for another 1000 steps

`ntpr=100` – write down output every 100 cycles

`ntwx=0` – do not write coordinate trajectory file

`ntr=1` – apply restraints of the group of atoms specified below to the reference coordinates

`cut=9999.0` – Cutoff distance of nonbonded interaction calculation in angstroms. The higher the number the more interacting atoms are considered and the more accurate and computationally expensive the calculaition is. For implicit solvent simulation huge number is usually used.

`gbion=3` – turn on GBION model

`nmropt=1` – turn on distance restraints for ions

`intdiel=1` – internal dielectric of the solute molecule

`gbsa=3` – take into account the energy of the surface tension

**Parameters implemented into GB approximation of interaction energy of different atom pairs:**

`gi_coef_1_p=1,` – $K_{GB}$ for pair solute atom – cation

`gi_coef_1_n=0.05,` – $K_{GB}$ for pair solute atom – anion

`gi_coef_2_pp=1,` – $K_{GB}$ for pair cation – cation

`gi_coef_2_pn=0.05,` – $K_{GB}$ for pair cation – anion

`gi_coef_2_nn=1,` – $K_{GB}$ for pair anion – anion

`intdiel_ion_1_p=36,` – $K_{\epsilon}$ for pair solute atom – cation

`intdiel_ion_1_n=10,` – $K_{\epsilon}$ for pair solute atom – anion

`intdiel_ion_2_pp=36,` – $K_{\epsilon}$ for pair cation – cation

`intdiel_ion_2_pn=10,` – $K_{\epsilon}$ for pair anion – cation

`intdiel_ion_2_nn=10,` – $K_{\epsilon}$ for pair anion – anion

`gb_neckscale_ion_1_p=1,` – $K_{NS}$ for pair solute atom – cation

`gb_neckscale_ion_1_n=1,` – $K_{NS}$ for pair solute atom – anion

`gb_neckscale_ion_2_pp=1,` – $K_{NS}$ for pair cation – cation

`gb_neckscale_ion_2_pn=1,` – $K_{NS}$ for pair anion – cation

`gb_neckscale_ion_2_nn=1,` – $K_{NS}$ for pair anion – anion

**Parameters of restraints**

`&wt type='END'` – no conditions are varied during the simulation

`DISANG=disang_KCl.txt` – read restraints for ions from file

`RESTRAIN DNA` – specifying restraints for DNA

`0.1` – restraint constant for DNA

`RES 1 1264` – specifying residues included into nucleosome, restraints will be applied to these residues.

To run the energy minimization change working directory to `Nucleosome_simulation_files` and type in command line:

`pmemd.cuda -O -i min.in -o min.out -p dna.top -c dna.crd -r min.ncrst -inf min.mdinfo -ref dna.crd`

Here flag `-O` induces overwriting the output files, `-i min.in` specifies file with input parameters, `-o min.out` specifies file, where output values will be written, `-p dna.top` specifies topology file, `-c dna.crd` – file with initial coordinates, `-r min.ncrst` – file with final coordinates, `-inf min.mdinfo` – file with intermediate values of energies and performance metrics, `-ref dna.crd` – reference coordinates for restraints for DNA atoms.

The output file of the simulation (`min.out`) should look like this:
```
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Mon May 23 02:33:10 2022
| Run on 01/05/2024 at 20:19:55

|   Executable path: pmemd
| Working directory: /home/YOUR_DIRECTORY
|          Hostname: strugatsky.cbb.lan

  [-O]verwriting output
```

and so on. If the simulation goes as should, there will be a section with results of the simulation:

```
--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------
  
  
  
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1       3.5040E+07     3.8015E+02     2.7334E+04     OP1     19435
  
 BOND    =    18843.4733  ANGLE   =    13343.7741  DIHED      =    13417.4807
 UB      =        0.0000  IMP     =        0.0000  CMAP       =      970.3763
 VDWAALS = 34999554.9114  EEL     = -1636975.1804  EGB        =  1613953.5606
 1-4 VDW =    11240.2218  1-4 EEL =    -1315.4912  RESTRAINT  =     4631.7865
 ESURF   =     2124.0188
 EAMBER  = 35035157.1453
 NMR restraints: Bond = 4631.787   Angle =     0.000   Torsion =     0.000
===============================================================================
```

This goes on upto NSTEP of 2000.

After the simulation files `min.out`, `min.ncrst` should appear, the last one will be used as a starting point for further simulations.
### 3.3 Heating
In this step the system will be heated from 0 K to 300 K linearly. The input file for this step `heat.in` includes the lines below:

```
Heat
 &cntrl
  imin=0,
  igb=8,
  ntx=1,
  irest=0,
  nstlim=60000,
  dt=0.002,
  ntf=2,
  ntc=2,
  tempi=0.0,
  temp0=300.0,
  ntpr=100,
  ntwx=100,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  ntr=1,
  gamma_ln=1,
  nmropt=1,
  ig=-1,
  gbion=3,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
&wt type='TEMP0', istep1=0, istep2=59999, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=59999, istep2=60000, value1=300.0, value2=300.0 /
&wt type='END' /
DISANG=disang_KCl.txt
&end
RESTRAIN CORE
0.1
RES 40 133
END
RESTRAIN CORE
0.1
RES 161 237
END
RESTRAIN CORE
0.1
RES 254 354
END
RESTRAIN CORE
0.1
RES 398 485
END
RESTRAIN CORE
0.1
RES 527 620
END
RESTRAIN CORE
0.1
RES 648 724
END
RESTRAIN CORE
0.1
RES 741 841
END
RESTRAIN CORE
0.1
RES 885 972
END
END
```

The parameters of the simulation are:

`imin=0` – MD simulation without minimization

`nstlim=60000` – length of the simulation in time steps

`dt=0.002` – time step of simulation in ps

`ntf=2` – turning calculation of the force for SHAKE constrained bonds of

`ntc=2` – Enable SHAKE to constrain all bonds involving hydrogen

`tempi=0.0` – initial temperature of the system in K

`temp0=300.0` – final temperature of the system in K

`ntpr=100` – write values to out file every 100 steps

`ntwx=100` – add snapshot to trajectory file every 100 steps

`ntb=0` – no periodic boundary conditions

`ntp=0` – turning off barostat

`ntt=3` – turning on Langevin thermostat

`gamma_ln=1` – Langevin thermostat collision frequency. In case of implicit water also controls speed of atoms

`ig=-1`– random seed for Langevin dynamics

`&wt type='TEMP0', istep1=0, istep2=59999, value1=0.0, value2=300.0` – defining heating of the system from 0 K to 300 K

To run the heating of the system type in command line being in the directory  `Nucleosome_simulation_files`:

`pmemd.cuda -O -i heat.in -o heat.out -p dna.top -c min.ncrst -ref dna.crd -r heat.ncrst -x heat.nc -inf heat.mdinfo`

AMBER would produce the next files after this simulation: 
* `heat.out` – contains aggregate parameters of the system, such as temperature and energies
* `heat.nc` – contains time series of the atom coordinates of the system
* `heat.ncrst` – final coordinates of the system

### 3.4 Equilibration of the system

Now let ions relax around DNA at 300 K. The input file for this step `equil.in` consists of the lines below:

```
equilibration
 &cntrl
  imin=0,
  igb=8,
  ntx=5,
  irest=1,
  nstlim=1500000,
  dt=0.002,
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=1000,
  ntwx=1000,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  gamma_ln=0.05,
  ig=-1,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

`irest=1`, `ntx=5`: read coordinates+velocities from `heat.ncrst`.

`nstlim=1 500 000`, `dt=0.002`: 1.5 ns total.

No `RESTRAIN DNA` (unless you wish to restrain DNA lightly; here we allow DNA to sample freely).

Ions remain constrained by `DISANG`.

To run the simulation, type in the command line being in  `Dickerson_Drew_Dodecamer_files`: 

`pmemd.cuda -O -i equil.in -o equil.out -p dna.top -c heat.ncrst -r equil.ncrst -x equil.nc -inf equil.mdinfo -ref dna.crd`

Outputs:  
• `equil.out`  
• `equil.nc` (3 ns trajectory)  
• `equil.ncrst` (final coordinates+velocities)

### 3.5 Production run

After equilibration, run production MD to sample DNA conformation. See in `Nucleosome_simulation_files` or create file `prod.in`:

```
Production
 &cntrl
  imin=0,
  igb=8,
  ntx=5,
  irest=1,
  nstlim=5000000,
  dt=0.002,
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=500,
  ntwx=500,
  cut=9999.0,
  ntb=0,
  ntp=0,
  ntt=3,
  gamma_ln=0.05,
  ig=-1,
  gbion=3,
  nmropt=1,
  intdiel=1,
  gbsa=3,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=54,
  intdiel_ion_1_n=10,
  intdiel_ion_2_pp=54,
  intdiel_ion_2_pn=10,
  intdiel_ion_2_nn=10,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

* `nstlim=5 000 000`, `dt=0.002`: 10 ns production run.
* `ntpr=500`, `ntwx=500`: print/write output once per 1 ps.
* DNA is unrestrained; ions remain semi-restrained by `DISANG`.

Run production simulation:

`pmemd.cuda -O -i prod.in -o prod.out -p dna.top -c equil.ncrst -r prod.ncrst -x prod.nc -inf prod.mdinfo -ref dna.crd`

Outputs:  
* `prod.out` (energies)
* `prod.nc` (trajectory, every 1 ps)
* `prod.ncrst` (final snapshot)

### 3.6 Analysis of DNA stability

Use **CPPTRAJ** to compute RMSD of DNA heavy atoms over the production trajectory. Input file `cpptraj.in` for this analysis should consist of the lines below:

```
trajin prod.nc
reference nucleosome.pdb
autoimage
rms :40-135,161-237,254-354,398-485,527-622,648-724,741-841,885-974&@CA,N,C nofit reference out rmsd_core.txt
run
quit
```

To run the analysis, type:

`cpptraj -p nucleosome.top < cpptraj.in`

To visualize the data, use python script:

```
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('rms_core.txt', skiprows=2)
time = data[:,0]  # ps
rms  = data[:,1]  # Å

plt.figure(figsize=(10,6))
plt.plot(time, rms, linewidth=2)
plt.xlabel('Time (ps)', fontsize=14)
plt.ylabel('RMSD (Å)', fontsize=14)
plt.title('Histone core RMSD Over 10 ns Production', fontsize=16)
plt.axhline(np.mean(rms), linestyle='--', color='gray')
plt.tight_layout()
plt.savefig('rmsd_core.png', dpi=300)
```

To run the provided script, type:

`python graph.py`

The file `rmsd_core.png` would contain the graph of RMSD vs. time. It should look like this:
![rmsd_core](Pictures/rmsd_core.png)
### 3.7 Visualization of trajectory using ChimeraX

1. **Open** ChimeraX.
2. Drag and drop file **nucleosome.pdb**
3. **Trajectory**: Drag and drop **prod.nc** (MD trajectory).
4. Click **Play** to observe nucleosome + ions’ dynamics.
You can apply the ChimeraX commands below for better visual representation of histone tails and nucleosome core:
```
select :40-135,161-237,254-354,398-485,527-622,648-724,741-841,885-974,975-1264
hide sel cartoon
hide sel atoms
show sel surface
set bgColor white
select :40-135,161-237,254-354,398-485,527-622,648-724,741-841,885-974
color sel dark grey
select :975-1264
color sel #f68419
select :1-39,136-160,238-253,355-365,366-397,488-526,623-647,725-740,842-852,853-884
show sel car
color sel #6a137a
car style t 5
select clear
```
You should observe something like this:![Prod_frame_nucleosome](Pictures/Prod_frame_nucleosome.png)

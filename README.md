## Simulation of large systems using implicit solvent/explicit ions model

This tutorial aims to demonstrate simulation of a large system containing DNA (nucleosome) in implicit solvent in combination with explicit ions. We assume that user has basic skills of running MD simulations using AMBER package. It also requires installing python with NumPy and Matplotlib libraries and CHIMERA([link for downloading](https://www.cgl.ucsf.edu/chimera/download.html)).

The outline of this tutorial:

1. Preparing the system with nucleosome for the simulations
2. Preparing input files for the simulation of nucleosome
3. Simulation of the nucleosome using implicit water/explicit ions model
4. Analysis of the nucleosome stability
5. Visualization of ion distribution around DNA using CHIMERA

### SECTION 1. Brief introduction
Large systems, such as nucleosomes, require 5- or 6-digits numbers of water molecules for the proper simulation with explicit water model. Another approach to treating water (implicit water model) allows taking water into account as a continuous environment around the molecule. The solvation free energy of the solvated molecules defines their behavior in the solution. The solvation energy $\Delta G_{solv}$ can be represented as a sum of electrostatic $\Delta G_{el}$ and non-polar $\Delta G_{np}$ contributions:

$\Delta G_{solv} = \Delta G_{el} + \Delta G_{np}$ 

which are estimated independently in most cases. One of the most popular approximation to calculating $\Delta G_{np}$ is based on the assumption that $\Delta G_{np}$ is proportional to the solvent-accessible surface area (SASA).

One of the most widely used approximation for calculation $\Delta G_{el}$ is the generalized Born (GB) model. The canonical GB approximation is based on the equation originally proposed by Still *et al.*:

$\Delta G_{el} = -\frac{1}{2} \left(\frac{1}{\epsilon_{in}} -\frac{1}{\epsilon_{out}}\right)\sum_{i,j} \frac{q_{i}q_{j}}{\sqrt{d^2_{ij} + R_{i}R_{j}\exp\left(-d^2_{ij}/(4 R_i R_j)\right)}}$,

where $\epsilon_{in}$ and $\epsilon_{out}$ are the dielectric constants of the solute and the solvent, respectively, $d_{ij}$ is the distance between solute atoms $i$ and $j$, and $q_{i}$ are the atomic charges. The key parameters modulating the interaction energy are *the effective Born radii* $R_i$. $R_i^{-1}$ characterizes the average degree of solvent exposure of atom $i$.

The GB model describes behaviour of the topologically connected structures really well, but it does not take into account descrete ions around solute. The ISEXI model was developed as an extension of GB model to simulate ions around DNA with implicit solvent. It implies additional coefficients to the GB equation:

$\Delta G_{el} = -\frac{1}{2} \left(\frac{1}{\mathbf{\epsilon_{in}(a,b)}} -\frac{1}{\epsilon_{out}}\right)\sum_{i,j} \frac{q_{i}q_{j}}{\sqrt{d^2_{ij} + R_{i}R_{j}\exp\left(-d^2_{ij}/(\mathbf {\gamma(a,b)} R_i R_j)\right)}}$

The expression above emphasizes the main idea of the ISEXI model that the functional form of charge-charge interaction is different for charges that are connected through the solute or the solvent. This is achieved by variations of $\mathbf{\epsilon_{in}}$ and $\mathbf{\gamma(a,b)}$  for different pairs of interacting atoms separately: solute-solute, solute-ion and ion-ion.

### SECTION 2. Preparation of the system for simulations

#### 1. Preparation of the nucleosome structure

Here we use the previously developed protocol of nucleosome simulation, reworked slightly to be suitable for simulations using implicit solvent/explicit ions model. To get the initial structures for simulations, go to the site ([link](https://zenodo.org/records/8315307)) and download archive. You can find the nucleosome structure via path: `md_setup/md_protocol/OPC/ff99SB/R3A/01_equil_histone_tails/1_build/nucleosome.pdb`, also available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/nucleosome.pdb)). You can visualize the structure using CHIMERA, just open the `nucleosome.pdb` file using the program. You should see the picture like below:

![(./Pictures/nucleosome_stretched.png)](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/nucleosome_stretched.png)

For structure preparation, we use the leap script (file tleap.script, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/tleap.script))).

```
source leaprc.DNA.OL15
source leaprc.protein.ff14SB
set default PBradii mbondi3
loadAmberParams frcmod.ionsjc_tip4pew
source leaprc.water.opc
mol = loadpdb nucleosome.pdb
addions mol Na+ 5223
addions mol Cl- 5093
saveamberparm mol dna.top dna.crd
savepdb mol dna.pdb
quit
```

details of the script:

`source leaprc.DNA.OL15` – loading the force field for DNA

`source leaprc.protein.ff14SB` – loading force field for histone proteins

`set default PBradii mbondi3` – setting the default PBradii to mbondi3 – this fits implicit water model we use here

`loadAmberParams frcmod.ionsjc_tip4pew`  – loading ion parameters, that were used for optimization of implicit water model parameters

`source leaprc.water.opc `– loading the library that contains ions

`mol = loadpdb nucleosome.pdb` – loading structure of nucleosome

`addions mol Na+ 5223` – adding cations to the nucleosome structure

`addions mol Cl- 5093` – adding anions to the nucleosome structure

`saveamberparm mol nuc.top nuc.crd` – saving topology and initial coordinates of the structure for simulation

`savepdb mol nuc.pdb` – saving the system to PDB file

To run the script, type in command line:

`tleap -f tleap.script`

After running the script you should get 3 new files in the working directory: nuc.top nuc.crd and nuc.pdb. First one is a topology file of the system you created, it contains parameters required for MD simulation that define energies of all the interactions, second one contains the initial coordinates of the system. You can visualize the result using CHIMERA. Just open the `dna.pdb` file using the program. You should see a picture like this:
![./Pictures/nucleosome_stretched_ions_all.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/nucleosome_stretched_ions_all.png) ![./Pictures/nucleosome_stratched_ions_part.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/nucleosome_stretched_ions_part.png)
Both are the fugures of the same system, at the second some ions are hidden for clarity.

#### 2. Preparing file with restraints

Since the system is simulated in implicit water there is no periodic boundary conditions for the system. Thus, ions will diffuse away from the nucleosome, being not anchored to the molecule simulated. Here we use standard practice implemented into AMBER package, defining the restraints that are applied to a pair of atoms. Use the file `disang.py` (available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/disang.py))) to prepare the restraints.

To run the generation of the file, type in command line:

`python disang.py`

After this you will get a file with restraints. It contains of a number of the repetitive lines:
```
&rst
iresid=0,
iat=-1,-1,r1=0.0,r2=0.0,r3=240,r4=250,rk2=0.0, rk3=20.0, igr1=5444,5443,5439,3751,3749,3747,3746,3745,2169,2168,igr2=24946
/
```

Here flag `iresid=0` allows to select individual atoms in the molecule. Flags `iat=-1,-1` make the program read groups of files `igr1` and `igr2` and calculate restraints between their centers of mass. The first group of atoms are 10 closest ones to the center of mass of the whole nucleosome. Distances `r1`, `r2`, `r3` and `r4` define the restraint force graph form:
![./Pictures/restraints.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/restraints.png)
From 0 to `r1` force linearly depends on the distance, from `r1` to `r2` parabolically, `r2-r3` is a flat region, from `r3` to `r4`  – parabolically, from `r4` to $\infty$  – linearly. In our case 0-`r3` is a flat region:
![restraints_flat.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/restraints_flat.png)

### SECTION 3. Energy minimization

In the molecule, some clashes may appear during assembling of the system. The energy minimization step is necessary to remove bad contacts. Without minimization, the energy of contacting atoms may be high enough to crash the simulation. During minimization, atoms will be moved to find the closest structure with acceptable energy.

#### 1. Create input file for minimization step

For minimization process we use pmemd program of AMBER. The input file (`min.in`, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/min.in))) for this step consists of the lines presented below.

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
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=8,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=8,
  intdiel_ion_2_nn=8,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
  gbsa=3,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
RESTRAIN NUCLEOSOME
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

`cut=9999.0` – Cutoff distance of nonbonded interaction calculation in angstroms. The higher the number the more interacting atoms are considered and the more accurate and computationally expensive the calculaition is.

`gbion=3` – turn on ISEXI model

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

`intdiel_ion_1_n=8,` – $K_{\epsilon}$ for pair solute atom – anion

`intdiel_ion_2_pp=36,` – $K_{\epsilon}$ for pair cation – cation

`intdiel_ion_2_pn=8,` – $K_{\epsilon}$ for pair anion – cation

`intdiel_ion_2_nn=8,` – $K_{\epsilon}$ for pair anion – anion

`gb_neckscale_ion_1_p=1,` – $K_{NS}$ for pair solute atom – cation

`gb_neckscale_ion_1_n=1,` – $K_{NS}$ for pair solute atom – anion

`gb_neckscale_ion_2_pp=1,` – $K_{NS}$ for pair cation – cation

`gb_neckscale_ion_2_pn=1,` – $K_{NS}$ for pair anion – cation

`gb_neckscale_ion_2_nn=1,` – $K_{NS}$ for pair anion – anion

**Parameters of restraints**

`&wt type='END'` – no conditions are varied during the simulation

`DISANG=disang_NaCl.txt` – read restraints for ions from file

`RESTRAIN NUCLEOSOME` – specifying restraints for nucleosome

`20.0` – restraint constant for nucleosome

`RES 1 1264` – specifying residues included into nucleosome, restraints will be applied to these restraints.

#### 2. Run the energy minimization

To run the energy minimization type in command line:

`pmemd.cuda -O -i min.in -o min.out -p dna.top -c dna.crd -r min.ncrst -inf min.mdinfo -ref dna.crd`

Here flag `-O` induces overwriting the output files, `-i min.in` specifies file with input parameters, `-o min.out` specifies file, where output values will be written, `-p dna.top` specifies topology file, `-c dna.crd` – file with initial coordinates, `-r min.ncrst` – file with final coordinates, `-inf min.mdinfo` – file with intermediate values of energies, `-ref dna.crd` – reference coordinates for restraints for nucleosome atoms.

The output file of the simulation (`min.out`) should look like this:
```
          -------------------------------------------------------
          Amber 22 PMEMD                              2022
          -------------------------------------------------------

| PMEMD implementation of SANDER, Release 22

|  Compiled date/time: Mon May 23 02:33:10 2022
| Run on 01/05/2024 at 20:19:55

|   Executable path: pmemd
| Working directory: /data/kolesnikov/E_0677_nucleosome_conses_GBions/GBion_K_gbsa
|          Hostname: strugatsky.cbb.lan

  [-O]verwriting output
```
and so on. If the simulation goes as should, there will be a section with results of the simulation:
```
--------------------------------------------------------------------------------
   4.  RESULTS
--------------------------------------------------------------------------------
  
  
  
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1       3.5000E+07     2.3299E+06     4.6750E+08     HE3       543
  
 BOND    =    18843.4733  ANGLE   =    13343.7741  DIHED      =    20653.3169
 VDWAALS = 34997534.8245  EEL     = -2426200.6972  EGB        =  2365964.3321
 1-4 VDW =    11240.4517  1-4 EEL =    -1315.4897  RESTRAINT  =        0.0000
 ESURF   =      265.0066
 NMR restraints: Bond =    0.000   Angle =     0.000   Torsion =     0.000
===============================================================================
```
This goes on upto NSTEP of 2000.

### SECTION 4. Heating

In this step the system will be heated from 0 K to 300 K linearly. The input file for this step (`heat.in`, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/heat.in))) consists of the lines below:
```
Heat
 &cntrl
  imin=0,
  igb=8,
  ntx=1,
  irest=0,
  nstlim=10000,
  dt=0.0002,
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
  gamma_ln=0.05,
  nmropt=1,
  ig=-1,
  gbion=3,
  intdiel=1,
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=8,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=8,
  intdiel_ion_2_nn=8,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
  gbsa=3,
 /
&wt type='TEMP0', istep1=0, istep2=9999, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=9999, istep2=10000, value1=300.0, value2=300.0 /
&wt type='END' /
DISANG=disang_NaCl.txt
&end
```
The parameters of the simulation are:

`imin=0` – MD simulation without minimization

`nstlim=10000` – length of the simulation in time steps

`dt=0.0002` – time step of simulation in ps

`ntf=2` – turning calculation of the force for SHAKE constrained bonds of

`ntc=2` – Enable SHAKE to constrain all bonds involving hydrogen

`tempi=0.0` – initial temperature of the system in K

`temp0=300.0` – final temperature of the system in K

`ntpr=100` – write values to out file every 100 steps

`ntwx=100` – add snapshot to trajectory file every 100 steps

`ntb=0` – no periodic boundary conditions

`ntp=0` – turning off barostat

`ntt=3` – turning on Langevin thermostat

`gamma_ln=0.05` – Langevin thermostat collision frequency. In case of implicit water also controls speed of atoms

`ig=-1`– random seed for Langevin dynamics

`&wt type='TEMP0', istep1=0, istep2=9999, value1=0.0, value2=300.0` – defining heating of the system from 0 K to 300 K

To run the heating of the system type in command line:

`pmemd.cuda -O -i heat.in -o heat.out -p dna.top -c min.ncrst -ref dna.crd -r heat.ncrst -x heat.nc -inf heat.mdinfo`

### SECTION 5. Equilibration of the histone tails

To simulate the nucleosome properly, one should start with equilibrating of the initially stretched histone tails. To do so, one should run MD simulation with initially stretched histone tails for long enough time to let the tails to condense to nucleosome. The input file for this step (`equil.in`, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/equil.in))) consists of the lines below:
```
equilibration
 &cntrl
  imin=0,
  igb=8,
  ntx=5,
  irest=1,
  nstlim=1500000,
  dt=0.001,
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
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=8,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=8,
  intdiel_ion_2_nn=8,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
  gbsa=3,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

To run the simulation, type in the command line: 

`pmemd.cuda -O -i equil.in -o equil.out -p dna.top -c heat.ncrst -r equil.ncrst -x equil.nc -inf equil.mdinfo -ref dna.crd`

### SECTION 6. Production run

The nucleosome with equilibrated tails may be can now. The input file for production run differs from the one for equilibrating of the histone tails by time of simulation and frequency of output. Parameters of productioon run (file `prod.in`, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/prod.in))) are listed below:
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
  gi_coef_1_p=1,
  gi_coef_1_n=0.05,
  gi_coef_2_pp=1,
  gi_coef_2_pn=0.05,
  gi_coef_2_nn=1,
  intdiel_ion_1_p=36,
  intdiel_ion_1_n=8,
  intdiel_ion_2_pp=36,
  intdiel_ion_2_pn=8,
  intdiel_ion_2_nn=8,
  gb_neckscale_ion_1_p=1,
  gb_neckscale_ion_1_n=1,
  gb_neckscale_ion_2_pp=1,
  gb_neckscale_ion_2_pn=1,
  gb_neckscale_ion_2_nn=1,
  gbsa=3,
 /
 &wt type='END',
 /
DISANG=disang_NaCl.txt
&end
```

To run the simulation, type in command line:

`pmemd.cuda -O -i prod.in -o prod.out -p dna.top -c equil.ncrst -r prod.ncrst -x prod.trj -inf prod.mdinfo -ref dna.crd`

### SECTION 7. Analysis and visualization of the results

#### 1. Analysis of RMSD

To analyse the stability of the nucleosome, we use RMSD analysis. CPPTRAJ program implemented into AMBER allows to calculate the RMSD of a molecule during the simulation. To run CPPTRAJ with input topology file of the system we simulate, type in command line:

`cpptraj -p dna.top`

You should see the output that looks like showed below:
  
```
CPPTRAJ: Trajectory Analysis. V5.1.0
    ___  ___  ___  ___
     | \/ | \/ | \/ | 
    _|_/\_|_/\_|_/\_|_
  
| Date/time: 01/11/24 15:30:39
| Available memory: 15.560 GB
  
Reading 'dna.top' as Amber Topology
Radius Set: ArgH and AspGluO modified Bondi2 radii (mbondi3)
Loading previous history from log 'cpptraj.log'
```

To load initial structure for reference type:

`trajin dna.crd`

Then to load the trajectory for further analysis type:

`trajin prod.trj

After loading the trajectory type in command line:

`rms ToFirst :40-133,161-237,254-354,398-485,527-620,648-724,741-841,885-972,975-1264 out rms_Nucleosome_no_tails.txt`

This command defines atoms that are taken into account in calculation. These are residues of the nucleosome, histone tails excluded.

To run the calculation you defined, type:

`run`

And to quit from the CPPTRAJ program, type:

`quit`

To visualize the results we have here, we use python libraries NumPy and Matplotlib. One can use any other desired method for graph visualization. The python script we use is (file `rmsd.py`, available via ([link](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Files/rmsd.py))) :

```
import numpy as np
import matplotlib.pyplot as plt
rmsd = np.loadtxt('rms_Nucleosome_no_tails.txt', skiprows=2)
plt.figure(figsize=(12,8))
plt.plot(rmsd[:,0]/1000, rmsd[:,1], label = 'Major', color='#6a137a', linewidth=2)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0,10)
plt.xlabel('Time (ns)', fontsize=30)
plt.ylabel('RMSD nucleosome ($\AA$)', fontsize=30)
plt.axhline(np.average(rmsd[:,1]), color='#2E7D32', linewidth=3)
plt.savefig('rmsd_nucleosome.png', bbox_inches='tight')
print(np.average(rmsd[:,1]))
```

To run the script, type in command line:

`python rmsd.py`

After running the script you should get the graph (file rmsd_nucleosome.png) similar to this one:

![./Pictures/rmsd_nucleosome.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/rmsd_nucleosome.png)

The green line emphasizes the average RMSD value. Values are not expected to be more than 5 Å.

#### 2. Visualization of the trajectory.

To visualize the trajectory itself we use program CHIMERA. Open it and choose in upper menu `Tools > MD/Ensemble Analysis > MD movie` For prmtop file choose `dna.top` and for trajectory `prod.trj`. You should see the picture like this:
![./Pictures/image_DNA_traj.png](https://github.com/EgorBiophys/ISEXI_tutorial/blob/main/Pictures/image_DNA_traj.png)

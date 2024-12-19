# Polymer-looping-under-tension-FJC-
Code used for publication: Polymer dynamics under tension: mean first passage time for looping

## LAMMPS
There are 2 main files: 1) main.lmp and 2) main_customised_fix.lmp. The former checks every 1000 timesteps whether the polymer has formed a loop and prints the timestep at which a loop is formed to the file bond_info.${seed}.txt. Due to this if statement however, this code might be slow, especially when one wants to check loop formation with a higher accurary than every 1000 timesteps. This issue is fixt in the latter main file, for which we have slightly changed the LAMMPS source code. We changed fix_bond_react.cpp as follows: https://github.com/wouterel/lammps/commit/314e9ca4e797eabb391b2d5b575cb6f5e259d4e0. With this fix, main_customised_fix.lmp prints the exact timestep a loop is formed, without requiring an if statement. 

The code is now written for a polymer containing 50 beads. When one wants to change this, the data file polymer50.data has to be modified, as well as the main file. In the main file, one should change the id of the right bead to the new number of beads. Furthermore, the code applies of force of 1 (LJ units) on both ends of the polymer, but this can also be changed to a different number. 

Finally, the capture radius is 1 (LJ units) in the code here. When one wants to change this, one has to adapt the third number in the first reaction step. For example, now it says react rxn1_stp1 	 all 1 0.0 1, while for a capture radius of 1.5 this would become react rxn1_stp1 	 all 1 0.0 1.5. When using a capture radius of 1.5, one has to add neighbor 5.0 bin (just before equilibration) in main.lmp/main_customised_fix.lmp and change the simulation box in the data file to -1000 1000 for every dimension.

## MATLAB

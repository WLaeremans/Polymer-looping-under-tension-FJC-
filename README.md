# Polymer-looping-under-tension-FJC-
Code used for publication: Polymer dynamics under tension: mean first passage time for looping

## LAMMPS
There are 2 main files: 1) main.lmp and 2) main_customised_fix.lmp. The former checks every 1000 timesteps whether the polymer has formed a loop and prints the timestep at which a loop is formed to the file bond_info.${seed}.txt. Due to this if statement however, this code might be slow, especially when one wants to check loop formation with a higher accurary than every 1000 timesteps. This issue is fixt in the latter main file, for which we have slightly changed the LAMMPS source code. 


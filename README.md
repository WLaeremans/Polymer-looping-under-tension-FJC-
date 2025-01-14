# Polymer-looping-under-tension-FJC-
Code used for publication in Physical Review E: "Polymer dynamics under tension: mean first passage time for looping", by Wout Laeremans, Anne Floor den Ouden, Jef Hooyberghs and Wouter G. Ellenbroek.

Code and data available on Zenodo: https://doi.org/10.5281/zenodo.14614951

## LAMMPS
For the LAMMPS code, LAMMPS (2 Aug 2023) was used. There are 3 main files: 1) main.lmp and 2) main_customised_fix.lmp and 3) main_crossing.lmp. The former checks every 1000 timesteps whether the polymer has formed a loop and prints the timestep at which a loop is formed to the file bond_info.${seed}.txt. Due to this if statement however, this code might be slow, especially when one wants to check loop formation with a higher accurary than every 1000 timesteps. This issue is fixed in the latter main file, for which we have slightly changed the LAMMPS source code. We changed fix_bond_react.cpp as follows: https://github.com/wouterel/lammps/commit/314e9ca4e797eabb391b2d5b575cb6f5e259d4e0. With this fix, main_customised_fix.lmp prints the exact timestep a loop is formed, without requiring an if statement. Finally, main_crossing.lmp checks every 1000 timesteps whether the x-coordinates have crossed.

The code is now written for a polymer containing 50 beads. When one wants to change this, the data file polymer50.data has to be modified, as well as the main file. In the main file, one should change the id of the right bead to the new number of beads. Furthermore, the code applies a force of 1 (LJ units) on both ends of the polymer, but this can also be changed to a different number. 

Finally, the capture radius is 1 (LJ units) in the code here. When one wants to change this, one has to adapt the third number in the first reaction step. For example, now it says react rxn1_stp1 	 all 1 0.0 1, while for a capture radius of 1.5 this would become react rxn1_stp1 	 all 1 0.0 1.5. When using a capture radius of 1.5, one has to add neighbor 5.0 bin (just before equilibration) in the main file and change the simulation box in the data file to -1000 1000 for every dimension.

The 3 reaction files rxn.map, rxn_pre.data_template and rxn_post.data_template have no further use than simply defining a reaction to efficiently determine when a loop is formed. 

## MATLAB
For the MATLAB code, MATLAB R2023a was used. The MATLAB files for the analysis are catogarised in 3 main subfolders, containing code and associated data. The first is Looping_without_force, which is used for Figure 2 in the article. To run the code, one should open de file FJC_plot1.m and load the folder Data. This folder should be added to the main path. The Data files are named as follows: Nx_fy_rcz.txt, where x is the number of joints (so x+1 beads), y the force and z the capture radius. However, when the capture radius is one (which we see as the default, this is not specified, for example N14_f1.txt. For capture radii 0.5 and 1, main.lmp was used to generate. Hence, they only have an accuracy of +- 500 timesteps. Therefore, in the matlab code we subtract 500 timesteps to average this out. Every data file contains 100 numbers, as 100 simulations were performed in LAMMPS for every set. 

The second subfolder contains the code and data to produce Figure 3 in the article. Again, the Data folder should be added to the main path. The data with N = 19 and rc = 1 is accurate up to +- 500 timesteps. All the other data files were produced with the customised fix. 

The final subfolder contains 3 MATLAB files, FJC_plot3 used to produce Figure 4 in the article, FJC_plot4 used to produce Figure 6 in the article and FJC_plot5 used to produce Figure 7 in the article. Concerning the final file, code is added to perform the nummerical integration using the exact distribution for N = 19 (dotted line in Figure), which takes quite some time to run. However, in the Data folder, we have added as well data_exact_solution_N19.mat, which can be loaded and plotted as plot(f_range,log(tau),'LineWidth',1,'LineStyle', '--', 'Color', 'k');, then the user does not have to run the part solving the integral nummerically.  



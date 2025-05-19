Required Packages
-----------------------

numpy

pandas

scipy

matplotlib

lammps

os

sys

random

subprocess

stringio

latex

scienceplots

csv

Installing Lammps
-----------------------

A LAMMPS executible for Windows can be installed on the following site: https://packages.lammps.org/windows.html

This will also install the Python library.

Input File
-----------------------

The input file is a CSV file, to be formatted identically to the given params.csv file. ```t_low``` and ```t_high``` represent the boundaries for the temperature range, in K. Pressure represents the desired pressure of the simulation, in bar. Pair style represents the type of interaction parameter simulated, which must be LJ unless more interaction parameters are added. Intraspecies represents the inclusion of intraspecies interactions, which will default to none unless H2 is entered for H2 intraspecies interactions. Diffusion type represents the diffuction calculation performed by LAMMPS, which defaults to MSD (Mean Square Diffusion) and can be changed to VACF.
Below is an example of how the input file should be formatted:

```
Species 1 Species 2  Particle Number 1  Particle Number 2  Pressure  Pair Style  Diffusion Type
Sn        H2         200                200                0.1       LJ          MSD
        
t_low   t_high  Intraspecies  num_sim
300     1000    NA            10
```

The following are important values that may need to be changed within the code. The found_tolerance is the tolerance at which the simulation converges to. The default value is 7e-5 and works for most simulations. You may need to increase this value to 7e-4 or higher when testing at greater temperatures. The neighbor_skins values are found from multiple simulations for specific pressures. This value effects the speed of the simulation as well as the dangerous, non-physical builds present in the simulation. All of these values are for given pressures at 200x200 particles. If you are changing the particle count or want to simulate a new pressure at 200x200 particles, you must run the simulation at the new parameters and analyze the LAMMPS MPI task timing breakdown in terminal. In this breakdown, the key sections are Pair and Neigh. A testament to an efficent simulation will have the total percentage of task time for these sections be about the same. They should both take about 10-15% of the total task time with the Modify section taking the majority of the time. If the simulation is still taking long, you can decrease the Neigh task time to about 1% while keeping Pair at about 10-15%. It is also crucial to analyze the Dangerous build just under the LAMMPS MPI task timing breakdown. It is important that no dangerous builds are created during the run. To decrease the Neigh task time, you must increase the neighbor_skin value for that pressure. If you increase it by too much, the Pair task time will skyrocket and you will get dangerous builds. If you decrease the neighbor_skin value for that pressure by too much, the Neigh task time will skyrocket and the simulation will take days to complete. The negihbor_skins values are crucial to running this code effectively, so each value must be custom for a differing volume size, particle count, and, therefore, pressure. The neighbor_skin values provided are a good guide to follow for future neighbor skin values, although they may be quite different for higher particle counts. I much lower neighbor skin length may be needed for more particles to keep the system physical and without dangerous builds. 

Running the Simulation
-----------------------

In console:

```
$ python3 diffusionSnH2.py diffusionparams.csv
```

Output
-----------------------

Outputs a graph of Pressure (Pa) times the Diffusion Coefficient ($m^{-2}s^{-1}$) vs the temperature in Kelvin, for each LAMMPS simulation. Plots a graph for Tin-Hydrogen interactions as found in literature. Below is an example of a generic output, using an input file with 10 temperature steps.

![Updated sample plot](graphs/Error_Graph_Lennard-Jones_1.0bar.png)

Below is an example of a graph developed from 10 tests. It plots the average value from each test at every point along with its standard deviation to two sigma. The following is a guide to get this error graph. At the end of a single test, an output file named "final_data_*pressure*_*fileseed*.txt" is created. After any number of tests of the same parameters, put all of the final_data files into a single directory and run geterrorgraph.py from within that directory. It also includes functionality for a fitted graph where it can be fit the final data up to a desired temperature. The values given_pressure and parameter_type need to be changed within the code for the titles on the graph. 

![Updated sample plot](graphs/Sample_Error_Max_Graph_1.0bar)

Below is an example of a graph developed from multiple error graphs all plotted on one graph from different pressures. This is done similarly to geterrorgraph.py but with getgradientgraph.py.

![Updated sample plot](graphs/All_Pressures.png)

Theory
-----------------------

diffusionSnH2.py generates a given amount of $H_2$ molecules and $Sn$ atoms randomly within a cube of a specific size at a given pressure with periodic boundaries. They are assigned velocities corresponding to the temperature of the system. A Lennard-Jones potential is used to simulate interactions - intraspecies interactions are neglected unless given otherwise while interspecies interactions are considered. A Nose-Hoover thermostat is used to fix the temperature of the ensemble.

LAMMPS first undergoes an minimzation run that checks for overallping particles. The data gathering run is a variable number of steps of 1 femtosecond. LAMMPS outputs the diffusion, temperature, pressure and timestep, every 5,000 steps.

The first method to calculate the diffusion coefficient used is through the velocity autocorrelation function. This can be set by appending the vacf variable to the total list while creating the input file. The diffusion coefficient is calculated through the following steps:

$v_i=v_x(t),v_y(t),v_z(t)$

$\langle v_i(0) \cdot v_i(t) \rangle = \frac{1}{N} \sum ^N _{i=1} v_i(t=t_0) \cdot v_i(t=t_0+\Delta t)$

Where N is the total amount of molecules.

$D = \frac{1}{3}\int^t_0\langle v_i(0) \cdot v_i(t) \rangle dt$


This integral is evaluated numerically using a trapezoidal Reimann sum. After a long time, found to be about 500 ps, it stops increasing and becomes constant. For each simulation, the diffusion coefficient is averaged between t = 3ns and t = 5ns. This value is then multiplied by the average pressure across this period of time.

This process is repeated for each temperature in the temperature range specified. Data obtained through this method is very noisy although a general trend is visibile in the data. This code can be modified to run with more atoms, which may solve this issue.

The mean square displacement is also used to find the diffusion coefficient. This is the default but can be set by appending the msd variable to the total list while creating the input file. The relation for finding $D$ is shown below.

$\langle r^2 \rangle = 6Dt+C$

The slope of the MSD graph plotted against time is used to find the diffusion coefficient; this value is found to converge at the given tolerance. Data is averaged between the final 250,000 timesteps before convergence, and is found the be accurate to previous results. It is recommended to use msd over vacf.


Useful Links
-----------------------
https://docs.lammps.org/compute_msd.html

https://docs.lammps.org/compute_vacf.html

http://utkstair.org/clausius/docs/mse614/pdf/diffusion_intro_v02.pdf

https://c4science.ch/source/lammps/browse/master/examples/DIFFUSE/


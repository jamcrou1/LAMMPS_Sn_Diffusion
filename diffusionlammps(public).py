import numpy as np
import pandas as pd
import scipy.stats as stat
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lammps import lammps
import os
import sys
import random

#new imports
import subprocess
from io import StringIO

#to run lammps simulation
lmp = lammps()

#CONSTANTS
#POTENTIAL PARAMETERS
#SnH2 (Elliot & Weaver) Lennard-Jones Paramters
EPSILON_SNH2 = 0.102378
SIGMA_SNH2 = 2.62975

#H2H2 (Weaver) Lennard_Jones Parameters
EPSILON_H2 = 0.0055148
SIGMA_H2 = 2.7828

#Add additional paramters of your choosing here, Make sure to add a new LJ pair style to acceptable_LJ_pair and edit the pair style in diffusionparams.csv

#PHYSICAL PARAMETERS
MU_SN = 118.71
MU_H2 = 2.01568

#setting initial velocity from temperature
def get_velocity(mu,T):
    #Converts atomic mass(mu) to energy units
    m = 931494697.25613*mu #ev/c^2

    kb = 8.6173303e-5 #boltzmann (eV)

    #Calculates root mean square velocity
    v = np.sqrt(3*kb*T/m)*2.998*10**6 #angstrom/picosecond
    return v


#Writing LAMMPS input file
#Is called to write a new file for each simulation where the temperature varies between simulations
def write_input(filename):
    with open(filename,'w') as f:
        log = 'log output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt\n'
        intro = '''units metal
        dimension 3
        atom_style atomic
        boundary p p p
        timestep 0.001\n'''
        acceptable_LJ_pair = ['LJ']
        if pair_style in acceptable_LJ_pair:
            intro += 'pair_style lj/cut 13\n'
        else:
            raise ValueError('Invalid Pair Style\nAccepted Pair Styles: ' + acceptable_LJ_pair)
        #pair_style hybrid lj/cut 13 morse 13 #This changes the pair_style from only Lennard-Jones to both Lennard-Jones and Morse Potential Parameters

        #Eliminates all particle interactions, make sure to uncomment pairs pair_coeff * * below
        #intro = 'pair_style zero 5 nocoeff\n' 
	    
        #Setup the simulation box as well as two halves of the box where one can spawn H2 in the bottom box and Sn in the top box
        regions = ('region simulation_box block -' + box_size + ' ' + box_size + ' -' + box_size + ' ' + box_size + ' -' + box_size + ' ' + box_size + '\n'
        + 'create_box 2 simulation_box\n'
        + 'region top_box block -' + box_size + ' ' + box_size + ' 0 ' + box_size + ' -' + box_size + ' ' + box_size + '\n'
        + 'region bottom_box block -' + box_size + ' ' + box_size + ' -' + box_size + ' 0 -' + box_size + ' ' + box_size + '\n')

        #Creates atoms with particle number, velocity, and mass parameters
        create = ('create_atoms 1 random '+str(particle_number1)+' '+seed1+' simulation_box \n'
        + 'group ' + name1 + ' id 1:' + str(particle_number1)+' \n' 
        + 'create_atoms 2 random ' + str(particle_number2) + ' ' + seed2 + ' simulation_box \n'
        + 'group ' + name2 + ' id ' + str(particle_number1+1) + ':' + str(particle_number1+particle_number2) + '\n')

        setup = ('mass 1 ' + str(mu1) + '\nmass 2 ' + ' ' + str(mu2) + '\n'
        + 'velocity ' + name1 + ' create ' +str(v1) + ' ' + seed1 + ' dist gaussian\n'
        + 'velocity ' + name2 + ' create ' +str(v2) + ' ' + seed2 + ' dist gaussian\n')
        
        #Defines the intraspecies and interspecies interactions. Decides whether or not to use intraspecies interactions based on the intraspecies parameter
        #Uses Lennard-Jones parameters found from literature
        pairs = '' 
        #Sn intraspecies interactions are currently not present but can be added on request.
        pairs += 'pair_coeff 1 1 0 0 0\n'

        if pair_style == 'LJ':
            pairs += 'pair_coeff 1 2 ' + str(EPSILON_SNH2) + ' ' + str(SIGMA_SNH2) + '\n'
        else:
            raise ValueError('Invalid Pair Style')

        if intraspecies == 'H2' or intraspecies == 'SnH2':
            pairs += 'pair_coeff 2 2 ' + str(EPSILON_H2) + ' ' + str(SIGMA_H2) + '\n'
        else:
            pairs += 'pair_coeff 2 2 0 0 0\n'
        
        #Eliminates all particle interactions, make sure to uncomment intro pair_style zero above as well
        #pairs = 'pair_coeff * *\n'

        #Fixes the temperature of the simulation to be constant throughout a run. Minimizes the energy at the start of the run so that no particles overlap.
        fixes = 'fix nvt all nvt temp ' + str(temperature) + ' ' + str(temperature) + ' $(100.0*dt)\n' #Changed from 100.0*dt to 10.0
        minimize = 'minimize 1.0e-4 1.0e-6 1000 10000\n'
        neighbors = 'neighbor ' + str(neighbor_skins[pair_style][desired_pressure]) + ' bin\nneigh_modify every 10 delay 0 check yes exclude type 2 2 exclude type 1 1 \n'
        reset = 'reset_timestep 0\n'

        #Computes the mean square displacement of all of the Sn particles in every direction. Takes the slope of the MSD to calculate the diffusion coefficient.
        diffusion_sim = ''
        msd = '''compute msd_Sn all msd
        variable diffusion equal c_msd_Sn[4]/6/((step*dt)+1.0e-15)\n'''
        vacf = '''compute vacf_Sn Sn vacf
        fix 5 all vector 1 c_vacf_Sn[4]
        variable diffusion equal 0.333*dt*trap(f_5)\n'''
        if diffusion_type == 'VACF':
            diffusion_sim = vacf
        elif diffusion_type == 'MSD':
            diffusion_sim = msd
        else:
            raise ValueError('The diffusion type must be MSD or VACF')

        #Secondary Method of computing the mean square displacement
        #fix 9 all vector 10 c_msd_Sn[4]
        #equal slope(f_9)/6/10/dt\n'
        
        #Prints the diffusion coefficient value every 5000 timesteps to a new file for residual testing. Thermo decides what variables are printed on screen
        thermo = 'thermo_style custom step temp v_diffusion press\n'

        residuals = ('fix 2 all print ' + str(thermo_timesteps) + ' "${diffusion}" file diffusion_output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt screen no\n')
        #+ 'variable simTemperature equal temp\n' 
        #+ 'fix 3 all print ' + str(thermo_timesteps) + ' "${simTemperature}" file temperature_output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt screen no\n') 
        #Above lines are needed to create a temperature vs time plot for debugging purposes


        #Creates dump files for visualization in VMD. First runs for 1,000 timesteps taking an image for every timestep, then resets. 
        #Afterwards, applies a dump to take an image every 50,000 timesteps for the whole run.

        # dump = 'dump myDump all atom 5 dump_' + str(round(temperature)) + 'K_' + str(file_seed) + '.msd\n'
        # dump_run = 'run 50000\n' 
        # undump = 'undump myDump\n'

        #Uncomment below and add to 'total' after 'residuals' to get a dump file for the whole run with frames every 250000 timesteps
        #dump_all = 'dump myDumpAll all atom 250000 dump_all_' + str(round(temperature)) + 'K_' + str(file_seed) + '.msd\n'

        #Sets the tolerance to where the simulation converges
        tolerance = ('variable toleranceFinal equal ' + str(found_tolerance) + '\n' 
        + 'variable toleranceTimestep equal 0.001\n'
        + 'variable startDiffusion equal 1.0e-15\n'
        + 'label loop\n')

        #Generate_output_from_csv reads from the end of the data backwards until it hits a run command.
        #This can cause errors for reading the output, so be careful adding code after run!
        #Runs for 250,000 timesteps, checks if the fractional difference of the start and the end is more than the tolerance. 
        #If true, jumps back up to run an extra 250,000 timesteps. If false, ends the run.
        run = 'thermo ' + str(thermo_timesteps) + '\nrun ' + str(run_timesteps) + '\n' 
        loop_run = '''variable endDiffusion equal ${diffusion}
        variable difference equal abs(${endDiffusion}-${startDiffusion})/(${startDiffusion})
        variable startDiffusion equal ${endDiffusion}
        print "TOLERANCE: ${toleranceFinal}"
        print "DIFFERENCE: ${difference}"
        if "(${difference} > ${toleranceFinal})" then "jump SELF loop"\n'''

        clear = 'clear'
        #Uncomment below for dump files
        #total  = [log,intro,regions,create,setup,pairs,fixes,diffusion_sim,minimize,reset,neighbors,dump,dump_run,undump,reset,dump_all,thermo,residuals,tolerance,run,loop_run,clear]
        total  = [log,intro,regions,create,setup,pairs,fixes,diffusion_sim,neighbors,minimize,reset,thermo,residuals,tolerance,run,loop_run,clear]
        f.writelines(total)
        f.close()


#Generates the Pandas dataframe from the output file, starts from the end of the output file and reads backwards until it gets the lines between
#"Loop time" and "run ". This ensures that the output dataframe gets the correct variables from the output even with changes to the LAMMPS code
def generate_output_from_csv(logpath):
    with open(logpath,'r') as f:
        lines = f.readlines()
        f.close()

    start_idx = None
    end_idx = None
    end_found = False

    for i, line in enumerate(reversed(lines)):
        if "run " in line and end_found:
            start_idx = len(lines)-i+2
            break
        if "Loop time" in line:
            end_idx = len(lines)-i-1
            end_found = True
    
    relevant_lines = lines[start_idx:end_idx]
    relevant_content = ''.join(relevant_lines)

    df = pd.read_csv(StringIO(relevant_content),sep = '\s+')

    return df


#Generates a custom output with specific data for diagnostic purposes. Is not used for any other data collection or manipulation.
#Easily changeable if one wants to write specific data here.
def generate_specified_ouptut_from_csv(logpath, filename):
    with open(logpath,'r') as f:
        lines = f.readlines()
        f.close()

    timestep_lines = []
    neighbor_lines = []

    previous_line = None
    for i, line in enumerate(lines):
        if "Loop time" in line:
            timestep_lines.append((previous_line.split())[0])
        if "Total # of neighbors" in line:
            neighbor_lines.append((line.split()[5]))
        previous_line = line

    open_as = ''
    if sim_number == 1:
        open_as = 'w'
    else:
        open_as = 'a'

    with open(filename,open_as) as f:
        if open_as == 'w':
            f.writelines(filename)

        simulation = '\nSimulation #' + str(sim_number) + ' Seed #' + seed1 + "\n"
        timestep_line = ''
        for i in range(len(timestep_lines)):
            timestep_line += str(timestep_lines[i]) + ': ' + str(neighbor_lines[i]) + ' neighbors\n'

        total  = [simulation,timestep_line]
        f.writelines(total)
        f.close()


#Generates a custom error output with fitslope and pressure data. This file is used to create an average graph with error for multiple runs in geterrorgraph.py
def generate_error_output_data(logpath, filename):
    open_as = ''
    if sim_number == 1:
        open_as = 'w'
    else:
        open_as = 'a'

    df = generate_output_from_csv(logpath)
    temp = df['Temp'].astype(float)
    pavg = np.average(temp[1:])*1.380649*10**(-23)*(particle_number1+particle_number2)/((side_length)**3*10**(-30))

    D_values = df['v_diffusion']*10**(-8)
    with open(filename,open_as) as f:
        line = str(pavg) + ','
        for value in D_values:
            line += str(value) + ','
        line = line[:-1]
        line += '\n'

        f.writelines(line)
        f.close()


#Reads the fitslope data printed onto a new file in the LAMMPS code. 
#Generates plots of the log of the residuals over time and fitslope over time. Used to find and prove the best tolerance value.
def plot_residuals(logpath):
    diffusion_list = []
    with open(logpath,'r') as f:
        for line in f:
            diffusion = line.strip()
            diffusion_list.append(diffusion)
        f.close()

    diffusion_list.pop(0)
    diffusion_list = [float(x)*10**(-8) for x in diffusion_list]
    diffusion_final_list = []
    [diffusion_final_list.append(x) for x in diffusion_list if x not in diffusion_final_list and x < 1000]
    timestep_fit_list = list(range(0, len(diffusion_final_list)))
    timestep_fit_list = [x * 5000 for x in timestep_fit_list]

    diffusion_old = 0
    difference_list = []
    for diffusion_new in diffusion_list:
        if float(diffusion_old) == 0:
            difference_list.append(0)
        elif float(diffusion_old) == float(diffusion_new):
            continue
        else:
            difference = (abs(float(diffusion_old) - float(diffusion_new)))/float(diffusion_old)
            if (difference > 10**-10 and difference < 10):
                difference_list.append(difference)
        diffusion_old = diffusion_new

    timestep_list = list(range(0, len(difference_list)))
    timestep_list = [x * 5000 for x in timestep_list]

    plt.figure(dpi=400)
    plt.plot(timestep_list,difference_list,'.',label='Residuals')
    plt.yscale('log')
    plt.xlabel('Timesteps (ps)')
    plt.ylabel('Fractional Difference')
    plt.title('Log plot of the residuals at '+ str(round(temperature,2)) + ' K')
    plt.legend()
    png = 'Residuals_' + str(particle_number1) + 'x' + str(particle_number2) + 'Particles_' + str(round(side_length)) + 'Å_' + intraspecies + 'Intraspecies_' + str(round(temperature)) + 'K_'+ str(num_sim) + 'Steps_' + str(file_seed) + '.png'
    plt.savefig(png)

    plt.figure(dpi=400)
    plt.plot(timestep_fit_list,diffusion_final_list,'.',label='Diffusion Coefficient')
    plt.xlabel('Timesteps (ps)')
    plt.ylabel('Diffusion Coefficient ($\\frac{m^2}{s}$)')
    plt.title('Diffusion Coefficient at ' + str(round(temperature,2)) + ' K')
    plt.legend()
    png2 = 'DiffusionCoefficient_' + str(particle_number1) + 'x' + str(particle_number2) + 'Particles_' + str(round(side_length)) + 'Å_' + intraspecies + 'Intraspecies_' + str(round(temperature)) + 'K_'+ str(num_sim) + 'Steps_' + str(file_seed) + '.png'
    plt.savefig(png2, bbox_inches='tight')
    plt.close()


def plot_temperature(logpath):
    temp_dataframe = pd.read_csv(logpath, header=None, names=['Temperatures'])
    temp_list = temp_dataframe['Temperatures'].tolist()
    temp_list = temp_list[2:]
    temp_list = [float(x) for x in temp_list]

    timestep_list = list(range(0, len(temp_list)))
    timestep_list = [x * 5000 for x in timestep_list]

    plt.figure(dpi=400)
    plt.plot(timestep_list,temp_list,'.',label='Temperature')
    plt.xlabel('Timesteps (ps)')
    plt.ylabel('Temperature (K)')
    plt.title('Temperature Variance at '+ str(round(temperature,2)) + ' K',y=1.05)
    plt.suptitle('Average Temperature: ' + str(np.round(np.average(temp_list),2)) + ' K, Standard Deviation: ' + str(np.round(np.std(temp_list),3)),fontsize=10,y=0.92)
    plt.legend()
    png = 'Temperature_' + str(particle_number1) + 'x' + str(particle_number2) + 'Particles_' + str(round(side_length)) + 'Å_' + intraspecies + 'Intraspecies_' + str(round(temperature)) + 'K_'+ str(num_sim) + 'Steps_' + str(file_seed) + '.png'
    plt.savefig(png)
    plt.close()


#Calculates the diffusion coefficient*pressure
def get_pD_msd(logpath):
    df = generate_output_from_csv(logpath)
    print("V_DIFFUSION AVG: " + str(np.average(df['v_diffusion'])))
    D = df['v_diffusion']*10**(-8) #l^2/t

    #Takes the average of all of the values in the last run. This would be 50 values over a 250,000 timestep run where convergence occured
    Davg = np.average(D)
    
    #Takes the temperature values throughout the run. They should always average to the input temperature since we hold temperature constant in LAMMPS.
    #As the run goes on, the temperatures vary much more and may get to +/-10 of the input temperature.
    #temp = df['Temp'].astype(float)
    #print("TEMPERATURE: " + str(np.average(temp[1:])))
    print("TEMPERATURE: " + str(temperature))


    #Average pressure equation is average temp * bultzmanns constant * number of particles / (volume)^3 * unit converter value
    #Pressure value given in LAMMPS can be negative due to the attractive vs repulsive forces between particles, so not used in this case
    pavg = temperature*1.380649*10**(-23)*(particle_number1+particle_number2)/((side_length)**3*10**(-30))
    
    pD = pavg*Davg
    return pD
    

#Defines all of the parameter values
#Runs the number of simulations by writing and running the lammps file for each temperature
if __name__ == '__main__':
    params = pd.read_csv(sys.argv[1])

    #General Particle Parameters    
    name1 = params['Species 1'][0]
    name2 = params['Species 2'][0]
    particle_number1 = int(params['Particle Number 1'][0])
    particle_number2 = int(params['Particle Number 2'][0])
    intraspecies = str(params['Intraspecies'][0])
    pair_style = str(params['Pair Style'][0])
    diffusion_type = str(params['Diffusion Type'][0])
    #Change mu if testing different species
    mu1 = MU_SN
    mu2 = MU_H2

    #General Simulation Parameters
    num_sim = params['num_sim'][0]
    side_length = 0
    box_size = 0
    median_temp = (float(params['t_high'][0]) + float(params['t_low'][0])) /2
    temperatures = np.linspace(float(params['t_low'][0]),float(params['t_high'][0]),num=int(params['num_sim'][0]))
    desired_pressure = float(params['Pressure'][0]) #In Barr
    run_timesteps = 250000
    thermo_timesteps = int(run_timesteps/50)
    found_tolerance = 7*10**-5

    #Neighbor Parameters
    #Need to add new neighbor skins here for every new pair style added, this includes adding specific values at deisred pressures. 
    neighbor_skins = ({'LJ': {5.0: 2, 1.0: 5, 0.5: 8, 0.1: 10, 0.05: 20, 0.025: 20, 0.01: 40, 0.005: 60, 0.001: 200}})
    if desired_pressure not in neighbor_skins[pair_style].keys():
        raise ValueError('Need to add a neighbor skin for desired pressure into neighbor_skins')

    pD = []
    sim_number = 1
    file_seed = str(random.randint(100000,999999))
    error_datafile = 'final_data_' + str(desired_pressure) + 'barr_' + file_seed + '.txt'
    specific_output_datafile = str(particle_number1) + 'x' + str(particle_number2) + 'Particles_' + str(desired_pressure) + 'Barr_' + intraspecies + 'Intraspecies_VariableTimesteps_' + str(num_sim) + 'Steps_' + str(file_seed) + '.txt'

    for temperature in temperatures:
        print(f'Simulation Number: {sim_number}')
        #We determine the side length and therefore pressure of the simulation here. The two methods are described more in the write up.
        #Our default method uses a changing volume and constant pressure
        side_length = ((temperature * 1.380649*10**(-23) * (particle_number1+particle_number2)) / (desired_pressure * 10**5 * 10**-30))**(1/3)
        #side_length = ((median_temp * 1.380649*10**(-23) * (particle_number1+particle_number2)) / (desired_pressure * 10**5 * 10**-30))**(1/3)

        box_size = str(side_length/2)

        seed1  = str(random.randint(100000,999999))
        seed2 = str(random.randint(100000,999999))
        input_datafile = 'input_' + str(round(temperature)) + 'K_' + str(file_seed) + '.lmp'
        v1 = get_velocity(mu1,temperature)
        v2 = get_velocity(mu2,temperature)
        write_input(input_datafile)

        #lmp.file executes the LAMMPS input script
        lmp.file(input_datafile)
        
        output_datafile = 'output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt'
        
        generate_error_output_data(output_datafile, error_datafile)
        pressure_diffusion = get_pD_msd(output_datafile)
        pD.append(pressure_diffusion)
        residual_datafile = 'diffusion_output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt'
        plot_residuals(residual_datafile)

        #Uncomment for temperature data 
        #temperature_datafile = 'temperature_output_' + str(round(temperature)) + 'K_' + str(file_seed) + '.txt'
        #plot_temperature(temperature_datafile)
        #Uncomment for specific output for easier debugging purposes
        #generate_specified_ouptut_from_csv(output_datafile, specific_output_datafile)

        sim_number += 1


    #Plots our found computed data 
    plt.figure(dpi=200)
    plt.plot(temperatures,pD,'.-',label='LAMMPS')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Diffusion Coefficient * Pressure ($\\frac{Pa{\\cdot}m^2}{s}$)')
    plt.title('P*D vs T for ' + str(particle_number1) + ' Sn and ' + str(particle_number2) + ' ' + str(name2) + ' at ' + str(desired_pressure) + ' Barr')
    plt.legend()
    png = str(particle_number1) + 'x' + str(particle_number2) + 'Particles_' + str(desired_pressure) + 'Barr_' + intraspecies + 'Intraspecies_' + str(num_sim) + 'Steps_' + str(file_seed) + '.png'
    plt.savefig(png)
    plt.close()


import matplotlib.pyplot as plt
from matplotlib import cm, colors
import scienceplots
import pandas as pd
import numpy as np
import glob
import sys
import csv

if __name__ == '__main__':
    pressures = [1.0,0.5,0.1,0.05,0.01,0.005,0.001] #Change for pressures tested
    plt.figure(dpi=200)
    plt.style.use(['science'])

    cmap = plt.get_cmap('viridis', len(pressures))
    p = 0
    all_PD = []
    all_Temps = []

    Temps = [300,377.77777778,455.55555556,533.33333333,611.11111111,688.88888889,766.66666667,844.44444444,922.22222222,1000]
    for pressure in pressures:
        file_pattern = 'final_data_' + str(pressure) + '*.txt'

        file_list = glob.glob(file_pattern)

        lines = [[] for _ in range(10)]

        for file_name in file_list:
            with open(file_name, 'r') as file:
                file_lines = file.readlines()
                for i in range(10):
                    pavg = 0
                    j = 0

                    for value in file_lines[i].split(','):
                        if j == 0:
                            pavg = float(value)
                        else: 
                            lines[i].append(float(value)*pavg)
                        j += 1

        pD = []
        Error = []
        for line in lines:
            pD.append(np.average(line))
            Error.append((np.std(line))*2)

        #CONTOUR PLOT
        # all_PD.append(pD)
        # all_Temps.append(Temps)

        #NORMAL PLOT
        plt.plot(Temps,pD,'.-',label=(str(pressure) + ' bar'),color=cmap(p))
        #plt.errorbar(Temps,pD,yerr = Error,fmt = '.',capsize=5,color=cmap(p)) #color=cmap(len(pressures)-p)

        p += 1 #Loop End


    parameter_type = 'LJ'
    

    # #CONTOUR PLOT
    # plt.contourf(Temps,pressures,all_PD)
    # plt.xlabel('Temperature (K)', fontsize=14)
    # plt.ylabel('Pressure (bar)', fontsize=14)
    # plt.yscale('log')
    # plt.title('P$\cdot$D modeled by Temperature and Pressure', fontsize=14)
    
    # Contour Plat Extras
    # norm = colors.Normalize(0,64)
    # cbar = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=plt.get_cmap('viridis')),ax=plt.gca())
    # cbar.set_label('P$\cdot$D ($\\frac{Pa{\\cdot}m^2}{s}$)\n'+parameter_type+' Parameters')
    # cbar.ax.get_yaxis().set_ticks([0,8,16,24,32,40,48,56,64])
    # png = 'Contour_'+parameter_type+str(len(pressures))+'.png'
    

    # #NORMAL PLOT
    plt.ylabel('P$\cdot$D ($\\frac{Pa{\\cdot}m^2}{s}$)', fontsize=14)
    plt.xlabel('Temperature (K)', fontsize=14)
    plt.title('P$\cdot$D vs T', fontsize=14)
    plt.legend(fontsize=13)
    png = 'Gradient_'+parameter_type+str(len(pressures))+'.png'

    plt.savefig(png)

import matplotlib.pyplot as plt
import latex
import scienceplots
plt.rcParams.update(plt.rcParamsDefault)
import pandas as pd
import numpy as np
import glob
import sys
import csv
from scipy.optimize import curve_fit


def fit_func(x,a,b):
    return a*x**b

def temperature_range(max_temp):
    temps_max = [300,377.77777778,455.55555556,533.33333333,611.11111111,688.88888889,766.66666667,844.44444444,922.22222222,1000] 
    if max_temp == '1700':
        temps_max.extend([1000,1077.77777778,1155.55555556,1233.33333333,1311.11111111,1388.88888889,1466.66666667,1544.44444444,1622.22222222,1700])
    elif max_temp == '2400':
        temps_max.extend([1000,1077.77777778,1155.55555556,1233.33333333,1311.11111111,1388.88888889,1466.66666667,1544.44444444,1622.22222222,1700,1777.77777778,1855.55555556,1933.33333333,2011.11111111,2088.88888889,2166.66666667,2244.44444444,2322.22222222,2400])
    elif max_temp == '3100':
        temps_max.extend([1000,1077.77777778,1155.55555556,1233.33333333,1311.11111111,1388.88888889,1466.66666667,1544.44444444,1622.22222222,1700,1777.77777778,1855.55555556,1933.33333333,2011.11111111,2088.88888889,2166.66666667,2244.44444444,2322.22222222,2400,2477.77777778,2555.55555556,2633.33333333,2711.11111111,2788.88888889,2866.66666667,2944.44444444,3022.22222222,3100])
    else:
        print('Default temperature range of 300-1000K used')
    return temps_max

if __name__ == '__main__':
    file_pattern = 'final_data_*.txt'
    plt.style.use(['science']) #AIP
    #Example styles: ieee, pgf, nature, notebook
    file_list = glob.glob(file_pattern)

    Temps = [0,300,377.77777778,455.55555556,533.33333333,611.11111111,688.88888889,766.66666667,844.44444444,922.22222222,1000]

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

    temps = [300,377.77777778,455.55555556,533.33333333,611.11111111,688.88888889,766.66666667,844.44444444,922.22222222,1000]
    given_pressure = input('What is the pressure? ')
    parameter_type = input('File name? ')

    fit_check = input('Do you want a fitted graph? ')

    cmap = plt.get_cmap('winter', 2)
    plt.figure(dpi=325)
    png = ''
    if fit_check.lower() == 'yes' or fit_check.lower() == 'y':
        range_check = input('What max temperature do you want to fit to? [1700], [2400], or [3100]')
        fit_temps = temperature_range(range_check)
        params = curve_fit(fit_func,temps,pD)
        [a,b] = params[0]
        fit_pD = fit_func(fit_temps,a,b)
        plt.plot(temps,pD,'.-',label=('Sn H2 Data'),color=cmap(0))
        plt.plot(fit_temps,fit_pD,'.-',label=('Fit Function: b = ' + str(round(b,3)) + '\na = ' + str(round(a,8))),color=cmap(1))
        png = 'Error_Max_Fitted_Graph_'+parameter_type+'_'+given_pressure+'bar.png'
    else:
        plt.plot(temps,pD,'.-',label='Sn H2 Data')
        png = 'Error_Max_Graph_'+parameter_type+'_'+given_pressure +'bar.png'
    
    plt.errorbar(temps,pD,yerr = Error,fmt = '.',capsize=5,color='rebeccapurple')

    #Uncomment for max relative error of point on graph
    # relative_error = np.array(Error) / np.array(pD)
    # max_error = max(relative_error)
    # max_error_index = (np.where(relative_error == max_error))[0]
    # print(relative_error)
    #plt.plot(np.array(temps)[max_error_index], np.array(pD)[max_error_index]+np.array(Error)[max_error_index],'.',zorder=100, label="Max Relative Standard Deviation: \ncv = " + str(round(max_error,8)),color='turquoise')

    plt.xlabel('Temperature (K)')
    plt.ylabel('P$\cdot$D ($\\frac{Pa{\\cdot}m^2}{s}$)')
    plt.title('P$\cdot$D vs T at ' + given_pressure + ' bar')
    plt.legend()
    plt.savefig(png)

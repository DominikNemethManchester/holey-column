# -*- coding: utf-8 -*-
"""
v2.0
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import csv



#GEOMETRY
L = 12.0 #extra material at ends
D = 8.5 # distance between holes
RADIUS = 4.0
E = 6*(123.4e3 + 19.23e3)
density = 1160

#E = 2e11
#density = 7500
speed_of_sound = np.sqrt(E/density)


#FITTING
def func(x, m, c):
    return m*x + c


#DATA

BASE_DIRECTORY = "C:/Users/f51139dn/UNIFORM_colum_compression_v6/"
SIZES = [10]
MESH_SIZE = 0.2
VISCOSITIES = np.array([0.1, 1, 2.5, 5])
VELOCITIES = np.array([0.01, 0.05, 0.1, 0.15, 0.2]) * 5164


for velocity in VELOCITIES:
    size=10
    num_holes = size
    for viscosity in VISCOSITIES:
        # subdir = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(100*velocity)) + "uniform"
        subdir = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(velocity)) + "_ELASTIC"
        # subdir = "Job" + "-h" + str(int(num_holes)) + "-v" + str(int(100*viscosity)) + "-hyper" + "-i" + str(int(100*velocity))
        input_directory = BASE_DIRECTORY + subdir
        os.chdir(input_directory)
        print(os.getcwd())
        file_name = "compression_wave_results_" + subdir + ".csv"
        compression_data = np.genfromtxt(file_name, delimiter=',', skip_header=True)
        if np.size(compression_data) < 4:
            continue
        compression_data = compression_data[compression_data[:, 1].argsort()]
        u2 = compression_data[:,0]
        t = compression_data[:,1]

        distances = [(L+ D/2 + n * D) for n in range(0, num_holes)]

        # TRANSIENCE
        # Threshold for the difference in consecutive y-values
        t_diff_threshold = 0.1* D/speed_of_sound  # Adjust this threshold based on your data

        # Find the index where the difference is larger than the threshold
        change_point = np.argmax(np.abs(np.diff(t)) > t_diff_threshold)

        # Skip transient data
        t_filtered = t[change_point + 1:]
        u2_filtered = u2[change_point + 1:]

        # Skip transient data
        t_filtered = t[change_point:]
        u2_filtered = u2[change_point:]

        xdata = t_filtered
        ydata = u2_filtered
        
        filtered_results = np.column_stack((t_filtered, u2_filtered))
        EXTENSION = '-RESULTS'
        results_path = BASE_DIRECTORY + subdir + EXTENSION
        if os.path.exists(results_path) == False:
            os.mkdir(results_path)
        os.chdir(results_path)
        np.savetxt("compression_wave_results_filtered_" + subdir + ".csv", filtered_results, delimiter=',')

        popt, pcov = curve_fit(func, xdata, ydata)
        m, c = popt

        #PLOTTING
        fig = plt.figure(figsize=(8,5))
        ax = fig.add_subplot(111)
        ax.tick_params(which='both', direction='in', labelsize=10)
        ax.plot(t, func(t, *popt), linewidth=2, ls='-', color='k')
        ax.plot(t, u2, ls='none', marker='d', markerfacecolor='orange', markeredgecolor='black')

        ax.set_xlabel('Time (s)', fontsize=14)
        ax.set_ylabel('Distance (m)', fontsize=14)
        ax.legend(fontsize=11, loc='upper left', bbox_to_anchor=(-0.01,1))
        ax.set_xlim(0)
        ax.set_ylim(0)

        output_path = 'compression_wave_speed_' + subdir
        output_pdf_path = output_path + '.pdf'
        plt.savefig(output_pdf_path, dpi = 200)
        plt.show()

        #WRITE FILE FOR SPEEDS
        output_csv_path = output_path + '.csv'
        with open(output_csv_path, 'w', newline='') as output_csv:
            writer = csv.writer(output_csv)

        # Write header
            writer.writerow(
                ['v_c (m/s)', 'v_c / c', 'gradient', 'offset'])

            writer.writerow([m, m/speed_of_sound, m, c])
            print(m / speed_of_sound)

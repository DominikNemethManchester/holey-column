# -*- coding: utf-8 -*-
"""
v2.0
"""

import numpy as np
import os
import matplotlib.pyplot as plt

#PLOTTING
#plt.style.use(['science', 'ieee'])
plt.rcParams.update({'figure.dpi': '200'})

BASE_DIRECTORY = "C:/Users/f51139dn/UNIFORM_colum_compression_v6/"
SIZES = [10]
MESH_SIZE = 0.2
VISCOSITIES = [5]
# INDENTER_SPEEDS = np.array([0.1, 1, 1.5, 2.0, 3, 4, 5])
VELOCITIES = np.array([0.01, 0.05, 0.1, 0.15, 0.2]) * 5164

# E = 6*(123.4e3 + 19.23e3)
# density = 1160
E = 2e11
density = 7500
speed_of_sound = np.sqrt(E/density)

# FIX TWO OF THESE
FIXED_SIZE = True
FIXED_VISCOSITY = True
FIXED_INDENTER_SPEED = False

# SET FIXED VALUES 
if FIXED_SIZE:
    size = 10  # number of holes
if FIXED_VISCOSITY:
    v = 5 # linear viscosity parameter
if FIXED_INDENTER_SPEED:
    i = 10/6 # Indenter Speed

# SET UP ARRAYS
velocities_array = np.array([])
compression_array = np.array([])
buckling_array = np.array([])
time_lag_array = np.array([])

# VARYING INDENTER SPEED 
if not FIXED_INDENTER_SPEED:
    for i in VELOCITIES:
        #JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper-i' + str(int(100 * i))
        #JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper10'
        JOB = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(i)) + "_ELASTIC"
        # JOB = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(i*100)) + "uniform"
        EXTENSION = '-RESULTS'
        try:
            os.chdir(BASE_DIRECTORY + JOB + EXTENSION)
        except FileNotFoundError:
            continue

        compression_filename = 'compression_wave_speed_' + JOB + ".csv"
        #buckling_filename = 'buckling_wave_speed_' + JOB + '.csv'
        #time_lag_filename = 'time_lag_' + JOB + '.csv'

        compression_data = np.genfromtxt(compression_filename, delimiter=',', skip_header=True)
        #buckling_data = np.genfromtxt(buckling_filename, delimiter=',', skip_header=True)
        #time_lag_data = np.genfromtxt(time_lag_filename, delimiter=',', skip_header=True)

        # WAVE SPEEDS normalised by speed of sound (not using v because I have been using that for un-normalised speeds)
        # n means normalised, v means velocity
        v_c = compression_data[0]
        # nv_b = buckling_data[1]
        #t_lag = time_lag_data
        ydata = (v_c - 0)/speed_of_sound
        velocities_array = np.append(velocities_array, i)
        compression_array = np.append(compression_array, ydata)
        #buckling_array = np.append(buckling_array, nv_b)
        #time_lag_array = np.append(time_lag_array, t_lag)
    
    # CHANGE DIRECTORY BEFORE SAVING PLOT
    EXTENSION = '-FINAL_PLOTS'
    RUN = '-s' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper'
    results_path = BASE_DIRECTORY +  "uniform_m_" + str(int(MESH_SIZE*100)) + '-FINAL_PLOTS1e-4'
    if not os.path.exists(results_path):
        os.mkdir(results_path)
    os.chdir(results_path)
    # COMPRESSION WAVE
    fig1 = plt.figure(figsize=(8, 5))
    ax1 = fig1.add_subplot(111)
    ax1.tick_params(which='both', direction='in', labelsize=10)
    ax1.set_xlabel('$v_0 / c$ ', fontsize=16)  #(ms$^{-1}$)
    ax1.set_ylabel('$v_c /c$', fontsize=16)
    #plt.ylim((1.0, 1.2))
    x = np.array(velocities_array) / speed_of_sound
    ax1.plot(x, compression_array, ls='None', marker='v', markerfacecolor='red', markeredgecolor='black')
    plt.savefig('indenter_speed_compression_default' + RUN + '.pdf')
    # BUCKLING WAVE
    # fig2 = plt.figure(figsize=(8, 5))
    # ax2 = fig2.add_subplot(111)
    # ax2.tick_params(which='both', direction='in', labelsize=10)
    # ax2.set_xlabel('$v_0$ (ms$^{-1}$)', fontsize=16)
    # ax2.set_ylabel('$v_b /c$', fontsize=16)
    # ax2.plot(INDENTER_SPEEDS, buckling_array, ls='None', marker='d', markerfacecolor='green', markeredgecolor='black')
    # plt.savefig('indenter_speed_buckling' + RUN + '.pdf')
    # # TIME LAG
    # fig3 = plt.figure(figsize=(8, 5))
    # ax3 = fig3.add_subplot(111)
    # ax3.tick_params(which='both', direction='in', labelsize=10)
    # ax3.set_xlabel('$v_0$ (ms$^{-1}$)', fontsize=16)
    # ax3.set_ylabel('Time lag (s)', fontsize=16)
    # ax3.plot(INDENTER_SPEEDS, time_lag_array, ls='None', marker='s', markerfacecolor='blue', markeredgecolor='black')

    # plt.savefig('indenter_speed_time_lag' + RUN + '.pdf')
    plt.show()

# # VARYING VISCOSITY
# if not FIXED_VISCOSITY:
#     for v in VISCOSITIES:
#         #JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper-i' + str(int(100 * i))
#         JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper10'
#         EXTENSION = '-RESULTS'
#         os.chdir(BASE_DIRECTORY + JOB + EXTENSION)

#         compression_filename = 'compression_wave_speed_' + JOB + ".csv"
#         buckling_filename = 'buckling_wave_speed_' + JOB + '.csv'
#         time_lag_filename = 'time_lag_' + JOB + '.csv'

#         compression_data = np.genfromtxt(compression_filename, delimiter=',', skip_header=True)
#         buckling_data = np.genfromtxt(buckling_filename, delimiter=',', skip_header=True)
#         time_lag_data = np.genfromtxt(time_lag_filename, delimiter=',', skip_header=True)

#         # WAVE SPEEDS normalised by speed of sound (not using v because I have been using that for un-normalised speeds)
#         # n means normalised, v means velocity
#         nv_c = compression_data[1]
#         nv_b = buckling_data[1]
#         t_lag = time_lag_data

#         compression_array = np.append(compression_array, nv_c)
#         buckling_array = np.append(buckling_array, nv_b)
#         time_lag_array = np.append(time_lag_array, t_lag)


#     # CHANGE DIRECTORY BEFORE SAVING PLOT
#     EXTENSION = '-FINAL_PLOTS'
#     RUN = '-s' + str(int(size))  + '-hyper-i' + str(int(100 * i))
#     results_path = BASE_DIRECTORY +  'run10per6'
#     if not os.path.exists(results_path):
#         os.mkdir(results_path)
#     os.chdir(results_path)
#     # COMPRESSION WAVE
#     fig1 = plt.figure(figsize=(8, 5))
#     ax1 = fig1.add_subplot(111)
#     ax1.tick_params(which='both', direction='in', labelsize=10)
#     ax1.set_xlabel('$b_1$ (Pa)', fontsize=16)
#     ax1.set_ylabel('$v_c /c$', fontsize=16)
#     ax1.plot(VISCOSITIES, compression_array, ls='None', marker='v', markerfacecolor='red', markeredgecolor='black')
#     plt.savefig('viscosity_compression' + RUN + '.pdf')
#     # BUCKLING WAVE
#     fig2 = plt.figure(figsize=(8, 5))
#     ax2 = fig2.add_subplot(111)
#     ax2.tick_params(which='both', direction='in', labelsize=10)
#     ax2.set_xlabel('$b_1$ (Pa)', fontsize=16)
#     ax2.set_ylabel('$v_b /c$', fontsize=16)
#     ax2.plot(VISCOSITIES, buckling_array, ls='None', marker='d', markerfacecolor='green', markeredgecolor='black')
#     plt.savefig('viscosity_buckling' + RUN + '.pdf')
#     # TIME LAG
#     fig3 = plt.figure(figsize=(8, 5))
#     ax3 = fig3.add_subplot(111)
#     ax3.tick_params(which='both', direction='in', labelsize=10)
#     ax3.set_xlabel('$b_1$ (Pa)', fontsize=16)
#     ax3.set_ylabel('Time lag (s)', fontsize=16)
#     ax3.plot(VISCOSITIES, time_lag_array, ls='None', marker='s', markerfacecolor='blue', markeredgecolor='black')
#     plt.savefig('viscosity_time_lag' + RUN + '.pdf')
#     plt.show()

# # VARYING SIZE
# if not FIXED_SIZE:
#     for size in SIZES:
#         #JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper-i' + str(int(100 * i))
#         JOB = 'Job-h' + str(int(size)) + '-v' + str(int(100 * v)) + '-hyper10'
#         EXTENSION = '-RESULTS'
#         os.chdir(BASE_DIRECTORY + JOB + EXTENSION)

#         compression_filename = 'compression_wave_speed_' + JOB + ".csv"
#         buckling_filename = 'buckling_wave_speed_' + JOB + '.csv'
#         time_lag_filename = 'time_lag_' + JOB + '.csv'

#         compression_data = np.genfromtxt(compression_filename, delimiter=',', skip_header=True)
#         buckling_data = np.genfromtxt(buckling_filename, delimiter=',', skip_header=True)
#         time_lag_data = np.genfromtxt(time_lag_filename, delimiter=',', skip_header=True)

#         # WAVE SPEEDS normalised by speed of sound (not using v because I have been using that for un-normalised speeds)
#         # n means normalised, v means velocity
#         nv_c = compression_data[1]
#         nv_b = buckling_data[1]
#         t_lag = time_lag_data
#         compression_array = np.append(compression_array, nv_c)
#         buckling_array = np.append(buckling_array, nv_b)
#         time_lag_array = np.append(time_lag_array, t_lag)
        
#     # CHANGE DIRECTORY BEFORE SAVING PLOT
#     EXTENSION = '-FINAL_PLOTS'
#     RUN = '-v' + str(int(100 * v)) + '-hyper-i' + str(int(100 * i))
#     results_path = BASE_DIRECTORY +  'run10per6'
#     if not os.path.exists(results_path):
#         os.mkdir(results_path)
#     os.chdir(results_path)

#     # COMPRESSION WAVE
#     fig1 = plt.figure(figsize=(8, 5))
#     ax1 = fig1.add_subplot(111)
#     ax1.tick_params(which='both', direction='in', labelsize=10)
#     ax1.set_xlabel('Number of holes', fontsize=16)
#     ax1.set_ylabel('$v_c /c$', fontsize=16)
#     ax1.plot(SIZES, compression_array, ls='None', marker='v', markerfacecolor='red', markeredgecolor='black')
#     plt.savefig('size_compression' + RUN + '.pdf')
#     # BUCKLING WAVE
#     fig2 = plt.figure(figsize=(8, 5))
#     ax2 = fig2.add_subplot(111)
#     ax2.tick_params(which='both', direction='in', labelsize=10)
#     ax2.set_xlabel('Number of holes', fontsize=16)
#     ax2.set_ylabel('$v_b /c$', fontsize=16)
#     ax2.plot(SIZES, buckling_array, ls='None', marker='d', markerfacecolor='green', markeredgecolor='black')
#     plt.savefig('size_buckling' + RUN + '.pdf')
#     # TIME LAG
#     fig3 = plt.figure(figsize=(8, 5))
#     ax3 = fig3.add_subplot(111)
#     ax3.tick_params(which='both', direction='in', labelsize=10)
#     ax3.set_xlabel('Number of holes', fontsize=16)
#     ax3.set_ylabel('Time lag (s)', fontsize=16)
#     ax3.plot(SIZES, time_lag_array, ls='None', marker='s', markerfacecolor='blue', markeredgecolor='black')

    

#     plt.savefig('size_time_lag' + RUN + '.pdf')
#     plt.show()

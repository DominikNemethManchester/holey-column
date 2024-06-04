import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from scipy.optimize import curve_fit
import os
import csv

#plt.style.use(['science', 'ieee'])
#plt.rcParams.update({'figure.dpi': '200'})

L = 12.0e-3 #extra material at ends
D = 8.5e-3 # distance between holes
RADIUS = 4.0e-3

E = 6*(123.4e3 + 19.23e3)
density = 1160
THEORETICAL_C = np.sqrt(E/density)

BASE_DIRECTORY = "C:/Users/f51139dn/VISCOSITY_SWEEP_HOLEY_HYPERELASTIC/"
OUTPUT_SUBDIR = "FINAL_PLOTS"
N_HOLES = 13
MESH_SIZE = 0.2e-3
COLUMN_LENGTH = 2*L + 2*RADIUS + (N_HOLES - 1)*D


VISCOSITIES = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
VELOCITIES = np.array([0.025, 0.05, 0.1, 0.15, 0.2]) * THEORETICAL_C

EXPERIMENTAL_COMPRESSION_SIGMA = 0.0025


ELLIPTICITY_PLOTS = False


# Data at constant indenter speed
data_at_i = []
speed_for_output = 0.1 * THEORETICAL_C


def func(x, m, c):
    return m*x + c


def filter_buckling_data(times, distances, v_compression):
    end_index = np.size(distances)
    for i in range(np.size(distances)-1):
        if times[i+1] < times[i]:
            compression_arrival_time = (2*COLUMN_LENGTH - distances[i]) / v_compression
            if times[i] > compression_arrival_time:
                end_index = i
                break
    
    filtered_times = times[:end_index]
    filtered_distances = distances[:end_index]
    filtered_indices = np.arange(end_index)
    
    return filtered_times, filtered_distances, filtered_indices


def get_peaks_ratios(e, t):
    first_peak_index = np.argmax(e)
    first_peak = [e[first_peak_index], t[first_peak_index]]
    
    for i in range(first_peak_index, np.size(e)):
        if e[i] < 1.01: break
    
    e_copy = np.delete(e, np.arange(0, i))
    second_peak_index = np.argmax(e_copy) + i
    second_peak = [e[second_peak_index], t[second_peak_index]]
    
    e_copy = e[first_peak_index:second_peak_index]
    trough_index = np.argmin(e_copy) + first_peak_index
    trough = [e[trough_index], t[trough_index]]
    
    return first_peak, second_peak, trough


ratios = np.zeros((np.size(VISCOSITIES), 1))


for v_index, v in enumerate(VISCOSITIES):
    compression_velocities = np.zeros(np.shape(VELOCITIES))
    self_contact_velocities = np.zeros(np.shape(VELOCITIES))
    buckling_wave_velocities = np.zeros(np.shape(VELOCITIES))
    
    for vel_index, velocity in enumerate(VELOCITIES):
        # buckling for v_0 / c <= 0.025 cannot be treated as a buckling wave
        # if velocity < 0.03 * THEORETICAL_C:
        #     continue
        
        RUN = "h" + str(int(N_HOLES)) + "_m" + str(int(100*1000*MESH_SIZE)) + "_v" + str(int(v*100)) +"_i"+ str(int(velocity*100)) + "_hyper"
        try:
            os.chdir(BASE_DIRECTORY + RUN)
        except FileNotFoundError:
            continue
        
        FILE_START = 'results_'+ "m" + str(int(100*1000*MESH_SIZE)) + "_h" + str(int(N_HOLES))  + "_v" + str(int(v*100)) +"_speed"+ str(int(velocity*100))
        data = np.genfromtxt('nodes_coords_circlenos_2d_'+'m' + str(int(100*1000*MESH_SIZE)) + '_h' + str(int(N_HOLES)) +'.csv', delimiter=',', skip_header=1)
    
        index_buckling_onset = np.array([])
        times_buckling_onset = np.array([])
        times_compression_onset = np.array([])
        index_self_contact_onset = np.array([])
        times_self_contact = np.array([])
        distances = np.array([])
        distances_c = np.array([])
        distances_sc = np.array([])
        
        for i in range(N_HOLES):
            lines = np.where(data[:, 3] == i)[0]
            index_offset = lines[0]
            x_values = data[lines, 1]
            delta_x = np.max(x_values) - np.min(x_values)
            x_average = np.average(np.unique(x_values))
            x_1 = x_average - 0.25 * delta_x
            x_2 = x_average + 0.25 * delta_x
            y_values = data[lines, 2]
            delta_y = np.max(y_values) - np.min(y_values)
            y_average = np.average(np.unique(y_values))
            y_1 = y_average - 0.25 * delta_y
            y_2 = y_average + 0.25 * delta_y
    
            # QUANDRANT ALIGNMENTS: y increases upwards, x increases to the right
            quadrant_left_indices = np.where(data[lines, 1] < x_1)[0] + index_offset  # left quadrant
            quadrant_right_indices = np.where(data[lines, 1] > x_2)[0] + index_offset  # right quadrant
            quadrant_top_indices = np.where(data[lines, 2] > y_2)[0] + index_offset  # top quadrant
            quadrant_bottom_indices = np.where(data[lines, 2] < y_1)[0] + index_offset  # bottom quadrant
    
            # knowing what lines corresponds to each quadrant we can index into those and find the node numbers
            quadrant_left = data[quadrant_left_indices, 0]
            quadrant_right = data[quadrant_right_indices, 0]
            quadrant_top = data[quadrant_top_indices, 0]
            quadrant_bottom = data[quadrant_bottom_indices, 0]
            
            quadrant_top_data = []
            for node in quadrant_top:
                node_data = np.genfromtxt(FILE_START + f'_node{node:.0f}_u2.csv', delimiter=',',
                                          skip_header=1)
                node_data_y = node_data[:, 1]
                quadrant_top_data.append(list(node_data_y))
            quadrant_top_data = np.mean(quadrant_top_data, axis=0)
    
            quadrant_bottom_data = []
            for node in quadrant_bottom:
                node_data = np.genfromtxt(FILE_START + f'_node{node:.0f}_u2.csv', delimiter=',',
                                          skip_header=1)
                node_data_y = node_data[:, 1]
                quadrant_bottom_data.append(list(node_data_y))
            quadrant_bottom_data = np.mean(quadrant_bottom_data, axis=0)
    
            quadrant_left_data = []
            for node in quadrant_left:
                node_data = np.genfromtxt(FILE_START + f'_node{node:.0f}_u1.csv', delimiter=',',
                                          skip_header=1)
                node_data_x = node_data[:, 1]
                quadrant_left_data.append(list(node_data_x))
            quadrant_left_data = np.mean(quadrant_left_data, axis=0)
    
            quadrant_right_data = []
            for node in quadrant_right:
                node_data = np.genfromtxt(FILE_START + f'_node{node:.0f}_u1.csv', delimiter=',',
                                          skip_header=1)
                node_data_x = node_data[:, 1]
                quadrant_right_data.append(list(node_data_x))
            quadrant_right_data = np.mean(quadrant_right_data, axis=0)
            
            times = node_data[:, 0]
            b = 0.5 * (2 * RADIUS + quadrant_right_data - quadrant_left_data)
            a = 0.5 * (2 * RADIUS + quadrant_top_data - quadrant_bottom_data)
            e = b / a
            
            # SELF CONTACT WAVE
            lower_threshold = MESH_SIZE / (2 * RADIUS)
            upper_threshold = 1 / lower_threshold
            for k, element in enumerate(e):
                if (element < lower_threshold) or (element > upper_threshold):
                    index_self_contact = k
                    time_self_contact = times[index_self_contact]
                    distance_sc = (L + D / 2 + i * D)
                    
                    index_self_contact_onset = np.append(index_self_contact_onset, index_self_contact)
                    times_self_contact = np.append(times_self_contact, time_self_contact)
                    distances_sc = np.append(distances_sc, distance_sc)
                    break
            
            # COMPRESSION WAVE
            normalised_e = e - 1
            for k, element in enumerate(normalised_e):
                if np.abs(element) > 5 * EXPERIMENTAL_COMPRESSION_SIGMA:
                    index_compression = k
                    time_compression = times[index_compression]
                    distance = (L + D / 2 + i * D)
                    
                    distances_c = np.append(distances_c, distance)
                    times_compression_onset = np.append(times_compression_onset, time_compression)    
                    break
                
            
            if velocity < 0.03 * THEORETICAL_C:
                continue
            
            
            # BUCKLING WAVE
            normalised_e = e - 0.99
            if i%2 == 1:
                continue
            for k, element in enumerate(normalised_e):
                if element < 0:
                    index_buckling = k
                    time_buckling = times[index_buckling]
                    distance = (L + RADIUS + i * D)
                    
                    distances = np.append(distances, distance)
                    index_buckling_onset = np.append(index_buckling_onset, index_buckling)
                    times_buckling_onset = np.append(times_buckling_onset, time_buckling)    
                    break
                
                
            # TEST PLOTS
            if ELLIPTICITY_PLOTS and i==6:
                plt.plot(times, e, color='k')
                # plt.plot(times, a, color='blue')
                # plt.plot(times, b, color='red')
                plt.scatter(time_compression, e[index_compression], s=20, color='k')
                plt.scatter(time_buckling, e[index_buckling], s=20, color='k')
                plt.scatter(time_self_contact, e[index_self_contact], s=20, color='k')
                plt.title(i)
                plt.xlim(0, 0.010)
                plt.ylim(0.8, 1.2)
                
                print(f"t_c: {time_compression:.4f}, e={e[index_compression]}")
                print(f"t_b: {time_buckling:.4f}, e={e[index_buckling]}")
                print(f"t_sc: {time_self_contact:.4f}, e={e[index_self_contact]}")
                
                # np.savetxt(BASE_DIRECTORY + "e_plot_hole10_" + RUN + ".csv", np.column_stack((times, e)), delimiter=',')
                
                first_peak, second_peak, trough = get_peaks_ratios(e, times)
                plt.scatter(first_peak[1], first_peak[0], color='r', s=20)
                plt.scatter(second_peak[1], second_peak[0], color='r', s=20)
                plt.scatter(trough[1], trough[0], color='r', s=20)
                ratios[v_index, 0] = (second_peak[0]- trough[0]) / (first_peak[0] - trough[0])
                
                plt.show()
        
                        
        # SAVE RESULTS IN NEW FOLDER
        EXTENSION = '-RESULTS'
        results_path = BASE_DIRECTORY + RUN + EXTENSION
        if os.path.exists(results_path) == False:
            os.mkdir(results_path)
        os.chdir(results_path)
        
        # PREPARE TEXT OUTPUT
        print(f"--- VISCOSITY = {v}, VELOCITY = {velocity/THEORETICAL_C}c ---")
        speed_strings = ['N/A', 'N/A']
        
        
        ####
        # START COMPRESSION WAVE
        ####
        
        if np.size(distances_c) > 2:
            popt_c, pcov_c = curve_fit(func, times_compression_onset, distances_c)
            m_c, c_c = popt_c
            # speed_strings[0] = f'v_sc = {m_sc / THEORETICAL_C:.3f}c'
            compression_velocities[vel_index] = m_c
        else:
            print("No compression wave")
            continue
        
        ####
        # END COMPRESSION WAVE
        ####
        
        
        ####
        # START SELF CONTACT
        ####
        PLOT_SELF_CONTACT = False
        
        if np.size(distances_sc) > 2:
            PLOT_SELF_CONTACT = True
            popt_sc, pcov_sc = curve_fit(func, times_self_contact, distances_sc)
            m_sc, c_sc = popt_sc
            speed_strings[0] = f'v_sc = {m_sc / THEORETICAL_C:.3f}c'
            self_contact_velocities[vel_index] = m_sc
        else:
            print("No self contact")
        
        ####
        # END SELF CONTACT
        ####
        
        
        ####
        # BUCKLING WAVE
        ####
        xdata = times_buckling_onset
        ydata = distances
        PLOT_BUCKLING = False
        
        if np.size(distances) >= 3:
            # results = np.column_stack((times_buckling_onset, distances))
            #np.savetxt('buckling_wave_results_' + JOB + '.csv', results, delimiter=',')
        
            # FILTER OUT OUTLIERS
            # xdata = xdata.reshape(-1, 1)
            # model = make_pipeline(PolynomialFeatures(1), LinearRegression())
            # model.fit(xdata, ydata)
    
            # y_pred = model.predict(xdata)
            # residuals = ydata - y_pred
            # threshold = 3 * np.std(residuals)
        
            # filtered_indices = np.where(np.abs(residuals) <= threshold)[0]
            # xdata = xdata.flatten()
            # filtered_x = xdata[filtered_indices]
            # filtered_y = ydata[filtered_indices]
            
            filtered_x, filtered_y, filtered_indices = filter_buckling_data(times_buckling_onset, distances, m_c)
            
            # FIT LINEAR TREND
            popt, pcov = curve_fit(func, filtered_x, filtered_y)
            m, c = popt
            buckling_wave_velocities[vel_index] = m
            
            # BUCKLING SPEED
            output_csv_path = 'buckling_wave_speed_' + RUN + 'test.csv'
            with open(output_csv_path, 'w', newline='') as output_csv:
                writer = csv.writer(output_csv)
                writer.writerow(['v_b (ms-1)', 'v_b / c', 'gradient', 'offset'])
                writer.writerow([m, m / THEORETICAL_C, m, c])
            
            speed_strings[1] = f"v_b = {m / THEORETICAL_C:.3f}c"
            
            
            PLOT_BUCKLING = True
        else:
            print("No buckling wave\n")

        
        print(f"{speed_strings[1]}  |   {speed_strings[0]} \n")
    
        # PLOT DATA AND FIT
        fig = plt.figure(figsize=(8, 5))
        ax = fig.add_subplot(111)
        
        ax.scatter(times_compression_onset, distances_c, marker='o', color='k',
                    label="Compression wave")
        ax.plot(times_compression_onset, func(times_compression_onset, m_c, c_c),
                linewidth=2, dashes=[2, 2], color='k')
        
        if(PLOT_SELF_CONTACT):
            ax.scatter(times_self_contact, distances_sc, marker='v', color='b',
                        label="Self contact wave")
            ax.plot(times_self_contact, func(times_self_contact, m_sc, c_sc),
                    linewidth=2, ls='-', color='k')
        
        if(PLOT_BUCKLING):
            ax.plot(xdata, func(xdata, *popt), linewidth=2, ls='--', color='k')
            ax.scatter(np.delete(xdata, filtered_indices), np.delete(ydata, filtered_indices),
                        marker='v', color='white', edgecolor='red')
            ax.scatter(filtered_x, filtered_y, marker='v', color='red',
                        label="Buckling wave")
            
            if velocity == speed_for_output:
                data_at_i.append(np.array([v, m / THEORETICAL_C, m_sc / THEORETICAL_C]))
        
        ax.set_title(f"Waves in holey column for viscosity={v}, speed={velocity/THEORETICAL_C}c")
        ax.set_xlabel('Time (s)', fontsize=16)
        ax.set_ylabel('Distance (m)', fontsize=16)
        ax.tick_params(which='both', direction='in', labelsize=10)
    
        # plt.xlim(0, 0.01)
        plt.legend(fontsize=12)
        plt.savefig('buckling_wave_results_' + RUN + 'test.pdf')
        plt.show()


    # FINAL SUMMARY PLOTS
    if not os.path.isdir(BASE_DIRECTORY + OUTPUT_SUBDIR):
        os.mkdir(BASE_DIRECTORY + OUTPUT_SUBDIR)
    os.chdir(BASE_DIRECTORY + OUTPUT_SUBDIR)
    
    good_buckling_indices = np.where(buckling_wave_velocities != 0)[0]
    good_sc_indices = np.where(self_contact_velocities != 0)[0]
    good_c_indices = np.where(compression_velocities != 0)[0]
    buckling_wave_velocities /= THEORETICAL_C
    self_contact_velocities /= THEORETICAL_C
    compression_velocities /= THEORETICAL_C
    
    fig = plt.figure(figsize=(8, 5))
    ax = fig.add_subplot(111)

    ax.scatter(VELOCITIES[good_buckling_indices] / THEORETICAL_C,
               buckling_wave_velocities[good_buckling_indices],
               marker='s', color='red', label='Buckling wave')
    ax.scatter(VELOCITIES[good_sc_indices] / THEORETICAL_C,
               self_contact_velocities[good_sc_indices],
               marker='v', color="blue", label='Self contact wave')
    # ax.scatter(VELOCITIES[good_c_indices] / THEORETICAL_C,
    #            compression_velocities[good_c_indices],
    #            marker='o', color="k", label='Compression wave')
    
    ax.set_title(f"Wave speeds in holey column, viscosity={v}", fontsize=16)
    ax.set_xlabel(r'$v_0/c$', fontsize=16)
    ax.set_ylabel(r'$v/c$', fontsize=16)
    ax.tick_params(which='both', direction='in', labelsize=10)
    
    plt.ylim(0, 0.6)
    plt.xlim(0)
    plt.legend(fontsize=12, loc='lower right')
    # plt.savefig(f'buckling_wave_results_v{v}.pdf')
    plt.show()
    
    
    if v == 50:
        output_data = np.column_stack((VELOCITIES[good_buckling_indices] / THEORETICAL_C,
                                       buckling_wave_velocities[good_buckling_indices],
                                       self_contact_velocities[good_sc_indices]))
        file_header = "v_0/c, v_b/c, v_sc/c"
        np.savetxt(BASE_DIRECTORY + "uniform_m_20-FINAL_PLOTS/v_vs_v0_viscosity50.csv",
                   output_data, header=file_header, delimiter=',')
    
    
    # if not os.path.isdir(BASE_DIRECTORY + OUTPUT_SUBDIR + "/wave speed data/"):
    #     os.mkdir(BASE_DIRECTORY + OUTPUT_SUBDIR + "/wave speed data/")
    # os.chdir(BASE_DIRECTORY + OUTPUT_SUBDIR + "/wave speed data/")
    
    # buckling_results = np.column_stack((VELOCITIES / THEORETICAL_C,
    #                                     compression_velocities,
    #                                     buckling_wave_velocities,
    #                                     self_contact_velocities))
    # np.savetxt(f"v{v}.csv", buckling_results, header="v_i/c, v_c/c, v_b/c, v_s/c", delimiter=",")


if VISCOSITIES[-1] < 30:
    ratios_output = np.column_stack((VISCOSITIES, ratios))
    np.savetxt(BASE_DIRECTORY + "peak_ratios.csv", ratios_output, delimiter=',')


data_at_i = np.row_stack(tuple(data_at_i))
np.savetxt(BASE_DIRECTORY + f"uniform_m_20-FINAL_PLOTS/speeds_at_i{speed_for_output/THEORETICAL_C:.2f}.csv",
           data_at_i, delimiter=',', header="viscosity, v_b/c, v_sc/c")

# -*- coding: utf-8 -*-
"""
v2.0:
    -Compression wave arrival now determined from gradient of
     displacement in time
"""

import numpy as np
from numpy.polynomial import Polynomial as P
import matplotlib.pyplot as plt
import os


E = 6*(123.4e3 + 19.23e3)
DENSITY = 1160
THEORETICAL_C = np.sqrt(E/DENSITY)

DIAMETER = 8.0e-3
HOLE_SEPARATION = 8.5e-3
ENDS_LENGTH = 12.0e-3
MESH_SIZE = 0.2e-3

SIZES = [13]
VISCOSITIES = [50.0]
VELOCITIES = np.array([0.025, 0.05, 0.075, 0.1, 0.125, 0.15]) * THEORETICAL_C


# CONFIGURATIONS
POLY_DEGREE = 35
DISPLACEMENT_CUTOFF = 0.5
REAL_ROOT_THRESHOLD = 1e-5
PADDING_RATIO = 2
SPEED_FRACTION = 0.01

BASE_DIRECTORY = "C:/Users/f51139dn/VISCOSITY_SWEEP_HOLEY_HYPERELASTIC/"

TEST_PLOTS = True

E = 6*(123.4e3 + 19.23e3)
DENSITY = 1160
THEORETICAL_C = np.sqrt(E/DENSITY)


def find_pulse_arrival_old(times, displacements, threshold_gradient, unpadded_cut_time):
    """
    Returns y coordinate, time and numpy.Polynomial object
    corresponding to point when the gradient of y vs t
    is equa to the threshold gradient
    """
    poly_series = P.fit(times, displacements, POLY_DEGREE, domain=[])
    poly_deriv = poly_series.deriv()

    gradient_poly = P([threshold_gradient])
    poly_deriv = poly_deriv - gradient_poly

    poly_roots = poly_deriv.roots()
    real_roots = poly_roots.real[np.abs(poly_roots.imag) < REAL_ROOT_THRESHOLD]
    real_roots = real_roots[np.where(real_roots > unpadded_cut_time)[0]]
    real_roots = np.sort(real_roots)
    if np.size(real_roots) == 0:
        print("No real roots in domain: no roots to polynomial")
        return None

    cutoff_time = real_roots[0]
    cutoff_y = poly_series(cutoff_time)
    if cutoff_time > np.max(times):
        print("No real roots in domain: time exceeded")
        return None

    return cutoff_time, cutoff_y, poly_series


def find_inflection_point(times, displacements, threshold_gradient, unpadded_cut_time):
    """
    Returns y coordinate, time and numpy.Polynomial object
    corresponding to point when the gradient of y vs t
    is equa to the threshold gradient
    """
    poly_series = P.fit(times, displacements, POLY_DEGREE, domain=[])
    poly_deriv = poly_series.deriv()
    poly_second_deriv = poly_deriv.deriv()

    # gradient_poly = P([threshold_gradient])
    # poly_deriv = poly_deriv - gradient_poly

    poly_roots = poly_second_deriv.roots()
    real_roots = poly_roots.real[np.abs(poly_roots.imag) < REAL_ROOT_THRESHOLD]
    real_roots = real_roots[np.where(real_roots > unpadded_cut_time)[0]]
    real_roots = np.sort(real_roots)
    if np.size(real_roots) == 0:
        print("No real roots in domain: no roots to polynomial")
        return None
    
    while(np.size(real_roots) > 0 and poly_deriv(real_roots[0]) < threshold_gradient):
        real_roots = np.delete(real_roots, 0)
    if np.size(real_roots) == 0:
        return None

    cutoff_time = real_roots[0]
    cutoff_y = poly_series(cutoff_time)
    if cutoff_time > np.max(times):
        print("No real roots in domain: time exceeded")
        return None
    
    cutoff_gradient = poly_deriv(cutoff_time)

    return cutoff_time, cutoff_y, poly_series, cutoff_gradient


def cut_data(time, displacement):
    cut_index = int(np.min(np.where(displacement > 1e-16)[0]) / PADDING_RATIO)
    unpadded_time = time[cut_index * PADDING_RATIO]
    time = time[cut_index:]
    displacement = displacement[cut_index:]
    
    # DIAMETER = 8.0e-3
    # HOLE_SEPARATION = 8.5e-3
    # ENDS_LENGTH = 12.0e-3
    # LIGAMENT_WIDTH = HOLE_SEPARATION - DIAMETER
    # COLUMN_WIDTH = 2 * DIAMETER + 2 * LIGAMENT_WIDTH
    # COLUMN_LENGTH = (2 * ENDS_LENGTH + 21 * DIAMETER
    #                  + (21 - 1) * LIGAMENT_WIDTH)
    
    # allowed_indices = np.where(time < COLUMN_LENGTH / THEORETICAL_C)
    return time, displacement, unpadded_time



def analyse_data(subdir, num_holes, nodes_coords_circlenos_file, velocity):
    LIGAMENT_WIDTH = HOLE_SEPARATION - DIAMETER
    COLUMN_LENGTH = (2 * ENDS_LENGTH + num_holes * DIAMETER
                     + (num_holes - 1) * LIGAMENT_WIDTH)
    COLUMN_HALF_LENGTH = COLUMN_LENGTH / 2
        
    directory = BASE_DIRECTORY + subdir
    node_file_data = np.genfromtxt(nodes_coords_circlenos_file, delimiter=',',
                                   skip_header=1)
    node_labels = np.array([])
    node_circlenos = np.array([])
    node_y_coordinates = np.array([])
    for data_line in node_file_data:
        node_circleno = data_line[3]
        if node_circleno in node_circlenos:
            continue
        node_labels = np.append(node_labels, data_line[0])
        node_circlenos = np.append(node_circlenos, data_line[3])
        node_y_coordinates = np.append(node_y_coordinates, data_line[2])
    node_y_coordinates += COLUMN_HALF_LENGTH
    
    # nodes_per_circle = np.size(np.where(node_circlenos==0)[0])
    # nodes = node_labels[::nodes_per_circle]
    nodes = node_labels
    
    signal_arrival_times = np.array([])
    new_signal_arrival_times = np.array([])
    signal_arrival_coordinates = np.array([])

    filenames = []
    for filename in os.listdir(directory):
        if 'results_' in filename and '_node' in filename and 'u2' in filename:
            filenames.append(directory + '/' + filename)

    for filename in filenames:
        start_cut = filename.rfind('e') + 1
        node = int(filename[start_cut:-7])
        if node not in nodes:
            continue

        node_index = np.where(node_labels == node)[0]
        circle_no = node_circlenos[node_index]
        node_y = node_y_coordinates[node_index]
        if circle_no <= 1:
            continue
        
        # if circle_no != 11:
        #     continue

        data = np.genfromtxt(filename, delimiter=',', dtype=float)
        times = data[:, 0]
        displacements = data[:, 1]
        
        reflected_wave_arrival_time = (2 * COLUMN_LENGTH - node_y) / THEORETICAL_C
        clean_indices = np.where(times < reflected_wave_arrival_time * 1.5)[0]
        times = times[clean_indices]
        displacements = displacements[clean_indices]
        
        
        times, displacements, unpadded_cut_time = cut_data(times, displacements)
        if np.size(times) < 10:
            print(f"Too few data points for i={velocity/THEORETICAL_C}c, circle={circle_no}")
            continue

        threshold_gradient = velocity * SPEED_FRACTION
        algorithm_outputs = find_inflection_point(times, displacements, threshold_gradient, unpadded_cut_time)
        if not algorithm_outputs:
            print(f"Failed to find pulse at i={velocity/THEORETICAL_C}c, circle={circle_no}")
            continue
        cutoff_time, cutoff_y, poly_series, threshold_gradient = algorithm_outputs
        
        if cutoff_time == None:
            continue
        
        print(circle_no)
        
        zero_crossing_time = cutoff_time - cutoff_y /threshold_gradient
        
        # PLOT TO CHECK CUTOFFS (FOR UNEXPECTED OUTPUTS)
        if TEST_PLOTS:
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111)
            ax.scatter(cutoff_time, cutoff_y, s=10, color='r')
            
            t = np.linspace(np.min(times), np.max(times), 200)
            data = np.genfromtxt(filename, delimiter=',', dtype=float)
            times = data[:, 0]
            displacements = data[:, 1]
            y = poly_series(t)
            ax.plot(times, displacements, label='data', color='k', dashes=[4, 4])
            ax.plot(t, y, label='interpolation', color='b')
            
            t = np.linspace(zero_crossing_time, cutoff_time)
            y = (t - zero_crossing_time) * threshold_gradient
            ax.plot(t, y, color='k', linewidth=1)
                  
            ax.scatter(zero_crossing_time, 0, s=40, color='r')
            ax.tick_params(which='both', direction='in', labelsize=10)
            ax.set_xlabel(r'time (s)', fontsize=13)
            ax.set_ylabel(r'$y$ coordinate (m)', fontsize=13)
            ax.set_title(f'{circle_no}', fontsize=13)
            plt.legend()
            plt.show()

        signal_arrival_coordinates = np.append(signal_arrival_coordinates, node_y)
        signal_arrival_times = np.append(signal_arrival_times, cutoff_time)
        new_signal_arrival_times = np.append(new_signal_arrival_times, zero_crossing_time)


    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    ax.scatter(signal_arrival_times, signal_arrival_coordinates, s=10, color='r')
    ax.set_xlabel('t')
    ax.set_ylabel('U2')
    # plt.xlim(0)
    # plt.ylim(0)
    plt.show()


    output = np.column_stack((signal_arrival_coordinates, signal_arrival_times, new_signal_arrival_times))
    output_filename = "compression_wave_results_" + subdir + ".csv"
    np.savetxt(output_filename, output, delimiter=',',
                header="u2, t, t_new")


def main():
    for num_holes in SIZES:
        nodes_coords_circlenos_file = BASE_DIRECTORY + "nodes_coords_circlenos_2d_m" + str(int(MESH_SIZE*100*1000)) + "_h" + str(int(num_holes)) + ".csv"
        for viscosity in VISCOSITIES:
            for velocity in VELOCITIES:
                subdir = "h" + str(int(num_holes)) + "_m" + str(int(MESH_SIZE*1e5)) + "_v" + str(int(viscosity*100)) +"_i"+ str(int(velocity*100)) + "_hyper"
                # subdir = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(100*velocity)) + "uniform"
                # subdir = "Job" + "-h" + str(int(num_holes)) + "-v" + str(int(100*viscosity)) + "-hyper" + "-i" + str(int(100*velocity))
                output_directory = BASE_DIRECTORY + subdir
                # nodes_coords_circlenos_file = "nodes_coords_circlenos_2d_m" + str(int(MESH_SIZE*100)) + "_h" + str(int(num_holes)) + ".csv"
                os.chdir(output_directory)
                print(os.getcwd())
                analyse_data(subdir, num_holes, nodes_coords_circlenos_file, velocity)


if __name__ == '__main__':
    main()

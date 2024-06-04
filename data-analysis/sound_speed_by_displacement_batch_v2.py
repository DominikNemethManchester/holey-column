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


DIAMETER = 8.0
HOLE_SEPARATION = 8.5
ENDS_LENGTH = 12.0
MESH_SIZE = 0.2

SIZES = [10]

VELOCITIES = np.array([0.01, 0.05, 0.1, 0.15, 0.2]) * 5164
VISCOSITIES = [5]



# CONFIGURATIONS
POLY_DEGREE = 35
DISPLACEMENT_CUTOFF = 0.5
REAL_ROOT_THRESHOLD = 1e-5
SPEED_FRACTION = 0.5

BASE_DIRECTORY = "C:/Users/f51139dn/UNIFORM_colum_compression_v6/"

TEST_PLOTS = False

E = 6*(123.4e3 + 19.23e3)
DENSITY = 1160
THEORETICAL_C = np.sqrt(E/DENSITY)



def cut_data(time, displacement):
    # cut_indices = np.where(displacement > DISPLACEMENT_CUTOFF)[0]
    # if np.size(cut_indices) > 0:
    #     cutoff_index = np.min(cut_indices)
    #     time = time[:cutoff_index]
    #     displacement = displacement[:cutoff_index]
    cut_indices = np.where(displacement < 1e-16)[0]
    if np.size(cut_indices) > 0:
        time = np.delete(time, cut_indices)
        displacement = np.delete(displacement, cut_indices)
    return time, displacement


def find_pulse_arrival(times, displacements, threshold_gradient):
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
    real_roots = real_roots[np.where(real_roots > np.min(times))[0]]
    real_roots = np.sort(real_roots)
    if np.size(real_roots) == 0:
        return None
    
    cutoff_time = real_roots[0]
    cutoff_y = poly_series(cutoff_time)
    if cutoff_time > np.max(times):
        return None
    
    return cutoff_y, cutoff_time, poly_series


def analyse_data(subdir, num_holes, nodes_coords_circlenos_file, velocity):
    LIGAMENT_WIDTH = HOLE_SEPARATION - DIAMETER
    COLUMN_LENGTH = (2 * ENDS_LENGTH + num_holes * DIAMETER
                     + (num_holes - 1) * LIGAMENT_WIDTH)
    COLUMN_HALF_LENGTH = COLUMN_LENGTH / 2
        
    directory = BASE_DIRECTORY + subdir
    node_file_data = np.genfromtxt(nodes_coords_circlenos_file, delimiter=',',
                                   skip_header=1)
    node_labels = node_file_data[:, 0]
    node_circlenos = node_file_data[:, 3]
    node_y_coordinates = node_file_data[:, 2]
    node_y_coordinates += COLUMN_HALF_LENGTH
    
    nodes_per_circle = np.size(np.where(node_circlenos==0)[0])
    # nodes = node_labels[::nodes_per_circle]
    nodes = node_labels
    
    signal_arrival_times = np.array([])
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

        data = np.genfromtxt(filename, delimiter=',', dtype=float)
        times = data[:, 0]
        displacements = data[:, 1]
        
        reflected_wave_arrival_time = (2 * COLUMN_LENGTH - node_y) / THEORETICAL_C
        clean_indices = np.where(times < reflected_wave_arrival_time)[0]
        times = times[clean_indices]
        displacements = displacements[clean_indices]
        
        times, displacements = cut_data(times, displacements)
        if np.size(times) < 10:
            continue

        threshold_gradient = velocity * SPEED_FRACTION
        algorithm_outputs = find_pulse_arrival(times, displacements, threshold_gradient)
        if not algorithm_outputs:
            continue
        cutoff_y, cutoff_time, poly_series = algorithm_outputs
        
        print(circle_no)

        # PLOT TO CHECK CUTOFFS (FOR UNEXPECTED OUTPUTS)
        if TEST_PLOTS:
            fig = plt.figure(figsize=(7, 5))
            ax = fig.add_subplot(111)
            ax.scatter(cutoff_time, cutoff_y, s=10, color='r')
            
            t = np.linspace(np.min(times), np.max(times), 200)
            data = np.genfromtxt(filename, delimiter=',', dtype=float)
            times = data[:, 0]
            displacements = data[:, 1]
            
            ax.plot(times, displacements)
            y = poly_series(t)
            ax.plot(t, y)
            ax.set_xlabel('t')
            ax.set_ylabel('U2')
            # if circle_no == 1:
            #     plt.savefig("v5_detection.pdf", dpi=200)
            plt.show()

        signal_arrival_coordinates = np.append(signal_arrival_coordinates, node_y)
        signal_arrival_times = np.append(signal_arrival_times, cutoff_time)

    fig = plt.figure(figsize=(7, 5))
    ax = fig.add_subplot(111)
    ax.scatter(signal_arrival_times, signal_arrival_coordinates, s=10, color='r')
    ax.set_xlabel('t')
    ax.set_ylabel('U2')
    plt.xlim(0)
    plt.ylim(0)
    plt.show()

    output = np.column_stack((signal_arrival_coordinates, signal_arrival_times))
    output_filename = "compression_wave_results_" + subdir + ".csv"
    np.savetxt(output_filename, output, delimiter=',',
                header="u2, t")


def main():
    for num_holes in SIZES:
        for viscosity in VISCOSITIES:
            for velocity in VELOCITIES:
                subdir = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(velocity)) + "_ELASTIC"
                # subdir = "m" + str(int(MESH_SIZE*100)) + "i" + str(int(100*velocity)) + "uniform"
                # subdir = "Job" + "-h" + str(int(num_holes)) + "-v" + str(int(100*viscosity)) + "-hyper" + "-i" + str(int(100*velocity))
                output_directory = BASE_DIRECTORY + subdir
                nodes_coords_circlenos_file = "UNIFORM_nodes_coords_circlenos_2d_m" + str(int(MESH_SIZE*100)) + "_h" + str(int(num_holes)) + ".csv"
                # nodes_coords_circlenos_file = "nodes_coords_circlenos_2d_m" + str(int(MESH_SIZE*100)) + "_h" + str(int(num_holes)) + ".csv"
                os.chdir(output_directory)
                print(os.getcwd())
                analyse_data(subdir, num_holes, nodes_coords_circlenos_file, velocity)


if __name__ == '__main__':
    main()

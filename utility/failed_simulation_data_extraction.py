# -*- coding: utf-8 -*-
"""
Example script for extracting data from a failed simulation.
"""
import numpy as np
import os
from abaqus import *
from abaqusConstants import *

import mesh
from mesh import ElemType
import regionToolset

MESH_SIZE = 0.2

# GEOMETRY OF PART
NUM_HOLES = 10

SAVED_VARIABLES = ("S", "E", "U")

VELOCITIES = np.array([0.1]) * 5164
viscosity = 5.0

base_directory = "C:/Users/f51139dn/elastic_poisson0.3/"

for velocity in VELOCITIES:
    job_name = "Job-h"+str(int(NUM_HOLES)) + "-i" + str(int(velocity)) + "-uniform"
    node_labels_path = (base_directory + "UNIFORM_nodes_coords_circlenos_2d_m"
                        + str(int(MESH_SIZE*100)) + "_h" + str(int(NUM_HOLES)) + ".csv")
    path = base_directory + "m" + str(int(100*MESH_SIZE)) +"i"+ str(int(velocity)) + "_ELASTIC"
    os.chdir(path)
    
    odb = session.openOdb(job_name + ".odb")
    node_labels = np.genfromtxt(node_labels_path, dtype=int, delimiter=",")
    for node_line in node_labels:
        node = node_line[0]
        
        ### KEY STEP ###
        region_key = 'Node PART-INSTANCE.' + str(node)
        data_u2 = np.array(odb.steps['Step-dynamic'].historyRegions[region_key]\
                            .historyOutputs['U2'].data)
        ### END KEY STEP ###
        
        np.savetxt(base_directory + 'results_m' + str(int(MESH_SIZE*100))
                    + "_h"+str(int(NUM_HOLES)) + "_v" + str(int(100*viscosity))
                    + "_speed"+ str(int(velocity))
                    + '_node' + str(node) + '_u2.csv',
                    data_u2, header='time,U2', delimiter=',')

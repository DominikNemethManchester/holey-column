# -*- coding: utf-8 -*-
"""
[See Johnson (2017) Fig. 7]

Draws N-hole flapped column, with the option to specify the various parameters
defined in Jonhson (2017).

Also does the usual setting up of the materials, sections, assembly, etc.

@author: ricra
"""

import numpy as np
import os
import shutil
from abaqus import *
from abaqusConstants import *
from viewerModules import *
from driverUtils import executeOnCaeStartup

import mesh
from mesh import ElemType
import regionToolset

E = 6*(123.4e3 + 19.23e3)
DENSITY = 1160
THEORETICAL_C = np.sqrt(E/DENSITY)

n_holes = [13]

MESH_SIZE = 0.2e-3
SECTION_THICKNESS = 1.0e-3

NUM_CPUS = 10

VELOCITIES = np.array([0.025, 0.05, 0.075, 0.1, 0.125, 0.15]) * THEORETICAL_C
VELOCITIES = np.flip(VELOCITIES)
VISCOSITIES = [1.0, 5.0, 10.0, 20.0, 50.0]

base_directory = "C:/Users/f51139dn/VISCOSITY_SWEEP_HOLEY_HYPERELASTIC/"


if not os.path.exists(base_directory):
    os.makedirs(base_directory)

for size in n_holes:
    for v in VISCOSITIES:
        for speed in VELOCITIES:
            os.chdir(base_directory)
            path = "h"+str(int(size)) + "_m" + str(int(1e5*MESH_SIZE)) + "_v" + str(int(100*v))+"_i" + str(int(100*speed)) + "_hyper"
            if os.path.exists(path):
                shutil.rmtree(path)
                os.makedirs(path)
            else:
                os.makedirs(path)
            output_dir = base_directory + path
            os.chdir(output_dir)
            # COMMONLY NEEDED VARIABLES
            MAX_STRAIN = 2 * speed / THEORETICAL_C
            
            # FOR DYNAMICAL ANALYSIS
            VELOCITY = speed
            EXCESS_TIME_FACTOR = 1.5
            
            # GEOMETRY OF PART
            NUM_HOLES = size
            
            DIAMETER = 8.0e-3
            HOLE_SEPARATION = 8.5e-3
            ENDS_LENGTH = 12.0e-3
            LIGAMENT_WIDTH = HOLE_SEPARATION - DIAMETER
            COLUMN_WIDTH = 2 * DIAMETER + 2 * LIGAMENT_WIDTH
            COLUMN_LENGTH = (2 * ENDS_LENGTH + NUM_HOLES * DIAMETER
                             + (NUM_HOLES - 1) * LIGAMENT_WIDTH)
            
            MAX_TIME = EXCESS_TIME_FACTOR * (COLUMN_LENGTH * MAX_STRAIN / VELOCITY)
            AMPLITUDE_TABLE = [[0, 0],
                               [MAX_TIME/EXCESS_TIME_FACTOR, 1]]
            
            SAVED_VARIABLES = ("S", "LE", "U")
            NUM_FIELD_OUTPUT_INTERVALS = 20
            
            
            # MODEL
            myModel = mdb.Model(name='Model-2d-flapped-column')
            
            OK_TO_RUN = True
            
            
            
            #---------------------PARTS-----------------------
            # PARTS
            def sketch_column():
                mySketch = myModel.ConstrainedSketch(name='Column Sketch',
                                                     sheetSize=COLUMN_LENGTH)
                left_x = -COLUMN_WIDTH / 2
                right_x = -left_x
                bottom_y = -COLUMN_LENGTH / 2
                mySketch.Line((left_x, bottom_y),
                              (right_x, bottom_y))
                mySketch.Line((left_x, -bottom_y),
                              (right_x, -bottom_y))
                mySketch.Line((left_x, bottom_y),
                              (left_x, bottom_y + ENDS_LENGTH))
                mySketch.Line((right_x, bottom_y),
                              (right_x, bottom_y + ENDS_LENGTH))
                mySketch.Line((left_x, -bottom_y),
                              (left_x, -bottom_y - ENDS_LENGTH))
                mySketch.Line((right_x, -bottom_y),
                              (right_x, -bottom_y - ENDS_LENGTH))
                
                circle_y_coordinates = np.array([])
                for n in range(NUM_HOLES):
                    center_x = 0
                    center_y = bottom_y + ENDS_LENGTH + DIAMETER / 2 + n * HOLE_SEPARATION
                    circle_y_coordinates = np.append(circle_y_coordinates, center_y)
                
                    mySketch.CircleByCenterPerimeter((center_x, center_y),
                                                     (center_x, center_y - DIAMETER / 2))
                    mySketch.ArcByCenterEnds((left_x, center_y),
                                             (left_x, center_y - DIAMETER / 2),
                                             (left_x, center_y + DIAMETER / 2))
                    mySketch.ArcByCenterEnds((right_x, center_y),
                                             (right_x, center_y + DIAMETER / 2),
                                             (right_x, center_y - DIAMETER / 2))
                    if n > 0:
                        mySketch.Line((left_x, center_y - DIAMETER / 2),
                                      (left_x, center_y - DIAMETER / 2 - LIGAMENT_WIDTH))
                        mySketch.Line((right_x, center_y - DIAMETER / 2),
                                      (right_x, center_y - DIAMETER / 2 - LIGAMENT_WIDTH))
                
                return mySketch, circle_y_coordinates
            
            column_sketch, circle_ys = sketch_column()
            part = myModel.Part(name='column', type=DEFORMABLE_BODY)
            part.BaseShell(sketch=column_sketch)
            
            # SETS
            end_coordinates = ((0, COLUMN_LENGTH / 2, 0),)
            ok_y_coord = COLUMN_LENGTH / 2 - ENDS_LENGTH / 2
            part.Set(name="set-all", faces=[part.faces.findAt(coordinates=((0, ok_y_coord, 0.0),))])
            part.Set(name="set-bottom", edges=part.edges.findAt(coordinates=((0, -COLUMN_LENGTH/2, 0),)))
            part.Set(name="set-top", edges=part.edges.findAt(coordinates=end_coordinates))
            part.Surface(name="surface-top", side1Edges=part.edges.findAt(coordinates=end_coordinates))
            
            # ASSEMBLY
            myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
            myModel.rootAssembly.Instance(name="part-instance", part=part, dependent=ON)
            
            # MATERIALS
            material_steel = myModel.Material("steel")
            material_steel.Density(table=((7500,),))
            material_steel.Elastic(table=((2.0e11, 0.3),))
            
            material_plastic = myModel.Material("hyperelastic")
            material_plastic.Density(table=((1160,),))
            material_plastic.Hyperelastic(((123.4e3, 19.23e3, 0),), type=MOONEY_RIVLIN,
                                          testData=OFF)
            
            # SECTION
            myModel.HomogeneousShellSection(name="part-section", material="hyperelastic",
                                            thickness=SECTION_THICKNESS)
            part.SectionAssignment(region=part.sets["set-all"], sectionName="part-section")
            
            
            
            #--------------------STEPS----------------------
            
            step = myModel.ExplicitDynamicsStep(
                name="Step-dynamic",
                previous="Initial",
                description="",
                nlgeom=True,
                timePeriod=MAX_TIME,
                linearBulkViscosity=v
                # quadBulkViscosity=0.12
            )
            
            
            
            #----------------------BCs------------------------
            
            front_rear_instance = myModel.rootAssembly.instances["part-instance"].sets["set-all"]
            front_rear_bc = myModel.DisplacementBC(name="BC-front-rear",
                                                   createStepName="Initial",
                                                   region=front_rear_instance, u3=SET)
            
            top_instance = myModel.rootAssembly.instances["part-instance"].sets["set-top"]
            top_displacement_bc = myModel.EncastreBC(name="BC-top",
                                                     createStepName="Initial",
                                                     region=top_instance)
            
            myModel.TabularAmplitude("Amp-displacement", data=tuple(AMPLITUDE_TABLE))
            
            bottom_instance = myModel.rootAssembly.instances["part-instance"].sets["set-bottom"]
            displacament_bc = myModel.DisplacementBC(name="displacement",
                                                     createStepName="Step-dynamic",
                                                     region=bottom_instance, u1=SET,
                                                     u2=MAX_STRAIN * COLUMN_LENGTH,
                                                     amplitude="Amp-displacement")
            
            
            
            #------------------INTERACTIONS-------------------
            
            normal_contact = myModel.ContactProperty("Normal-contact")
            general_contact = myModel.ContactExp("Self-Contact", createStepName="Step-dynamic")
            general_contact.contactPropertyAssignments.appendInStep("Step-dynamic",
                                                                    ((SELF, SELF, "Normal-contact"),))
            
            
            
            #--------------------MESHING AND NODE SETS----------------------
            # MESH
            part.setMeshControls((part.faces[0],), elemShape=TRI, technique=FREE)
            part.setElementType(part.sets["set-all"], (ElemType(elemCode=S3R,
                                                                elemLibrary=EXPLICIT, 
                                                                elemDeletion=ON),))
            part.seedPart(size=MESH_SIZE)
            part.generateMesh()
            
            # NODE SETS
            def define_node_set():
                all_nodes = part.nodes
                # Find the nodes for the first hole, then extrapolate to the other holes
                max_offset = MESH_SIZE
                center_x = 0
                center_z = 0
                
                center_y = circle_ys[0]
                top_nodes = []
                right_side_nodes = []
                left_side_nodes = []
                for node in all_nodes:
                    node_x = node.coordinates[0]
                    node_y = node.coordinates[1]
                    node_z = node.coordinates[2]
                
                    # TOP OF HOLE (i.e. -x-displaced node)
                    top_x = center_x - DIAMETER/2
                    if node_x < top_x + max_offset and node_x > top_x - max_offset:
                        if node_y < center_y + max_offset and node_y > center_y - max_offset:
                            if node_z < center_z + max_offset and node_z > center_z - max_offset:
                                top_nodes.append(node)
                
                    # LEFT OF HOLE (i.e. -y-displaced node)
                    left_y = center_y - DIAMETER/2
                    if node_x < center_x + max_offset and node_x > center_x - max_offset:
                        if node_y < left_y + max_offset and node_y > left_y - max_offset:
                            if node_z < center_z + max_offset and node_z > center_z - max_offset:
                                left_side_nodes.append(node)
                
                    # RIGHT OF HOLE (i.e. +y-displaced node)
                    right_y = center_y + DIAMETER/2
                    if node_x < center_x + max_offset and node_x > center_x - max_offset:
                        if node_y < right_y + max_offset and node_y > right_y - max_offset:
                            if node_z < center_z + max_offset and node_z > center_z - max_offset:
                                right_side_nodes.append(node)
                
                
                if len(top_nodes) > 1 or len(left_side_nodes) > 1:
                    top_center_coordinates = [center_x - DIAMETER/2, center_y, center_z]
                    left_center_coordinates = [center_x, center_y - DIAMETER/2, center_z]
                    top_required_offsets = np.array([])
                    for node in top_nodes:
                        node_coordinates = node.coordinates
                        max_offset = 0
                        for i in range(3):
                            offset = np.abs(top_center_coordinates[i] - node_coordinates[i])
                            if offset > max_offset:
                                max_offset = offset
                        top_required_offsets = np.append(top_required_offsets, max_offset)
                    left_required_offsets = np.array([])
                    for node in left_side_nodes:
                        node_coordinates = node.coordinates
                        max_offset = 0
                        for i in range(3):
                            offset = np.abs(left_center_coordinates[i] - node_coordinates[i])
                            if offset > max_offset:
                                max_offset = offset
                        left_required_offsets = np.append(left_required_offsets, max_offset)
                    left_required_offsets = np.sort(left_required_offsets)
                    top_required_offsets = np.sort(top_required_offsets)
                    max_offset = np.max([(top_required_offsets[0] + top_required_offsets[1]) / 2,
                                         (left_required_offsets[0] + left_required_offsets[1]) / 2])
                    max_offset += MESH_SIZE / 8
                
                
                node_set = []
                circle_numbers = []
                for circle_no, center_y in enumerate(circle_ys):
                    for node in all_nodes:
                        node_x = node.coordinates[0]
                        node_y = node.coordinates[1]
                        node_z = node.coordinates[2]
                
                        # TOP OF HOLE (i.e. -x-displaced node)
                        top_x = center_x - DIAMETER/2
                        if node_x < top_x + max_offset and node_x > top_x - max_offset:
                            if node_y < center_y + max_offset and node_y > center_y - max_offset:
                                if node_z < center_z + max_offset/2 and node_z > center_z - max_offset/2:
                                    node_set.append(node)
                                    circle_numbers.append(circle_no)
                
                        # BOTTOM OF HOLE (i.e. +x-displaced node)
                        top_x = center_x + DIAMETER/2
                        if node_x < top_x + max_offset and node_x > top_x - max_offset:
                            if node_y < center_y + max_offset and node_y > center_y - max_offset:
                                if node_z < center_z + max_offset/2 and node_z > center_z - max_offset/2:
                                    node_set.append(node)
                                    circle_numbers.append(circle_no)
                
                        # LEFT OF HOLE (i.e. -y-displaced node)
                        left_y = center_y - DIAMETER/2
                        if node_x < center_x + max_offset and node_x > center_x - max_offset:
                            if node_y < left_y + max_offset and node_y > left_y - max_offset:
                                if node_z < center_z + max_offset/2 and node_z > center_z - max_offset/2:
                                    node_set.append(node)
                                    circle_numbers.append(circle_no)
                
                        # RIGHT OF HOLE (i.e. +y-displaced node)
                        right_y = center_y + DIAMETER/2
                        if node_x < center_x + max_offset and node_x > center_x - max_offset:
                            if node_y < right_y + max_offset and node_y > right_y - max_offset:
                                if node_z < center_z + max_offset/2 and node_z > center_z - max_offset/2:
                                    node_set.append(node)
                                    circle_numbers.append(circle_no)
                print(str(len(node_set)) + " nodes in node set")
                print(str(int(len(node_set) / NUM_HOLES)) + " nodes per hole")
                
                return node_set, circle_numbers
            
            node_set, circle_numbers = define_node_set()
            node_array = mesh.MeshNodeArray(node_set)
            part.Set(nodes=node_array, name='set-nodes')
            
            # SAVE NODE LABELS AND COORDINATES
            node_labels = np.array([])
            node_xvals = np.array([])
            node_yvals = np.array([])
            for node in node_set:
                node_labels = np.append(node_labels, node.label)
                node_xvals = np.append(node_xvals, node.coordinates[0])
                node_yvals = np.append(node_yvals, node.coordinates[1])
                
            circle_numbers = np.array(circle_numbers)
            node_labels_path = "nodes_coords_circlenos_2d_m" + str(int(MESH_SIZE*100*1000)) +"_h" +str(int(size))+".csv"
            np.savetxt(node_labels_path,
                       np.column_stack((node_labels, node_xvals, node_yvals, circle_numbers)),
                       delimiter=",", header='node,x,y, circle no.')
            
            
            
            #--------------------OUTPUTS----------------------
            ### NOTE THAT ALL THE OUTPUTS ARE NOT NECESSARILY NEEDED
            
            # FIELD OUTPUTS
            field = myModel.FieldOutputRequest("F-Output-1", createStepName="Step-dynamic",
                                               variables=SAVED_VARIABLES,
                                               numIntervals=NUM_FIELD_OUTPUT_INTERVALS)
            
            # INTEGRATED OUTPUT SECTION
            # force_section = myModel.IntegratedOutputSection("I-Section-Force", surface=myModel.rootAssembly.\
            #                                                 instances["part-instance"].surfaces['surface-top'])
            
            # HISTORY OUTPUT REQUESTS
            # history_u = myModel.HistoryOutputRequest('H-Output-Displacement',
            #                                          'Step-dynamic',
            #                                          region=myModel.rootAssembly.instances["part-instance"].sets['set-bottom'],
            #                                          variables=['U2'])
            myModel.HistoryOutputRequest('H-Output-Nodes', 'Step-dynamic',
                                         region=myModel.rootAssembly.instances["part-instance"].sets['set-nodes'],
                                         variables=['U1', 'U2'])
            # myModel.HistoryOutputRequest('H-Output-Force', 'Step-dynamic',
            #                              integratedOutputSection="I-Section-Force", variables=['SOF'])
            
            
            
            #--------------------RUN JOB----------------------
            # OK_TO_RUN = False
            if(OK_TO_RUN):
                # JOB
                job = mdb.Job(name=path,
                              model="Model-2d-flapped-column",
                              numCpus=NUM_CPUS, numDomains=NUM_CPUS,
                              explicitPrecision=DOUBLE_PLUS_PACK)
                job.writeInput()
            
                # SUBMIT JOB
                job.submit()
                job.waitForCompletion()
            
            
            
            #--------------------EXTRACT RESULTS----------------------
            node_labels = np.genfromtxt(node_labels_path, dtype=int, delimiter=",")
            odb = session.openOdb(path + ".odb")
            for node_line in node_labels:
                node = node_line[0]
                region_key = 'Node PART-INSTANCE.' + str(node)
                data_u1 = np.array(odb.steps['Step-dynamic'].historyRegions[region_key]\
                                    .historyOutputs['U1'].data)
                np.savetxt(output_dir + '/results_m' + str(int(MESH_SIZE*100*1000))
                            + "_h"+str(int(NUM_HOLES)) + "_v" + str(int(100*v))
                            + "_speed"+ str(int(100*speed))
                            + '_node' + str(node) + '_u1.csv',
                            data_u1, header='time,U1', delimiter=',')
                
                data_u2 = np.array(odb.steps['Step-dynamic'].historyRegions[region_key]\
                                    .historyOutputs['U2'].data)
                np.savetxt(output_dir + '/results_m' + str(int(MESH_SIZE*100*1000))
                            + "_h"+str(int(NUM_HOLES)) + "_v" + str(int(100*v))
                            + "_speed"+ str(int(100*speed))
                            + '_node' + str(node) + '_u2.csv',
                            data_u2, header='time,U2', delimiter=',')
    
    
    
            #--------------------ANIMATION----------------------
            # ANIMATION
            executeOnCaeStartup()
            o2 = session.openOdb(name=path +".odb")
            session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=339.333343505859, 
                             height=208.144454956055)
            session.viewports['Viewport: 1'].setValues(displayedObject=o2)
            session.viewports['Viewport: 1'].makeCurrent()
            session.viewports['Viewport: 1'].maximize()
            session.viewports['Viewport: 1'].view.rotate(xAngle=0, yAngle=0, zAngle=-90, 
                                                         mode=TOTAL)
            session.viewports['Viewport: 1'].view.rotate(xAngle=0, yAngle=0, zAngle=-90, 
                                                         mode=TOTAL)
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
                DEFORMED, ))
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
                UNDEFORMED, ))
            session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
                CONTOURS_ON_DEF, ))
            session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
                visibleEdges=FREE)
            session.viewports['Viewport: 1'].animationController.setValues(
                animationType=TIME_HISTORY)
            session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
            session.viewports['Viewport: 1'].animationController.setValues(
                animationType=NONE)
            session.viewports['Viewport: 1'].view.fitView()
            session.viewports['Viewport: 1'].animationController.setValues(
                animationType=TIME_HISTORY)
            session.viewports['Viewport: 1'].animationController.play(duration=UNLIMITED)
            session.graphicsOptions.setValues(backgroundStyle=SOLID)
            session.graphicsOptions.setValues(backgroundColor='#FFFFFF')
            session.viewports['Viewport: 1'].view.fitView()
            session.viewports['Viewport: 1'].view.fitView()
            session.imageAnimationOptions.setValues(vpDecorations=ON, vpBackground=OFF, 
                                                    compass=OFF)
            session.writeImageAnimation(fileName=output_dir + '/' + path + '-animation', 
                                        format=AVI, canvasObjects=(session.viewports['Viewport: 1'], ))
            session.viewports['Viewport: 1'].minimize()
    
            print("DONE!")

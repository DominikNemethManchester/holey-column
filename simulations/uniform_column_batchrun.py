# -*- coding: utf-8 -*-
"""
Draws UNIFORM column with same size as an N-hole flapped column.

Also does the usual setting up of the materials, sections, assembly, etc.
"""

import numpy as np
import os
import shutil
from abaqus import *
from abaqusConstants import *
from driverUtils import executeOnCaeStartup

import mesh
from mesh import ElemType
import regionToolset


# COMMONLY NEEDED VARIABLES
MESH_SIZE = 0.2
MAX_STRAIN = 0.15
NUM_CPUS = 8
SECTION_THICKNESS = 1.0


E = 6*(123.4e3 + 19.23e3)
DENSITY = 1160
THEORETICAL_C = np.sqrt(E/DENSITY)

# FOR DYNAMICAL ANALYSIS
EXCESS_TIME_FACTOR = 1

# GEOMETRY OF PART
NUM_HOLES = 10

DIAMETER = 8.0
HOLE_SEPARATION = 8.5
ENDS_LENGTH = 12.0
LIGAMENT_WIDTH = HOLE_SEPARATION - DIAMETER
COLUMN_WIDTH = 2 * DIAMETER + 2 * LIGAMENT_WIDTH
COLUMN_LENGTH = (2 * ENDS_LENGTH + NUM_HOLES * DIAMETER
                 + (NUM_HOLES - 1) * LIGAMENT_WIDTH)

SAVED_VARIABLES = ("E", "U")
NUM_FIELD_OUTPUT_INTERVALS = 10

# SIMULATION PARAMETERS
VELOCITIES = np.array([0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15]) * THEORETICAL_C
VISCOSITIES = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]
VELOCITIES = np.flip(VELOCITIES)


base_directory = "C:/Users/f51139dn/VISCOSITY_SWEEP_UNIFORM_HYPERELASTIC_PULSE2/"


# MODEL
myModel = mdb.Model(name='Model-2d-uniform-column')


#---------------------PARTS-----------------------

# PARTS
def sketch_column():
    mySketch = myModel.ConstrainedSketch(name='Column Sketch',
                                          sheetSize=COLUMN_LENGTH)
    left_x = -COLUMN_WIDTH / 2
    bottom_y = -COLUMN_LENGTH / 2
    mySketch.rectangle((left_x, bottom_y), (-left_x, -bottom_y))
    
    circle_y_coordinates = np.array([])
    for n in range(NUM_HOLES):
        center_y = bottom_y + ENDS_LENGTH + DIAMETER / 2 + n * HOLE_SEPARATION
        circle_y_coordinates = np.append(circle_y_coordinates, center_y)
        circle_y_coordinates = np.append(circle_y_coordinates, center_y + HOLE_SEPARATION / 2)
    
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
material_steel.Elastic(table=((2.0e11, 0),))

material_plastic = myModel.Material("hyperelastic")
material_plastic.Density(table=((1160,),))
material_plastic.Hyperelastic(((123.4e3, 19.23e3, 0),), type=MOONEY_RIVLIN,
                              testData=OFF)

# SECTION
myModel.HomogeneousShellSection(name="part-section", material="hyperelastic",
                                thickness=SECTION_THICKNESS)
part.SectionAssignment(region=part.sets["set-all"], sectionName="part-section")



#--------------------MESHING----------------------
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
    centre_nodes = []
    for node in all_nodes:
        node_x = node.coordinates[0] + max_offset/2
        node_y = node.coordinates[1]
        node_z = node.coordinates[2]
    
        # CENTRE OF HOLE
        if node_x < center_x + max_offset and node_x > center_x - max_offset:
            if node_y < center_y + max_offset and node_y > center_y - max_offset:
                if node_z < center_z + max_offset and node_z > center_z - max_offset:
                    centre_nodes.append(node)
    
    
    if len(centre_nodes) > 1:
        top_center_coordinates = [center_x - max_offset/2, center_y, center_z]
        top_required_offsets = np.array([])
        for node in centre_nodes:
            node_coordinates = node.coordinates
            max_offset = 0
            for i in range(3):
                offset = np.abs(top_center_coordinates[i] - node_coordinates[i])
                if offset > max_offset:
                    max_offset = offset
            top_required_offsets = np.append(top_required_offsets, max_offset)
        
        top_required_offsets = np.sort(top_required_offsets)
        max_offset = (top_required_offsets[0] + top_required_offsets[1]) / 2
        max_offset += MESH_SIZE / 8
    
    
    node_set = []
    circle_numbers = []
    for circle_no, center_y in enumerate(circle_ys):
        for node in all_nodes:
            node_x = node.coordinates[0] + max_offset/2
            node_y = node.coordinates[1]
            node_z = node.coordinates[2]
    
            # TOP OF HOLE (i.e. -x-displaced node)
            if node_x < center_x + max_offset and node_x > center_x - max_offset:
                if node_y < center_y + max_offset and node_y > center_y - max_offset:
                    if node_z < center_z + max_offset and node_z > center_z - max_offset:
                        node_set.append(node)
                        circle_numbers.append(circle_no)
    
    print(str(len(node_set)) + " nodes in node set")
    print(str(int(len(node_set) / NUM_HOLES)) + " nodes per hole")
    
    return node_set, circle_numbers


node_set, circle_numbers = define_node_set()

node_labels = np.array([])
node_xvals = np.array([])
node_yvals = np.array([])
for node in node_set:
    node_labels = np.append(node_labels, node.label)
    node_xvals = np.append(node_xvals, node.coordinates[0])
    node_yvals = np.append(node_yvals, node.coordinates[1])

circle_numbers = np.array(circle_numbers)
node_labels_path = (base_directory + "UNIFORM_nodes_coords_circlenos_2d_m"
                    + str(int(MESH_SIZE*100)) + "_h" + str(int(NUM_HOLES)) + ".csv")
np.savetxt(node_labels_path,
            np.column_stack((node_labels, node_xvals, node_yvals, circle_numbers)),
            delimiter=",", header='node, x, y, circle no.')

node_array = mesh.MeshNodeArray(node_set)
part.Set(nodes=node_array, name='set-nodes')
    
   
for velocity in VELOCITIES:
    for viscosity in VISCOSITIES:
        MAX_STRAIN = 2 * velocity / THEORETICAL_C
        # if velocity < 0.02 * 5164:
        #     MAX_STRAIN = 0.02
        # elif velocity < 0.9 * 5164:
        #     MAX_STRAIN = 0.05
        # elif velocity > 19 * 5164:
        #     MAX_STRAIN = 0.4
        # else:
        #     MAX_STRAIN = 0.25
        path = base_directory + "m" + str(int(100*MESH_SIZE)) + "v" + str(int(viscosity*100)) +"i"+ str(int(velocity*100)) + "_HYPERELASTIC"
        if os.path.exists(path):
            shutil.rmtree(path)
            os.makedirs(path)
        elif os.path.exists(path)==False:
            os.makedirs(path)
        os.chdir(path)
        output_dir = path
        OK_TO_RUN = True
        
        MAX_TIME = EXCESS_TIME_FACTOR * (COLUMN_LENGTH * MAX_STRAIN / velocity)
        AMPLITUDE_TABLE = [[0, 0],
                           [MAX_TIME / EXCESS_TIME_FACTOR*0.1, 1*0.1],
                           [MAX_TIME / EXCESS_TIME_FACTOR*0.2, 1*0.1],
                           [MAX_TIME / EXCESS_TIME_FACTOR*0.3, 0]]
        
        #--------------------STEPS----------------------
        
        step = myModel.ExplicitDynamicsStep(
            name="Step-dynamic",
            previous="Initial",
            description="",
            nlgeom=True,
            timePeriod=MAX_TIME,
            linearBulkViscosity=viscosity
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
        
        myModel.TabularAmplitude("Amp-displacement", data=tuple(AMPLITUDE_TABLE), smooth=SOLVER_DEFAULT)
        
        bottom_instance = myModel.rootAssembly.instances["part-instance"].sets["set-bottom"]
        displacament_bc = myModel.DisplacementBC(name="displacement",
                                                 createStepName="Step-dynamic",
                                                 region=bottom_instance, u1=SET,
                                                 u2=MAX_STRAIN * COLUMN_LENGTH,
                                                 amplitude="Amp-displacement")
        
        
        
        #------------------INTERACTIONS-------------------
        # not necessary since uniform column
        
        
        #--------------------OUTPUTS----------------------
        # OUTPUT REQUEST
        field = myModel.FieldOutputRequest("F-Output-1", createStepName="Step-dynamic",
                                           variables=SAVED_VARIABLES,
                                           numIntervals=NUM_FIELD_OUTPUT_INTERVALS)
        field = myModel.FieldOutputRequest("F-Output-1", createStepName="Step-dynamic",
                                           variables=('S',),
                                           numIntervals=100)
        
        # INTEGRATED OUTPUT SECTION
        force_section = myModel.IntegratedOutputSection("I-Section-Force", surface=myModel.rootAssembly.\
                                                        instances["part-instance"].surfaces['surface-top'])
        
        # HISTORY OUTPUT REQUESTS
        history_u = myModel.HistoryOutputRequest('H-Output-Displacement',
                                                 'Step-dynamic',
                                                 region=myModel.rootAssembly.instances["part-instance"].sets['set-bottom'],
                                                 variables=['U2'])
        myModel.HistoryOutputRequest('H-Output-Nodes', 'Step-dynamic',
                                     region=myModel.rootAssembly.instances["part-instance"].sets['set-nodes'],
                                     variables=['U1', 'U2', 'S22'])
        myModel.HistoryOutputRequest('H-Output-Force', 'Step-dynamic',
                                     integratedOutputSection="I-Section-Force", variables=['SOF'])
        
        
        
        if(OK_TO_RUN):
            # JOB
            job_name = "Job-h"+str(int(NUM_HOLES)) + "-i" + str(int(velocity*100)) + "-uniform"
            job = mdb.Job(name=job_name,
                          model="Model-2d-uniform-column",
                          numCpus=NUM_CPUS, numDomains=NUM_CPUS,
                          explicitPrecision=DOUBLE_PLUS_PACK)
            job.writeInput()
            job.submit()
            job.waitForCompletion()
            
            odb = session.openOdb(job_name + ".odb")
            
            node_labels = np.genfromtxt(node_labels_path, dtype=int, delimiter=",")
            for node_line in node_labels:
                node = node_line[0]
                region_key = 'Node PART-INSTANCE.' + str(node)            
                data_u2 = np.array(odb.steps['Step-dynamic'].historyRegions[region_key]\
                                    .historyOutputs['U2'].data)
                np.savetxt(output_dir + '/results_m' + str(int(MESH_SIZE*100))
                            + "_h"+str(int(NUM_HOLES)) + "_v" + str(int(100*viscosity))
                            + "_speed"+ str(int(100*velocity))
                            + '_node' + str(node) + '_u2.csv',
                            data_u2, header='time,U2', delimiter=',')
                
            # ANIMATION
            session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=339.333343505859, 
                height=208.144454956055)
            session.viewports['Viewport: 1'].makeCurrent()
            session.viewports['Viewport: 1'].maximize()
            
            executeOnCaeStartup()
            session.viewports['Viewport: 1'].setValues(displayedObject=odb)
            session.viewports['Viewport: 1'].makeCurrent()
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
            session.writeImageAnimation(
                fileName=path+'/'+job_name+'-animation', 
                format=AVI, canvasObjects=(session.viewports['Viewport: 1'], ))
            # END ANIMATION
            
            os.chdir(base_directory)

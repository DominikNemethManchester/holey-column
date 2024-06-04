# holey-column

## Running the codes
Any code with ```import abaqus ``` in its header needs to be run from Abaqus CAE to execute, using ```File > Run Script```.

## Compression Wave Analysis
1. Run sound_speed_by_displacement_batch_v2.py
   Determines the compression wave arrival times according to the final method described in our report (by finding the inflection point of the U2 vs time curve).
2. Run compression_wave_analysis_v2.py
   Determines compression wave speeds and produces figures by plotting the node coordinates and the corresponding compression wave arrival times.


## Buckling and Self-Contact Wave Analysis
  Run bucklingwave_v2.py
  Calculates wave speeds and produces figures by plotting the node coordinates and the corresponding wave arrival times. Also produces figures that compare all three wave phenomena (in terms of viscosity and indenter speed)


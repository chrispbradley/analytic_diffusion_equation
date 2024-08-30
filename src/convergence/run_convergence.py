#!/usr/bin/env python

import sys,os

if len(sys.argv) > 1:
    if len(sys.argv) > 2:
        sys.exit('ERROR: too many arguments- currently only accepting up to 1 options: executablename')
    executable = sys.argv[1]
else:
    sys.exit('ERROR: must supply the exectuable name as an option.')

numberOfElements = (10,20,40,80)
interpolationTypes = (1,2,3,4)
timeSteps = (0.10,0.05,0.01,0.005)
dynamicDegrees = (1,2)

for dynamicDegree in dynamicDegrees:
    for timeStep in timeSteps:
        for interpolationType in interpolationTypes:
            for numberOfElement in numberOfElements:            
                command = "%s %d %d %f %i" % (executable,numberOfElement,interpolationType,timeStep,dynamicDegree)
                os.system(command)
                

import sys
import re 
import os 
import glob

import numpy as np
import math 

from scipy.spatial.distance import cdist



if 'RZ_PYTHON_TOOLS' not in os.environ :
    print('environment variable RZ_PYTHON_TOOLS not set')
    sys.exit(1)

sys.path.append(os.environ['RZ_PYTHON_TOOLS'])

artemisweb = os.environ['ARTEMISWEB']

sys.path.append(f'{artemisweb}/artemisdiscovery-server/pyscripts')

from gromacs_tools import *

from artemisserver_tools import *

from tools import dihedralForCoords

from copy import deepcopy


sys.path.append(os.environ['RZ_PYTHON_TOOLS'])

artemisweb = os.environ['ARTEMISWEB']

sys.path.append(f'{artemisweb}/artemisdiscovery-server/pyscripts')


import shutil 
import gc

import json
import pickle
import codecs

path = './settings.txt'
settings = importSettings(path)
usersettings = settings 

path = './baseDATA.0.base64.txt'

fH =  open(path, 'rb')
txt = fH.read()
d = codecs.decode(txt, "base64")
coordsource = pickle.loads(d)

coords = np.array([ x[6] for x in coordsource['atomsInORDER']])
radiiInORDER = coordsource['radiiInORDER']

X = coords[:,0]
Y = coords[:,1]
Z = coords[:,2]
#
probeRad = usersettings['proberad']
xMIN = np.min(X - radiiInORDER - 3.*probeRad)
xMAX = np.max(X + radiiInORDER + 3.*probeRad)
yMIN = np.min(Y - radiiInORDER - 3.*probeRad)
yMAX = np.max(Y + radiiInORDER + 3.*probeRad)
zMIN = np.min(Z - radiiInORDER - 3.*probeRad)
zMAX = np.max(Z + radiiInORDER + 3.*probeRad)
#
center = np.array([np.mean(X),np.mean(Y),np.mean(Z)])
#
# adjusted by atom and probe radius inside the function now
limits = ((xMIN,xMAX),(yMIN,yMAX),(zMIN,zMAX))

plane = 'X'
planeIDX = {'X':0, 'Y':1, 'Z':2}[plane]

levelSpacing = usersettings['level_spacing']
minLevel = round(limits[planeIDX][0]/levelSpacing)*levelSpacing - levelSpacing
maxLevel = round(limits[planeIDX][1]/levelSpacing)*levelSpacing + levelSpacing
#
level = minLevel
levels = []
circlesForLevels = [] 
contoursForLevels = []
#
while level <= maxLevel :
    levels.append(level)
    #
    contours = []
    #
    circles = atomCirclesForLevel(plane, level, coords, radiiInORDER, probeRad )
    #
    if len(circles) == 0 : 
        circlesForLevels.append(circles)
        contoursForLevels.append(contours)
        level += levelSpacing
        continue
    #
    circlesForLevels.append(circles)
    #
    edges, circles, circleNeighbors, components = neighborsForCircles(circles)
    #
    sing_contours = []
    nonsing_contours = []
    #
    for component in components :
        if len(component) == 1 :
            # add singleton
            sing_contours.append([ Arc(circles[list(component)[0]], None, None, None, None, singleton=True) ])
    #
    nonsing_contours = Arc.contoursFromCirclePairs( circles, edges, circleNeighbors )
    #
    #contours = sing_contours
    contours += [c for c in nonsing_contours if contourClockwise(c)]
    contours += [c for c in sing_contours ]
    #
    contoursForLevels.append(contours)
    level += levelSpacing

# write out data

fh = open('data.txt', 'w')

x = fh.write('# coordinates\n')
x = fh.write(f'{coords.shape[0]}\n')

for c in coords :
    x = fh.write('%12.8f\t%12.8f\t%12.8f\n' % (c[0],c[1],c[2]))

x = fh.write('# radii\n')
x = fh.write(f'{radiiInORDER.shape[0]}\n')

for r in radiiInORDER :
    x = fh.write('%12.8f\n' % r )

x = fh.write('# probeRad\n')
x = fh.write('%12.8f\n' % probeRad )

x = fh.write('# atom circles\n')

x = fh.write('%d\t#nlevels\n' % len(levels))

for lidx in range(len(levels)) :
    x = fh.write('%12.8f\t%d\t#level,ncircles\n' % (levels[lidx],len(circlesForLevels[lidx])))
    for circle in circlesForLevels[lidx] :
        x = fh.write('%d\t#atom\n' % circle.atoms[0])
        x = fh.write('%12.8f\t%12.8f\t%12.8f\t%12.8f\t%d\t#center,radius,index\n' % 
                     (circle.center[0],circle.center[1],circle.center[2],circle.radius,circle.index))


x = fh.write('# contours\n')
x = fh.write('%d\t#nlevels\n' % len(levels))

for lidx in range(len(levels)) :
    x = fh.write('%12.8f\t%d\t#level,ncontours\n' % (levels[lidx],len(contoursForLevels[lidx])))
    for contour in contoursForLevels[lidx] :
        x = fh.write('%d\t#n arcs next contour\n' % len(contour))
        for arc in contour :
            if arc.singleton :
                x = fh.write('%d\t#singletonarc,circleidx\n' % (arc.circle.index))
            else :
                x = fh.write('%d\t#nonsingletonarc,circleidx\n' % (arc.circle.index))
                x = fh.write('%12.8f\t%12.8f\t%12.8f\t%d\t#intersectI,circleIidx\n' % 
                    (arc.intersectI[0],arc.intersectI[1],arc.intersectI[2],arc.neighborI.index))
                x = fh.write('%12.8f\t%12.8f\t%12.8f\t%d\t#intersectJ,circleJidx\n' % 
                    (arc.intersectJ[0],arc.intersectJ[1],arc.intersectJ[2],arc.neighborJ.index))

fh.close() 


            





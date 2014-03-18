#!/usr/bin/env python

import sys, os
from math import *

def mathListToPyList(s):
    s = s.replace('{', '[').replace('}', ']')
    s = s.replace('*^', 'E')
    try:
        l = eval(s)
    except:
        return []

    return l

def computeNFlux(wTilde, s2, m, astar, greybody):
    if greybody == 0: return 0

    val = 4*pi*((1+astar**2)*wTilde - m*astar)/((nDim-4+1)+(nDim-4-1)*astar**2)
    ## Avoid overflow error when val > 709.78
    ## We can safely assume nFlux = 0 nevertherless of other factors
    if val > 709.78: return 0

    val = exp(val)
    if s2 == 1: val += 1
    else: val -= 1

    ## By definition, greybody factor = 0 when thermal factor is 0
    if val == 0: return 0;

    return max(0, greybody/val)

def findDataFiles(d):
    l = []
    for fileName in os.listdir(d):
        filePath = os.path.join(d, fileName)
        if os.path.isdir(filePath):
            l.extend(findDataFiles(filePath))
        elif 'greybody_D5_s' in fileName:
            l.append(os.path.join(d, fileName))
    return l

#spinStrToS2 = {"s0":0, "s12":1, "s1":2}

#srcDir = "/users/jhgoh/Dropbox/BH_SKKU/greybody code/5D_calculation"
dataFiles  = findDataFiles(os.path.expanduser("~/Dropbox/BH_SKKU/fast_s1s2"))
dataFiles += findDataFiles(os.path.expanduser("~/Dropbox/BH_SKKU/fast_s0"))
#dataFiles = findDataFiles(os.path.expanduser("~/Dropbox/BH_SKKU/greybody code/example_Hyun/fast_140314"))

cNFluxData = {}

for filePath in dataFiles:
    fileName = os.path.basename(filePath)
    print "Processing", fileName
    modeStr = fileName[len('greybody_D5_'):].split('_')

    s2 = int(modeStr[0][1])#spinStrToS2[modeStr[0]]
    a10  = int(modeStr[1][1:])

    contents = open(filePath).read()
    contents = contents.replace('\n', '').replace('\r', '').strip()
    for line in contents.replace('}}', '}}\n').split('greybodyTable'):
        line = line.strip()
        if line == "": continue
        mode, table = line.split(' = ')
        nDim, s, l, m, a = eval(mode)
        l2, m2 = l*2, m*2
        table = mathListToPyList(table)

        nFluxes = []
        for i in range(len(table)):
            x, y = table[i]
            nFlux = computeNFlux(x, s2, m2, a10, y)
            nFluxes.append((x, nFlux))

        cNFluxes = [(0,0)]
        for i in range(1, len(nFluxes)):
            x1, y1 = nFluxes[i-1]
            x2, y2 = nFluxes[i]
            area = (y2+y1)/2*(x2-x1)
            cNFluxes.append((x2, cNFluxes[-1][1]+area))

        cNFluxData[(nDim, s2, l2, m2, a10)] = cNFluxes

## Store all results into files
outDir = "."
f = open("%s/D%d/cFlux.dat" % (outDir, 5), "w")
print>>f, "#Cumulative flux data tables"
print>>f, "#I nDim s2 l2 m2 a10"
print>>f, "#X x1 x2 x3 x4 ...."
print>>f, "#Y y1 y2 y3 y4 ...."
for key in cNFluxData.keys():
    nDim, s2, l2, m2, a10 = key
    cNFlux = cNFluxData[key]
    print>>f, "I %d %d %d %d %d" % (nDim, s2, l2, m2, a10)
    print>>f, ("X "+(" ".join(["%13.9e" % x for x, y in cNFlux])))
    print>>f, ("Y "+(" ".join(["%13.9e" % y for x, y in cNFlux])))

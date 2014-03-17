#!/usr/bin/env python

import sys, os
from math import *

def mathListToPyList(s):
    s = s.replace('{', '[').replace('}', ']')
    s = s.replace('*^', 'E')
    l = eval(s)

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
dataFiles  = findDataFiles("/users/jhgoh/Dropbox/BH_SKKU/fast_s1s2")
dataFiles += findDataFiles("/users/jhgoh/Dropbox/BH_SKKU/fast_s0")

cNFluxData = {}

for filePath in dataFiles:
    fileName = os.path.basename(filePath)
    print "Processing", fileName
    modeStr = fileName[len('greybody_D5_'):-1].split('_')

    s2 = int(modeStr[0][1])#spinStrToS2[modeStr[0]]
    a10  = float(modeStr[1][1:])

    contents = open(filePath).read()
    contents = contents.replace('\n', '').replace('\r', '').strip()
    for line in contents.replace('}}', '}}\n').split('greybodyTable'):
        line = line.strip()
        if line == "": continue
        mode, table = line.split(' = ')
        nDim, s, l, m, a = eval(mode)
        table = mathListToPyList(table)

        nFluxes = []
        for i in range(len(table)):
            x, y = table[i]
            nFlux = computeNFlux(x, s2, m, a, y)
            nFluxes.append((x, nFlux))

        cNFluxes = [(0,0)]
        for i in range(1, len(nFluxes)):
            x1, y1 = nFluxes[i-1]
            x2, y2 = nFluxes[i]
            area = (y2+y1)/2*(x2-x1)
            cNFluxes.append((x2, cNFluxes[-1][1]+area))

        cNFluxData[(nDim, s, l, m, a)] = cNFluxes

## Store all results into files
outDir = "../data/flux"
fouts = []
for s2 in range(3):
    f = open("%s/D%d/cFlux_ss%d.dat" % s2, "w")
    fouts.append(f)
for key in cNFluxData.keys():
    nDim, s2, l, m, a = key
    cNFlux = cNFluxData[key]
    fout = open("%s/D%d/ss%1d/l%02d_m%03d_aa%02d.dat" % (outDir, nDim, s2, l, m, a), "w")

#!/usr/bin/env python

from ROOT import *
from array import array

gROOT.ProcessLine(".x ../rootlogon.C")

p1, p2 = 0., 8.
p = 4.

def findNearestIndex(x, vector):
    if x == vector[-1]: return len(vector)-1
    for i in xrange(len(vector)-1):
        if vector[i] <= x < vector[i+1]: return i
    return -1

def interpolate(x, xValues, yValues):
    i = findNearestIndex(x, xValues)
    if i < 0 or i == len(xValues)-1: i=len(xValues)-2
    x1, x2 = xValues[i:i+2]
    y1, y2 = yValues[i:i+2]
    if x1 == x2: return y1
    slope = (y2-y1)/(x2-x1)
    return slope*(x-x1) + y1

def differentiate(cgrp):
    grp = TGraph()
    grp.SetMarkerColor(cgrp.GetMarkerColor())
    grp.SetLineColor(cgrp.GetLineColor())
    grp.SetPoint(0, 0, 0)
    for i in xrange(0, cgrp.GetN()-1):
        x1 = cgrp.GetX()[i]
        x2 = cgrp.GetX()[i+1]
        y1 = cgrp.GetY()[i]
        y2 = cgrp.GetY()[i+1]
        dx = x2-x1
        dy = y2-y1
        grp.SetPoint(grp.GetN(), x2, dy/dx)
    return grp

print "Loading data table..."
dataTable = {}
mode, xValues, yValues = None, None, None
for line in open("../../../data/flux/cFlux_D5.dat").readlines():
    line = line.strip()
    if len(line) == 0: continue
    if line[0] == '#': continue

    key = line[0]
    values = eval(line[2:].replace(' ', ','))

    if line[0] == 'I':
        mode = values
        xValues = []
        yValues = []
    elif line[0] == 'X':
        xValues = array('d', values)
    elif line[0] == 'Y':
        yValues = array('d', values)
        dataTable[tuple(mode)] = (xValues, yValues)

print "Selecting test dataset..."
data1 = dataTable[(5, 0, 0, 0, p1)]
data2 = dataTable[(5, 0, 0, 0, p2)]
nData1 = len(data1[0])
nData2 = len(data2[0])
maxY1 = max(data1[1])
maxY2 = max(data2[1])

## Reference for direct comparison
data3 = dataTable[(5, 0, 0, 0, p)]
maxY3 = max(data3[1])
normY3 = [y/maxY3 for y in data3[1]]
grp3 = TGraph(len(data3[0]), array('d', data3[0]), array('d', normY3))

## Do the morphing
print "Starting to morph..."
normY1 = [y/maxY1 for y in data1[1]]
normY2 = [y/maxY2 for y in data2[1]]
newY = sorted(list(set(normY1+normY2)))
newX = []
for y in newY:
    x1 = interpolate(y, normY1, data1[0])
    x2 = interpolate(y, normY2, data2[0])
    x = interpolate(p, [p1, p2], [x1, x2])
    newX.append(x)

## Similar thing with random numbers
hInterp = TH1F("hInterp", "hInterp", 1000, 0, 1)
for i in xrange(1000000):
    y = gRandom.Uniform(0, 1)
    x1 = interpolate(y, normY1, data1[0])
    x2 = interpolate(y, normY2, data2[0])
    x = interpolate(p, [p1, p2], [x1, x2])
    hInterp.Fill(x)
hInterp.Scale(1./hInterp.Integral()/hInterp.GetBinWidth(1))

grp1 = TGraph(nData1, array('d', data1[0]), array('d', normY1))
grp2 = TGraph(nData2, array('d', data2[0]), array('d', normY2))
grp = TGraph(len(newX), array('d', newX), array('d', newY))
grp1.SetLineColor(kRed)
grp1.SetMarkerColor(kRed)
grp2.SetLineColor(kBlue)
grp2.SetMarkerColor(kBlue)
grp.SetLineColor(kGreen+1)
grp.SetMarkerColor(kGreen+1)

c = TCanvas("c", "c", 500, 500)
hFrameCDF = TH1F("hFrameCDF", "hFrame", 100, 0, 0.8)#max(max(data1[0]), max(data2[0])))
hFrameCDF.SetMinimum(0)
#hFrameCDF.SetMaximum(1.1*max(max(data1[1]), max(data2[1])))
hFrameCDF.SetMaximum(1.1)
hFrameCDF.Draw()
grp1.Draw("LP")
grp2.Draw("LP")
grp3.Draw("LP")
grp.Draw("LP")
legCDF = TLegend(0.5, 0.2, 0.9, 0.4)
legCDF.SetBorderSize(0)
legCDF.SetFillStyle(0)
legCDF.AddEntry(grp1, "a10=%d" % p1, "lp")
legCDF.AddEntry(grp2, "a10=%d" % p2, "lp")
legCDF.AddEntry(grp3, "a10=%d" % p, "lp")
legCDF.AddEntry(grp, "Interpolated: a10=%d" % p, "lp")
legCDF.Draw()

cDiff = TCanvas("cDiff", "cDiff", 500, 500)
hFrameDiff = TH1F("hFrameDiff", "hFrameDiff", 100, 0, 0.8)
hFrameDiff.SetMinimum(0)
hFrameDiff.SetMaximum(5)
hFrameDiff.Draw()
grpDiff1 = differentiate(grp1)
grpDiff2 = differentiate(grp2)
grpDiff = differentiate(grp)
grpDiff1.Draw("LP")
grpDiff2.Draw("LP")
hInterp.Draw("same")

grpDiff3 = differentiate(grp3)
grpDiff3.Draw("LP")
grpDiff.Draw("LP")

legDiff = TLegend(0.5, 0.7, 0.9, 0.9)
legDiff.SetBorderSize(0)
legDiff.SetFillStyle(0)
legDiff.AddEntry(grp1, "a10=%d" % p1, "lp")
legDiff.AddEntry(grp2, "a10=%d" % p2, "lp")
legDiff.AddEntry(grp3, "a10=%d" % p, "lp")
legDiff.AddEntry(grp, "Interpolated: a10=%d" % p, "lp")
legDiff.Draw()



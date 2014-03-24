#!/usr/bin/env python

import sys, os
from array import array
from math import *
from ROOT import *
gROOT.ProcessLine(".x ../rootlogon.C")
gStyle.SetPadRightMargin(0.1);


def mathListToPyList(s):
    s = s.replace('{', '[').replace('}', ']')
    s = s.replace('*^', 'E')

    try:
        l = eval(s)
    except:
        return []

    return l

def computeNFlux(wTilde, s2, m, astar, greybody):
    if abs(greybody) <= 1e-7: return 0

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
        if fileName == 'backup': continue
        filePath = os.path.join(d, fileName)
        if os.path.isdir(filePath):
            l.extend(findDataFiles(filePath))
        elif 'greybody_D5_s' in fileName:
            l.append(os.path.join(d, fileName))
    return l

c = TCanvas("c", "c", 900, 600)
c.Divide(3,2)
hFrameG, hFrameF = [], []
for s2 in range(3):
    h = TH1F("hFrameG%d" % s2, "hFrameG;#tilde#omega;Greybody factor", 100, 0, 20)
    h.SetMinimum(-1)
    c.cd(s2+1)
    h.Draw()
    hFrameG.append(h)

    h = TH1F("hFrameF%d" % s2, "hFrameF;#tilde#omega;Number flux", 100, 0, 20)
    h.SetMinimum(0)
    c.cd(s2+4)
    h.Draw()
    hFrameF.append(h)

colors = [kRed, kOrange, kGreen, kGreen+1, kAzure+1, kBlue, kBlue+1]

bins = array('d', [0.1*(x-1) for x in range(0, 20, 2)])
lbins = array('d', range(0, 30))
hSumNFlux = [
    TH1F("hSumNFlux0", "Scalar;a*;Number flux", len(bins)-1, bins),
    TH1F("hSumNFlux1", "Spinor;a*;Number flux", len(bins)-1, bins),
    TH1F("hSumNFlux2", "Vector;a*;Number flux", len(bins)-1, bins),
]
hSumNFlux[0].SetLineColor(kRed)
hSumNFlux[1].SetLineColor(kBlue)
hSumNFlux[2].SetLineColor(kGreen+1)

grpGs = [[], [], []]
grpFs = [[], [], []]
hF2Ds = [
    TH2F("hF2D0", "Scalar;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
    TH2F("hF2D1", "Spinor;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
    TH2F("hF2D2", "Vector;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
]

hNDOFs = [
    TH2F("hNDOF0", "Scalar;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
    TH2F("hNDOF1", "Spinor;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
    TH2F("hNDOF2", "Vector;a*;floor(l mode)", len(bins)-1, bins, len(lbins)-1, lbins),
]
#spinStrToS2 = {"s0":0, "s12":1, "s1":2}

for s2 in range(3):
    hF2Ds[s2].SetName("hF2D_ss%d" % s2)
    hF2Ds[s2].SetTitle("s2=%d" % s2)
    hNDOFs[s2].SetName("hNDOF_ss%d" % s2)
    hNDOFs[s2].SetTitle("s2=%d" % s2)

#srcDir = "/users/jhgoh/Dropbox/BH_SKKU/greybody code/5D_calculation"
dataFiles  = findDataFiles("/users/jhgoh/Dropbox/BH_SKKU/fast_s1s2")
dataFiles += findDataFiles("/users/jhgoh/Dropbox/BH_SKKU/fast_s0")

for filePath in dataFiles:
    fileName = os.path.basename(filePath)
    print "Processing", fileName
    modeStr = fileName[len('greybody_D5_'):].split('_')

    s2 = int(modeStr[0][1])#spinStrToS2[modeStr[0]]
    a10  = float(modeStr[1][1:])

    contents = open(filePath).read()
    contents = contents.replace('\\', '').replace('\n', '').replace('\r', '').strip()
    for line in contents.replace('}}', '}}\n').split('greybodyTable'):
        line = line.strip()
        if line == "": continue
        mode, table = line.split(' = ')
        nDim, s, l, m, a = eval(mode.replace('/2', '/2.'))
        table = mathListToPyList(table)
        if len(table) == 0:
            print "!!! Empty data table", mode
            continue

        nFluxTab = []
        grpG = TGraph()
        maxYG = 0
        for i in range(len(table)):
            x, y = table[i]
            maxYG = max(maxYG, y)
            grpG.SetPoint(i, x, y)
            nFlux = computeNFlux(x, s2, m, a, y)
            nFluxTab.append((x, nFlux))

        ## Apply smoothing
        for i in range(1, len(nFluxTab)-1):
            x, y = nFluxTab[i]
            if y != 0: continue
            y = (nFluxTab[i-1][1]+nFluxTab[i+1][1])/2
            nFluxTab[i] = (x, y)

        ## Draw
        maxYF = 0
        grpF = TGraph()
        for i in range(len(nFluxTab)):
            x, nFlux = nFluxTab[i]
            grpF.SetPoint(i, x, nFlux)
            maxYF = max(maxYF, nFlux)
        hFrameG[s2].SetMaximum(max(hFrameG[s2].GetMaximum(), 1.1*maxYG))
        hFrameF[s2].SetMaximum(max(hFrameF[s2].GetMaximum(), 1.1*maxYF))

        sumArea = 0.
        for i in range(1, grpF.GetN()):
            x1, x2 = grpF.GetX()[i-1], grpF.GetX()[i]
            y1, y2 = grpF.GetY()[i-1], grpF.GetY()[i]
            area = (y2+y1)/2*(x2-x1)
            sumArea += area
        hSumNFlux[s2].Fill(a, sumArea)
        hF2Ds[s2].Fill(a, l, sumArea)
        hNDOFs[s2].Fill(a, l, 1./(2*l+1))

        grpG.SetLineColor(colors[len(grpGs[s2])%len(colors)])
        grpG.SetName("grpG_D%d_ss%02d_l%02d_m%02d_a%02d" % (nDim, s2, l, m, a))
        grpG.SetTitle("D=%d s2=%02d l=%02d m=%02d a=%02d" % (nDim, s2, l, m, a))

        grpF.SetLineColor(colors[len(grpFs[s2])%len(colors)])
        grpF.SetName("grpF_D%d_ss%02d_l%02d_m%02d_a%02d" % (nDim, s2, l, m, a))
        grpF.SetTitle("D=%d s2=%02d l=%02d m=%02d a=%02d" % (nDim, s2, l, m, a))

        c.cd(s2+1)
        grpG.Draw("L")
        grpG.SetEditable(False)

        c.cd(s2+4)
        grpF.Draw("L")
        grpF.SetEditable(False)

        grpGs[s2].append(grpG)
        grpFs[s2].append(grpF)

c.Update()
c.Print("cGreybody.png")
c.Print("cGreybody.pdf")

cFlux = TCanvas("cFlux", "cFlux", 300, 300)
opt = ""
for h in reversed(hSumNFlux):
    h.Draw(opt)
    opt = "same"
legFlux = cFlux.BuildLegend()
legFlux.SetFillStyle(0)
legFlux.SetBorderSize(0)
cFlux.Print("cIntegratedFlux.png")
cFlux.Print("cIntegratedFlux.pdf")

cFlux2D = TCanvas("cFlux2D", "cFlux2D", 900, 300)
cFlux2D.Divide(3,1)
for i in range(3):
    cFlux2D.cd(i+1)
    hF2Ds[i].Draw("COLZ")

cNDOF2D = TCanvas("cNDOF2D", "cNDOF2D", 900, 300)
cNDOF2D.Divide(3,1)
for i in range(3):
    cNDOF2D.cd(i+1)
    hNDOFs[i].Draw("COLZ")


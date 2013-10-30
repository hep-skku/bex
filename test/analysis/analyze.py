#!/usr/bin/env python

import sys, os
from xml.dom import minidom
from ROOT import *
from math import *
gROOT.ProcessLine(".x rootlogon.C")

lheFile = sys.argv[1]
outFile = sys.argv[2]

fout = TFile(outFile, "RECREATE")
hN_all = TH1F("hN_all", "Particle multiplicity;Multiplicity;Events", 20, 0, 20)
hN_top = TH1F("hN_top", "Top;Multiplicity;Events", 20, 0, 20)
hN_bot = TH1F("hN_bot", "Bottom;Multiplicity;Events", 20, 0, 20)
hN_hvq = TH1F("hN_hvq", "Top+Bottom;Multiplicity;Events", 20, 0, 20)
hN_hig = TH1F("hN_hig", "Higgs;Multiplicity;Events", 20, 0, 20)
hN_glu = TH1F("hN_glu", "Gluon;Multiplicity;Events", 20, 0, 20)
hN_vec = TH1F("hN_vec", "Vector;Multiplicity;Events", 20, 0, 20)
hN_lep = TH1F("hN_lep", "Lepton;Multiplicity;Events", 20, 0, 20)
hN_oth = TH1F("hN_oth", "Others;Multiplicity;Events", 20, 0, 20)

hPt_top = TH1F("hPt_top", "Top;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_bot = TH1F("hPt_bot", "Bottom;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_hig = TH1F("hPt_hig", "Higgs;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_glu = TH1F("hPt_glu", "Gluon;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_vec = TH1F("hPt_vec", "Vector;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_lep = TH1F("hPt_lep", "Lepton;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)
hPt_oth = TH1F("hPt_oth", "Others;Transverse momentum p_{T} (GeV/c);Entries per 10GeV/c", 50, 0, 500)

lhexml = minidom.parse(lheFile)
print
for i, eventXml in enumerate(lhexml.getElementsByTagName("event")):
    print "\rEvent", i,

    eventContent = eventXml.firstChild.nodeValue.strip().split('\n')
    eventInfoText = eventContent[0].strip().split()
    nParticle, procId = [int(x) for x in eventInfoText[0:2]]
    xsec, qScale, aQCD, aQED = [float(x) for x in eventInfoText[2:]]
    eventInfoText = None

    n_all, n_top, n_bot, n_hig, n_glu, n_lep, n_vec = 0, 0, 0, 0, 0, 0, 0
    for j in range(nParticle):
        line = eventContent[j+1]
        l = line.strip().split()

        id = int(l[0])
        status = int(l[1])
        if status != 1: continue

#        mother1, mother2 = int(l[2]), int(l[3])
#        color1, color2 = int(l[4]), int(l[5])
        px, py, pz = float(l[6]), float(l[7]), float(l[8])
        pt = hypot(px, py)
        e, m = float(l[9]), float(l[10])
#        vt, spin = [float(x) for x in l[11:]]

        n_all += 1
        if abs(id) == 6:
            n_top += 1
            hPt_top.Fill(pt)
        elif abs(id) == 5:
            n_bot += 1
            hPt_bot.Fill(pt)
        elif abs(id) == 21:
            n_glu += 1
            hPt_glu.Fill(pt)
        elif abs(id) in (22, 23, 24):
            n_vec += 1
            hPt_vec.Fill(pt)
        elif abs(id) == 25:
            n_hig += 1
            hPt_hig.Fill(pt)
        elif abs(id) in (11, 13, 15, 12, 14, 16):
            n_lep += 1
            hPt_lep.Fill(pt)
        
    eventContent = None

    hN_all.Fill(n_all)
    hN_top.Fill(n_top)
    hN_bot.Fill(n_bot)
    hN_hvq.Fill(n_top+n_bot)
    hN_glu.Fill(n_glu)
    hN_vec.Fill(n_vec)
    hN_hig.Fill(n_hig)
    hN_lep.Fill(n_lep)
    hN_oth.Fill(n_all-n_top-n_bot-n_glu-n_vec-n_hig-n_lep)

for h in (hN_all, hN_top, hN_bot, hN_hvq, 
          hN_hig, hN_glu, hN_vec, hN_oth, hN_lep,
          hPt_top, hPt_bot, 
          hPt_hig, hPt_glu, hPt_vec, hPt_oth, hPt_lep):
    h.Write()

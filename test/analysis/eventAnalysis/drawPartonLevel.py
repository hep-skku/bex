#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x rootlogon.C")
#gStyle.ForceStyle()

scales = {
  "RS_6D_8TeV":0.003e-3,
  "RS_6D_14TeV":1.260,
  "RS_6D_100TeV":114.542e3,

  "ADD_10D_7TeV":2.432,
  "ADD_10D_8TeV":76.963,
  "ADD_10D_14TeV":140.273e3,
  "ADD_10D_100TeV":760.685e6,
}

cmEnergy = 8
#cmEnergy = 14
#cmEnergy = 100
modNameRS  = "RS_6D_%dTeV" % cmEnergy
modNameADD = "ADD_10D_%dTeV" % cmEnergy
fRS = TFile("hist_%s.root" % modNameRS)
fADD = TFile("hist_%s.root" % modNameADD)
scaleRS  = scales[modNameRS]
scaleADD = scales[modNameADD]

colorsFlavour = {
    "top":kAzure, "bot":kRed, "glu":kGreen+1,
    "hig":kGray, "vec":kOrange,
    "lep":kMagenta+1, "oth":kBlack, 
}
labelsFlavour = {
    "top":"Top", "bot":"Bottom", "glu":"Gluon",
    "hig":"Higgs", "vec":"W/Z/#gamma",
    "oth":"Others", "lep":"Leptons",
}

# Multiplicity plots
hN_partons = []
hStackN_ADD = THStack("hStackN_ADD", "ADD")
hStackN_RS  = THStack("hStackN_RS", "RS")
hStackPt_ADD = THStack("hStackPt_ADD", "ADD")
hStackPt_RS = THStack("hStackPt_RS", "RS")
legN_ADD = TLegend(0.65, 0.65, 0.92, 0.92)
legN_RS  = TLegend(0.65, 0.65, 0.92, 0.92)
for leg in (legN_ADD, legN_RS):
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
for flavour in ("top", "bot", "hig", "glu", "vec", "lep", "oth"):
    hN_ADD = fADD.Get("hN_%s" % flavour)
    hN_RS  = fRS.Get("hN_%s" % flavour)

    hN_ADD.Scale(scaleADD/1000)
    hN_RS.Scale(scaleRS/1000)
    hN_ADD.SetLineColor(colorsFlavour[flavour])
    hN_RS.SetLineColor(colorsFlavour[flavour])
    hN_ADD.SetLineWidth(2)
    hN_RS.SetLineWidth(2)

    hStackN_ADD.Add(hN_ADD)
    hStackN_RS.Add(hN_RS)
    legN_ADD.AddEntry(hN_ADD, labelsFlavour[flavour], "l")
    legN_RS.AddEntry(hN_RS, labelsFlavour[flavour], "l")

    hPt_ADD = fADD.Get("hPt_%s" % flavour)
    hPt_RS  = fRS.Get("hPt_%s" % flavour)

    hPt_ADD.Scale(scaleADD/1000)
    hPt_RS.Scale(scaleRS/1000)
    hPt_ADD.SetLineColor(colorsFlavour[flavour])
    hPt_RS.SetLineColor(colorsFlavour[flavour])
    hPt_ADD.SetLineWidth(2)
    hPt_RS.SetLineWidth(2)
    hPt_ADD.SetMaximum(scaleADD*10)
    hPt_ADD.SetMinimum(scaleADD*1e-3)
    hPt_RS.SetMaximum(scaleRS*10)
    hPt_RS.SetMinimum(scaleRS*1e-3)

    hStackPt_ADD.Add(hPt_ADD)
    hStackPt_RS.Add(hPt_RS)

    hN_partons.extend([hN_ADD, hN_RS, hPt_ADD, hPt_RS])

cN_ADD = TCanvas("cN_ADD", "cN_ADD", 500, 500)
hStackN_ADD.Draw("nostack")
hStackN_ADD.GetYaxis().SetTitle("Events per 1fb^{-1}")
hStackN_ADD.GetXaxis().SetTitle("Particle multiplicity")
hStackN_ADD.GetXaxis().SetRangeUser(0,10)
hStackN_ADD.Draw("nostack")
legN_ADD.Draw()
cN_ADD.Print("cN_%s.png" % modNameADD)

cN_RS = TCanvas("cN_RS", "cN_RS", 500, 500)
hStackN_RS.Draw("nostack")
hStackN_RS.GetYaxis().SetTitle("Events per 1fb^{-1}")
hStackN_RS.GetXaxis().SetTitle("Particle multiplicity")
hStackN_RS.GetXaxis().SetRangeUser(0,10)
hStackN_RS.Draw("nostack")
legN_RS.Draw()
cN_RS.Print("cN_%s.png" % modNameRS)

cPt_ADD = TCanvas("cPt_ADD", "cPt_ADD", 500, 500)
cPt_ADD.SetLogy()
hStackPt_ADD.Draw("nostack")
hStackPt_ADD.GetYaxis().SetTitle("Entries/20GeV/c per 1fb^{-1}")
hStackPt_ADD.GetXaxis().SetTitle("Transverse momentum p_{T} (GeV/c)")
hStackPt_ADD.Draw("nostack")
legN_ADD.Draw()
cPt_ADD.Print("cPt_%s.png" % modNameADD)

cPt_RS = TCanvas("cPt_RS", "cPt_RS", 500, 500)
cPt_RS.SetLogy()
hStackPt_RS.Draw("nostack")
hStackPt_RS.GetYaxis().SetTitle("Entries/20GeV/c per 1fb^{-1}")
hStackPt_RS.GetXaxis().SetTitle("Transverse momentum p_{T} (GeV/c)")
hStackPt_RS.Draw("nostack")
legN_RS.Draw()
cPt_RS.Print("cPt_%s.png" % modNameRS)




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
fRS = TFile("results_%s.root" % modNameRS)
fADD = TFile("results_%s.root" % modNameADD)
scaleRS  = scales[modNameRS]
scaleADD = scales[modNameADD]

hST_RS  = fRS.Get("hST")
hST_ADD = fADD.Get("hST")

hNLeptons_RS  = fRS.Get("hNLeptons")
hNLeptons_ADD = fADD.Get("hNLeptons")
hNJets_RS  = fRS.Get("hNJets")
hNJets_ADD = fADD.Get("hNJets")
hNBjets_RS  = fRS.Get("hNBjets")
hNBjets_ADD = fADD.Get("hNBjets")

hBjetPt_RS  = fRS.Get("hBjetPt")
hBjetPt_ADD = fADD.Get("hBjetPt")
hJetPt_RS  = fRS.Get("hJetPt")
hJetPt_ADD = fADD.Get("hJetPt")

for h in (hST_RS, hNLeptons_RS, hNJets_RS, hNBjets_RS, hJetPt_RS, hBjetPt_RS):
    h.Scale(scaleRS/1000)
    h.SetLineWidth(2)
for h in (hST_ADD, hNLeptons_ADD, hNJets_ADD, hNBjets_ADD, hJetPt_ADD, hBjetPt_ADD):
    h.Scale(scaleADD/1000)
    h.SetLineWidth(2)

cST_ADD = TCanvas("cST_ADD", "cST_ADD", 500, 500)
hST_ADD.GetXaxis().SetNdivisions(505)
hST_ADD.GetYaxis().SetTitle("Events/500GeV/c per 1fb^{-1}")
hST_ADD.GetXaxis().SetTitle("Scalar sum S_{T} (GeV/c)")
hST_ADD.Draw()
cST_ADD.Print("hST_%s.png" % modNameADD)
cST_RS = TCanvas("cST_RS", "cST_RS", 500, 500)
hST_RS.GetXaxis().SetNdivisions(505)
hST_RS.GetYaxis().SetTitle("Events/500GeV/c per 1fb^{-1}")
hST_RS.GetXaxis().SetTitle("Scalar sum S_{T} (GeV/c)")
hST_RS.Draw()
cST_RS.Print("hST_%s.png" % modNameRS)

## Particle multiplicity
cN_ADD = TCanvas("cN_ADD", "cN_ADD", 500, 500)
legN_ADD = TLegend(0.70, 0.75, 0.93, 0.93)
legN_ADD.SetFillStyle(0)
legN_ADD.SetBorderSize(0)
hNLeptons_ADD.SetLineColor(kOrange)
hNJets_ADD.SetLineColor(kBlack)
hNBjets_ADD.SetLineColor(kRed)
hStackN_ADD = THStack("hStackN_ADD", "hStackN_ADD;Multiplicity;Events per 1fb^{-1}")
hStackN_ADD.Add(hNLeptons_ADD)
hStackN_ADD.Add(hNJets_ADD)
hStackN_ADD.Add(hNBjets_ADD)
legN_ADD.AddEntry(hNLeptons_ADD, "Leptons", "l")
legN_ADD.AddEntry(hNJets_ADD, "Jets", "l")
legN_ADD.AddEntry(hNBjets_ADD, "B jets", "l")
hStackN_ADD.Draw("nostack")
legN_ADD.Draw()
cN_ADD.Print("cDetN_%s.png" % modNameADD)

cN_RS = TCanvas("cN_RS", "cN_RS", 500, 500)
legN_RS = TLegend(0.70, 0.75, 0.93, 0.93)
legN_RS.SetFillStyle(0)
legN_RS.SetBorderSize(0)
hNLeptons_RS.SetLineColor(kOrange)
hNJets_RS.SetLineColor(kBlack)
hNBjets_RS.SetLineColor(kRed)
hStackN_RS = THStack("hStackN_RS", "hStackN_RS;Multiplicity;Events per 1fb^{-1}")
hStackN_RS.Add(hNLeptons_RS)
hStackN_RS.Add(hNJets_RS)
hStackN_RS.Add(hNBjets_RS)
legN_RS.AddEntry(hNLeptons_RS, "Leptons", "l")
legN_RS.AddEntry(hNJets_RS, "Jets", "l")
legN_RS.AddEntry(hNBjets_RS, "B jets", "l")
hStackN_RS.Draw("nostack")
legN_RS.Draw()
cN_RS.Print("cDetN_%s.png" % modNameRS)

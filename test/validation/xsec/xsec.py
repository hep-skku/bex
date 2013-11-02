#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x ../rootlogon.C")

xsec_ADD6_cmEnergy = [
  (  7,   0.650*1e+0, 0.013*1e+0),
  (  8,  20.570*1e+0, 0.333*1e+0),
  (  9,   0.195*1e+3, 0.003*1e+3),
  ( 10,   0.978*1e+3, 0.012*1e+3),
  ( 11,   3.337*1e+3, 0.039*1e+3),
  ( 12,   8.847*1e+3, 0.095*1e+3),
  ( 13,  19.644*1e+3, 0.197*1e+3),
  ( 14,  38.223*1e+3, 0.362*1e+3),
#  ( 20, 481.572*1e+3, 3.589*1e+3),
#  ( 50,  26.168*1e+6, 0.128*1e+6),
#  (100, 208.349*1e+6, 0.824*1e+6),
]

xsec_ADD11_cmEnergy = [
  (  7,   2.950*1e+0, 0.072*1e+0),
  (  8,  92.482*1e+0, 1.874*1e+0),
  (  9,   0.875*1e+3, 0.015*1e+3),
  ( 10,   4.396*1e+3, 0.069*1e+3),
  ( 11,  15.024*1e+3, 0.216*1e+3),
  ( 12,  39.734*1e+3, 0.531*1e+3),
  ( 13,  88.045*1e+3, 1.102*1e+3),
  ( 14, 170.563*1e+3, 2.018*1e+3),
#  ( 20,   2.148*1e+6, 0.020*1e+6),
#  ( 50, 115.542*1e+6, 0.703*1e+6),
#  (100, 910.230*1e+6, 4.537*1e+6),
]

xsec_RS_cmEnergy = [
  (  7,  0.003*1e-3, 0.000*1e-3),
  (  8,  0.145*1e-3, 0.003*1e-3),
  (  9,  0.002*1e+0, 0.000*1e+0),
  ( 10,  0.014*1e+0, 0.000*1e+0),
  ( 11,  0.062*1e+0, 0.001*1e+0),
  ( 12,  0.204*1e+0, 0.002*1e+0),
  ( 13,  0.546*1e+0, 0.006*1e+0),
  ( 14,  1.260*1e+0, 0.013*1e+0),
#  ( 20, 33.697*1e+0, 0.294*1e+0),
#  ( 50,  7.386*1e+3, 0.045*1e+3),
#  (100,114.542*1e+3, 0.578*1e+3),
]

grpXsec_ADD6_cmEnergy = TGraphErrors()
grpXsec_ADD11_cmEnergy = TGraphErrors()
grpXsec_RS_cmEnergy = TGraphErrors()

for i in range(len(xsec_RS_cmEnergy)):
    beamEnergyADD6, xsecADD6, xsecErrADD6 = xsec_ADD6_cmEnergy[i]
    beamEnergyADD11, xsecADD11, xsecErrADD11 = xsec_ADD11_cmEnergy[i]
    beamEnergyRS, xsecRS, xsecErrRS = xsec_RS_cmEnergy[i]

    grpXsec_ADD6_cmEnergy.SetPoint(i, beamEnergyADD6, xsecADD6)
    grpXsec_ADD6_cmEnergy.SetPointError(i, 0, xsecErrADD6)

    grpXsec_ADD11_cmEnergy.SetPoint(i, beamEnergyADD11, xsecADD11)
    grpXsec_ADD11_cmEnergy.SetPointError(i, 0, xsecErrADD11)

    grpXsec_RS_cmEnergy.SetPoint(i, beamEnergyRS, xsecRS)
    grpXsec_RS_cmEnergy.SetPointError(i, 0, xsecErrRS)

cXsec_cmEnergy = TCanvas("cXsec_cmEnergy", "cXsec_cmEnergy", 500, 500)
cXsec_cmEnergy.SetLogy()
legXsec_cmEnergy = TLegend(.65, .20, .90, .40)
legXsec_cmEnergy.SetBorderSize(0)
legXsec_cmEnergy.SetFillStyle(0)
hFrameXsec_cmEnergy = TH1F("hFrameXsec_cmEnergy", "hFrame;Center of mass energy #sqrt{s} (TeV);Cross section (fb)", 100, 7, 14)
hFrameXsec_cmEnergy.SetMinimum(1e-6)
hFrameXsec_cmEnergy.SetMaximum(1e+7)
hFrameXsec_cmEnergy.GetYaxis().SetNdivisions(505)
hFrameXsec_cmEnergy.GetYaxis().SetLabelSize(0.04)
hFrameXsec_cmEnergy.Draw()
grpXsec_ADD11_cmEnergy.SetMarkerColor(kRed)
grpXsec_ADD6_cmEnergy.SetMarkerColor(kMagenta)
grpXsec_RS_cmEnergy.SetMarkerColor(kBlue)
grpXsec_ADD11_cmEnergy.SetLineColor(kMagenta)
grpXsec_ADD6_cmEnergy.SetLineColor(kRed)
grpXsec_RS_cmEnergy.SetLineColor(kBlue)
legXsec_cmEnergy.AddEntry(grpXsec_ADD11_cmEnergy, "ADD(11D)", "l")
legXsec_cmEnergy.AddEntry(grpXsec_ADD6_cmEnergy, "ADD(6D)", "l")
legXsec_cmEnergy.AddEntry(grpXsec_RS_cmEnergy, "RS", "l")
grpXsec_ADD11_cmEnergy.Draw("LP")
grpXsec_ADD6_cmEnergy.Draw("LP")
grpXsec_RS_cmEnergy.Draw("LP")
legXsec_cmEnergy.Draw()

cXsec_cmEnergy.Print("xsec_cmEnergy.png")
cXsec_cmEnergy.Print("xsec_cmEnergy.pdf")

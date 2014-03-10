#!/usr/bin/env python

from ROOT import *
gROOT.ProcessLine(".x ../rootlogon.C")

f = TFile("debug.root")
colors = [kRed, kBlue, kMagenta]

legend = TLegend(0.7, 0.2, 0.95, 0.35)
legend.SetBorderSize(0)
legend.SetFillStyle(0)

## Draw flux curves
grpFlux = []
cFlux = TCanvas("cFlux", "cFlux", 500, 500)
opt = "AP"
maxY = 0
for i in range(3):
    grp = f.Get("grpFlux_%d" % i)
    grp.SetTitle("S2=%d;Particle energy #omega;Flux" % i)
    grpFlux.append(grp)
    grp.SetMarkerStyle(0)
    grp.SetMarkerColor(colors[i])
    grp.SetLineColor(colors[i])
    grp.Draw(opt)
    opt = "P"
    legend.AddEntry(grp, "S2=%d" % i, "l")
    maxY = max(TMath.MaxElement(grp.GetN(), grp.GetY()), maxY)
for grp in grpFlux:
    grp.SetMaximum(1.2*maxY)
legend.Draw()

## Draw Wien's displacement law
grpTemVsPeakPos = []
cTemVsPeakPos = TCanvas("cTemVsPeakPos", "cTemVsPeakPos", 500, 500)
opt = "AP"
maxY = 0
for i in range(3):
    grp = f.Get("grpTemVsPeakPos_%d" % i)
    grp.SetTitle("S2=%d;BH temperature;Peak position" % i)
    grpTemVsPeakPos.append(grp)
    grp.SetMarkerStyle(0)
    grp.SetMarkerColor(colors[i])
    grp.SetLineColor(colors[i])
    grp.Draw(opt)
    opt = "P"
    maxY = max(TMath.MaxElement(grp.GetN(), grp.GetY()), maxY)
for grp in grpTemVsPeakPos:
    grp.SetMaximum(1.2*maxY)
legend.Draw()

## Draw Stephen-Boltzman's law
grpTemVsTotalFlux = []
cTemVsTotalFlux = TCanvas("cTemVsTotalFlux", "cTemVsTotalFux", 500, 500)
cTemVsTotalFlux.SetLogy()
cTemVsTotalFlux.SetLogx()
opt = "AP"
maxY = 0
for i in range(3):
    grp = f.Get("grpTemVsTotalFlux_%d" % i)
    grp.SetTitle("S2=%d;BH temperature;Total flux" % i)
    grpTemVsTotalFlux.append(grp)
    grp.SetMarkerStyle(0)
    grp.SetMarkerColor(colors[i])
    grp.SetLineColor(colors[i])
    grp.Draw(opt)
    opt = "P"
    maxY = max(TMath.MaxElement(grp.GetN(), grp.GetY()), maxY)
for grp in grpTemVsTotalFlux:
    grp.SetMaximum(1.2*maxY)
legend.Draw()

## Draw ordinary histograms
cPlots = []
hPlots = []
for plotName in ["NDecay", "EDecay", "MJLoss"]:
    c = TCanvas("c%s" % plotName, plotName, 500, 500)
    cPlots.append(c)

    h = f.Get("h%s" % plotName)
    if h.IsA().InheritsFrom("TH2"):
        h.Draw("COLZ")
    else:
        h.Draw()
    hPlots.append(h)

## Draw BH history plot
historyPlots = []
cHistory = TCanvas("cHistory", "BH history", 500, 500)
hFrameHistory = TH1F("hFrameHistory", "frame;Iteration;BH mass / initial mass", 100, 0, 40)
#hFrameHistory.SetMaximum(15000)
#hFrameHistory.SetMinimum(4500)
hFrameHistory.Draw()
for plotName in [x.GetName() for x in f.GetDirectory("MBHHistory").GetListOfKeys()]:
    print plotName
    grp = f.Get("MBHHistory/%s" % plotName)
    if grp == None: continue
    if grp.IsA().GetName() != "TGraph": continue
    grp.SetMarkerStyle(0)
    grp.Draw("LP")
    historyPlots.append(grp)

## Print them all
for c in [cFlux, cTemVsPeakPos, cTemVsTotalFlux, cHistory] + cPlots:
    c.Print("%s.png" % c.GetName())
    c.Print("%s.pdf" % c.GetName())

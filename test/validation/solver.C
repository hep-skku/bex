#include <iostream>
#include <cmath>
#include <TSystem.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TLegend.h>

using namespace std;

int nDim_;
const int minDim = 4;
const double maxA = 2;
const double rs = 1;

double computeRh(const double a)
{
  double x = rs;

/*
  if ( nDim_ == 4 )
  {
    const double det = rs*rs - 4*a*a;
    if ( det < 0 ) return -1;

    return (rs + sqrt(det))/2;
  }
  else if ( nDim_ == 5 )
  {
    const double det = rs*rs - a*a;
    if ( det < 0 ) return -1;

    return sqrt(det);
  }
  else
*/
  {  
    double dx = 1e9;
    const double p = 1./(nDim_-3);
    while ( std::abs(dx) > 1e-3*std::abs(x) )
    {
      const double det = x*x+a*a;
      const double f      = x*pow(det, p) - rs*pow(x, 2*p);
      const double fprime = pow(det, p) + 2*p*x*x*pow(det, p-1) - 2*p*rs*pow(x, 2*p-1);
      dx = f/fprime;
      x -= dx;
    }

    // Horizon radius must be positive definite
    if ( x < 0 ) x = 0;

    return x;
  }

  return -1;
}

void solver()
{
  gROOT->ProcessLine(".x rootlogon.C");

  TCanvas* c = new TCanvas("c", "c", 500, 500);
  TH1F* hFrame = new TH1F("hFrame", "hFrame;Rotation parameter a;r_{h}/r_{s}", 100, 0, maxA);
  hFrame->SetMinimum(0);
  hFrame->SetMaximum(rs*1.1);
  hFrame->Draw();

  TLegend* leg = new TLegend(0.22, 0.22, 0.45, 0.5);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  TGraph* graphs[10];
  int colors[10] = {kBlack, kRed, kOrange+1, kPink+1, kMagenta+1, kGreen+1, kAzure+1, kBlue};
  for ( nDim_ = minDim; nDim_ <= 11; ++nDim_ )
  {
    TGraph* grp = new TGraph();
    for ( int i=0; i<100; ++i )
    {
      const double a = maxA*i/100;
      const double y = computeRh(a);
      grp->SetPoint(i, a, y);
    }
    grp->SetLineColor(colors[nDim_-minDim]);
    grp->Draw("L");
    graphs[nDim_-minDim] = grp;
    leg->AddEntry(grp, Form("D%d", nDim_), "l");
  }
  leg->Draw();
}

#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

const double kL = 11.2;
const double L = 1;
const double pi = TMath::Pi();

double profile_s(double* x, double* p)
{
  const double y = x[0];
  if ( pow((y-L)/L, 2) < 1e-3 ) return 30;
  else return 0;
  return 0;
}

double profile_f(double* x, double* p)
{
  const double y = x[0];
  const double sign = p[0];
  const double nu = p[1];

  const double det = 1+2*nu*sign;
  const double norm = sqrt((kL/pi*det)/(exp(kL*det)-1));

  return 1./sqrt(L)*norm*exp( 0.5*det*kL/L*y );

}

double profile_v(double* x, double* p)
{
  //const double y = x[0];
  return 1./sqrt(L);
}

void particleProfile()
{
  const double minL = 0, maxL = L;

  TF1* fS = new TF1("fS", profile_s, minL, maxL, 0);
  TF1* fV = new TF1("fV", profile_v, minL, maxL, 0);
  TF1* fFQ1 = new TF1("fFQ1", profile_f, minL, maxL, 2);
  TF1* fFQ2 = new TF1("fFQ2", profile_f, minL, maxL, 2);
  TF1* fFQ3 = new TF1("fFQ3", profile_f, minL, maxL, 2);
  TF1* fFu1 = new TF1("fFu1", profile_f, minL, maxL, 2);
  TF1* fFu2 = new TF1("fFu2", profile_f, minL, maxL, 2);
  TF1* fFu3 = new TF1("fFu3", profile_f, minL, maxL, 2);
  TF1* fFd1 = new TF1("fFd1", profile_f, minL, maxL, 2);
  TF1* fFd2 = new TF1("fFd2", profile_f, minL, maxL, 2);
  TF1* fFd3 = new TF1("fFd3", profile_f, minL, maxL, 2);
  fFQ1->SetParameters(+1, -0.579);
  fFQ2->SetParameters(+1, -0.517);
  fFQ3->SetParameters(+1, -0.473);
  fFu1->SetParameters(-1, +0.742);
  fFu2->SetParameters(-1, +0.558);
  fFu3->SetParameters(-1, -0.339);
  fFd1->SetParameters(-1, +0.711);
  fFd2->SetParameters(-1, +0.666);
  fFd3->SetParameters(-1, +0.533);

  TCanvas* c = new TCanvas("cProfile", "cProfile", 500, 500);
  c->SetLogy();

  fS->SetLineColor(kBlack);
  fV->SetLineColor(kOrange+1);

  fFQ1->SetLineColor(kGreen+1);
  fFu1->SetLineColor(kAzure+1);
  fFd1->SetLineColor(kRed+1);

  fFQ2->SetLineColor(kGreen+2);
  fFu2->SetLineColor(kAzure+2);
  fFd2->SetLineColor(kRed+2);

  fFQ3->SetLineColor(kGreen+3);
  fFu3->SetLineColor(kAzure+0);
  fFd3->SetLineColor(kRed+0);

  fFQ3->SetLineWidth(3);
  fFu3->SetLineWidth(3);
  fFd3->SetLineWidth(3);

  TLegend* leg = new TLegend(0.17, 0.75, 0.60, 0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->AddEntry(fFQ1, "Quark Q_{1}", "l");
  leg->AddEntry(fFQ2, "Quark Q_{2}", "l");
  leg->AddEntry(fFQ3, "Quark Q_{3}", "l");

  leg->AddEntry(fFu1, "Quark u_{1}", "l");
  leg->AddEntry(fFu2, "Quark u_{2}", "l");
  leg->AddEntry(fFu3, "Quark u_{3}", "l");

  leg->AddEntry(fFd1, "Quark d_{1}", "l");
  leg->AddEntry(fFd2, "Quark d_{2}", "l");
  leg->AddEntry(fFd3, "Quark d_{3}", "l");

  leg->AddEntry(fV, "Boson", "l");


  TH1F* hFrame = new TH1F("hFrame", ";y position;Profile amplitude #chi(y)", 100, minL, maxL);
  hFrame->SetMinimum(1e-2);
  hFrame->SetMaximum(10);
  hFrame->Draw();

  fV->Draw("same");
  fFQ1->Draw("same");
  fFu1->Draw("same");
  fFd1->Draw("same");

  fFQ2->Draw("same");
  fFu2->Draw("same");
  fFd2->Draw("same");

  fFQ3->Draw("same");
  fFu3->Draw("same");
  fFd3->Draw("same");

  leg->Draw();
}

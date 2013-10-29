#include <iostream>
#include <TF1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <THStack.h>

const double pi = TMath::Pi();
const double kL = 11.3*pi;
const double L = 1;
const double mD = 1;

double profile_sca(double* x, double* p)
{
  const double y = x[0];
  if ( pow((y-L)/L, 2) < 1e-3 ) return 10;
  return 0;
}

double profile_fer(double* x, double* p)
{
  const double y = x[0];
  const double sign = p[0];
  const double nu = p[1];

  const double det = 1+2*nu*sign;
  const double norm = sqrt((kL/pi*det)/(exp(kL*det)-1));

  return 1./sqrt(L)*norm*exp( 0.5*det*kL/L*y );
}

double profile_vec(double* x, double* p)
{
  return 1./sqrt(L);
}

double profileSqr_sca(double* x, double* p)
{
  const double y = x[0];
  if ( pow((y-L)/L, 2) < 1e-3 ) return 100;
  return 0;
}

double profileSqr_fer(double* x, double* p)
{
  const double y = x[0];
  const double sign = p[0];
  const double nu = p[1];

  const double det = 1+2*nu*sign;
  const double norm2 = (kL/pi*det)/(exp(kL*det)-1);

  return 1./L*norm2*exp( det*kL/L*y );
}

double profileSqr_vec(double* x, double* p)
{
  return 1/L;
}

double dfactorFtn_fer(double* x, double* p)
{
  const double y = x[0];
  return exp(3*kL/L*(y-L))*profile_fer(x, p);
}

double dfactorFtn_vec(double* x, double* p)
{
  const double y = x[0];
  return exp(3*kL/L*(y-L))*profile_vec(x, p);
}

double cfactorFtn_ff(double* x, double* p)
{
  const double y = x[0];
  const double sign1 = p[0];
  const double nu1 = p[1];
  const double sign2 = p[2];
  const double nu2 = p[3];

  return exp(2*kL/L*(y-L))*profile_fer(x, &p[0])*profile_fer(x, &p[2]);
}

double cfactorFtn_gg(double* x, double* p)
{
  const double y = x[0];
  return exp(2*kL/L*(y-L))*profile_vec(x, p)*profile_vec(x, p);
}

double cfactorFtn_gf(double* x, double* p)
{
  const double y = x[0];
  return exp(2*kL/L*(y-L))*profile_vec(x, p)*profile_fer(x, p);
}

void particleProfile()
{
  const double minL = 0, maxL = L;

  const unsigned int nFermion = 9;
  const double fermionPar1[nFermion] = {1, 1, 1, -1, -1, -1, -1, -1, -1};
  const double fermionPar2[nFermion] = {-0.579, -0.517, -0.473, +0.742, +0.558, -0.339, +0.711, +0.666, +0.533};
  const char* fermionNames[nFermion] = {"Q1", "Q2", "Q3", "u1", "u2", "u3", "d1", "d2", "d3"};

  TF1* fSca = new TF1("fSca", profile_sca, minL, maxL, 0);
  TF1* fVec = new TF1("fVec", profile_vec, minL, maxL, 0);
  TF1* fFs[nFermion], * fSqrFs[nFermion], * fDfacFs[nFermion];
  for ( unsigned int i=0; i<nFermion; ++i )
  {
    fFs[i] = new TF1(Form("fF%s", fermionNames[i]), profile_fer, minL, maxL, 2);
    fSqrFs[i] = new TF1(Form("fSqrF%s", fermionNames[i]), profileSqr_fer, minL, maxL, 2);
    fDfacFs[i] = new TF1(Form("fDfacF%s", fermionNames[i]), dfactorFtn_fer, minL, maxL, 2);

    fFs[i]->SetParameters(fermionPar1[i], fermionPar2[i]);
    fSqrFs[i]->SetParameters(fermionPar1[i], fermionPar2[i]);
    fDfacFs[i]->SetParameters(fermionPar1[i], fermionPar2[i]);
  }

  TF1* fSqrSca = new TF1("fSqrSca", profileSqr_sca, minL, maxL, 0);
  TF1* fSqrVec = new TF1("fSqrVec", profileSqr_vec, minL, maxL, 0);
  TF1* fDfacVec = new TF1("fDfacVec", dfactorFtn_vec, minL, maxL, 0);

  double cbb = 0, cgb = 0, cll = 0, cgl = 0, cbl = 0;
  TF1* fCfactorFF = new TF1("fCfactorFF", cfactorFtn_ff, minL, maxL, 4);
  TF1* fCfactorGF = new TF1("fCfactorGF", cfactorFtn_gf, minL, maxL, 2);
  TF1* fCfactorGG = new TF1("fCfactorGG", cfactorFtn_gg, minL, maxL, 0);
  double cgg = fCfactorGG->Integral(minL, maxL);
  for ( unsigned int i=0; i<nFermion; ++i )
  {
    if ( i == 5 ) continue;
    fCfactorGF->SetParameters(fermionPar1[i], fermionPar2[i]);
    if ( i == 2 or i == 8 ) cgb += fCfactorGF->Integral(minL, maxL);
    else cgl += fCfactorGF->Integral(minL, maxL);
    for ( unsigned int j=0; j<nFermion; ++j )
    {
      if ( j == 5 ) continue;
      fCfactorFF->SetParameters(fermionPar1[i], fermionPar2[i], fermionPar1[j], fermionPar2[j]);
      const double cxx = fCfactorFF->Integral(minL, maxL);

      if ( (i == 2 or i == 8) and (j == 2 or j == 8) ) cbb += cxx;
      else if ( (i == 2 or i == 8) or (j == 2 or j == 8) ) cbl += cxx;
      else cll += cxx;
      cout << Form("%s%s ", fermionNames[i], fermionNames[j]) << cxx << endl;
    }
  }
  cout << "Cgg = " << cgg << endl;
  cout << "Cbb = " << cbb << endl;
  cout << "Cgb = " << cgb << endl;
  cout << "Cll = " << cll << endl;
  cout << "Cbl = " << cbl << endl;
  cout << "Cgl = " << cgl << endl;

  return;

/*
  TH1F* hProfIR = new TH1F("hProfIR", "Profile at IR;;Amplitude at IR brane", 10, 0, 10);
  hProfIR->GetXaxis()->SetBinLabel( 1, "Vector"); hProfIR->SetBinContent( 1, fVec->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 2, "Q_{1}"); hProfIR->SetBinContent( 2, fFQ1->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 3, "Q_{2}"); hProfIR->SetBinContent( 3, fFQ2->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 4, "Q_{3}"); hProfIR->SetBinContent( 4, fFQ3->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 5, "u_{1}"); hProfIR->SetBinContent( 5, fFu1->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 6, "u_{2}"); hProfIR->SetBinContent( 6, fFu2->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 7, "u_{3}"); hProfIR->SetBinContent( 7, fFu3->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 8, "d_{1}"); hProfIR->SetBinContent( 8, fFd1->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel( 9, "d_{2}"); hProfIR->SetBinContent( 9, fFd2->Eval(L));
  hProfIR->GetXaxis()->SetBinLabel(10, "d_{3}"); hProfIR->SetBinContent(10, fFd3->Eval(L));

  const double minDL = 0.8*L;
  TH1F* hSqrDSca = new TH1F("hSqrDSca", "Scator;BH radius;Relative D factor", 100, minDL, L);
  TH1F* hSqrDVec = new TH1F("hSqrDVec", "Vector;BH radius;Relative D factor", 100, minDL, L);
  TH1F* hSqrDFQ3 = new TH1F("hSqrDFQ3", "Q3;BH radius;Relative D factor", 100, minDL, L);
  TH1F* hSqrDFu3 = new TH1F("hSqrDFu3", "u3;BH radius;Relative D factor", 100, minDL, L);
  TH1F* hSqrDFd3 = new TH1F("hSqrDFd3", "d3;BH radius;Relative D factor", 100, minDL, L);
  TH1F* hSqrDFxy = new TH1F("hSqrDFxy", "others;BH radius;Relative D factor", 100, minDL, L);
  hSqrDSca->SetFillColor(kWhite);
  hSqrDVec->SetFillColor(kOrange+1);
  hSqrDFQ3->SetFillColor(kGreen+3);
  hSqrDFu3->SetFillColor(kAzure+0);
  hSqrDFd3->SetFillColor(kRed+0);
  hSqrDFxy->SetFillColor(kYellow);

  TH1F* hRelDSca = new TH1F("hRelDSca", "Scator;BH radius;Relative D factor", 1, 0, L);
  TH1F* hRelDVec = new TH1F("hRelDVec", "Vector;BH radius;Relative D factor", 1, 0, L);
  TH1F* hRelDFQ3 = new TH1F("hRelDFQ3", "Q3;BH radius;Relative D factor", 1, 0, L);
  TH1F* hRelDFu3 = new TH1F("hRelDFu3", "u3;BH radius;Relative D factor", 1, 0, L);
  TH1F* hRelDFd3 = new TH1F("hRelDFd3", "d3;BH radius;Relative D factor", 1, 0, L);
  TH1F* hRelDFxy = new TH1F("hRelDFxy", "others;BH radius;Relative D factor", 1, 0, L);
  hRelDSca->SetFillColor(kWhite);
  hRelDVec->SetFillColor(kOrange+1);
  hRelDFQ3->SetFillColor(kGreen+3);
  hRelDFu3->SetFillColor(kAzure+0);
  hRelDFd3->SetFillColor(kRed+0);
  hRelDFxy->SetFillColor(kYellow);
   
  for ( int i=1; i<=100; ++i )
  {
    const double x = hSqrDVec->GetXaxis()->GetBinCenter(i);
    const double dSca = 1;

    const double dSqrVec = fSqrVec->Integral(x, L);
    const double dSqrFQ3 = fSqrFQ3->Integral(x, L);
    const double dSqrFu3 = fSqrFu3->Integral(x, L);
    const double dSqrFd3 = fSqrFd3->Integral(x, L);
    double dSqrFxy = 0;
    dSqrFxy += fSqrFQ1->Integral(x, L);
    dSqrFxy += fSqrFu1->Integral(x, L);
    dSqrFxy += fSqrFd1->Integral(x, L);
    dSqrFxy += fSqrFQ2->Integral(x, L);
    dSqrFxy += fSqrFu2->Integral(x, L);
    dSqrFxy += fSqrFd2->Integral(x, L);
    const double dSqrSum = dSca+dSqrVec+dSqrFQ3+dSqrFu3+dSqrFd3+dSqrFxy;
    hSqrDSca->SetBinContent(i, dSca/dSqrSum);
    hSqrDVec->SetBinContent(i, dSqrVec/dSqrSum);
    hSqrDFQ3->SetBinContent(i, dSqrFQ3/dSqrSum);
    hSqrDFu3->SetBinContent(i, dSqrFu3/dSqrSum);
    hSqrDFd3->SetBinContent(i, dSqrFd3/dSqrSum);
    hSqrDFxy->SetBinContent(i, dSqrFxy/dSqrSum);
  }

  for ( int i=1; i<=hRelDVec->GetNbinsX(); ++i )
  {
    const double x = hSqrDVec->GetXaxis()->GetBinLowEdge(i);
    const double dSca = 1;

    const double dVec = pow(fDfacVec->Integral(x, L), 2);
    const double dFQ3 = pow(fDfacFQ3->Integral(x, L), 2);
    const double dFu3 = pow(fDfacFu3->Integral(x, L), 2);
    const double dFd3 = pow(fDfacFd3->Integral(x, L), 2);
    double dFxy = 0;
    dFxy += pow(fDfacFQ1->Integral(x, L), 2);
    dFxy += pow(fDfacFu1->Integral(x, L), 2);
    dFxy += pow(fDfacFd1->Integral(x, L), 2);
    dFxy += pow(fDfacFQ2->Integral(x, L), 2);
    dFxy += pow(fDfacFu2->Integral(x, L), 2);
    dFxy += pow(fDfacFd2->Integral(x, L), 2);
    const double dSum = dSca+dVec+dFQ3+dFu3+dFd3+dFxy;
    hRelDSca->SetBinContent(i, dSca/dSum);
    hRelDVec->SetBinContent(i, dVec/dSum);
    hRelDFQ3->SetBinContent(i, dFQ3/dSum);
    hRelDFu3->SetBinContent(i, dFu3/dSum);
    hRelDFd3->SetBinContent(i, dFd3/dSum);
    hRelDFxy->SetBinContent(i, dFxy/dSum);

  }

  THStack* hSqrD = new THStack();
  hSqrD->SetTitle("Relative D factor;Overlap size;D factor fraction");
  hSqrD->Add(hSqrDFxy);
  hSqrD->Add(hSqrDFd3);
  hSqrD->Add(hSqrDFu3);
  hSqrD->Add(hSqrDFQ3);
  hSqrD->Add(hSqrDVec);
  hSqrD->Add(hSqrDSca);
  TLegend* legSqrD = new TLegend(0.75, 0.65, 0.95, 0.90);
  legSqrD->SetFillStyle(0);
  legSqrD->SetBorderSize(0);
  legSqrD->AddEntry(hSqrDSca, "Higgs", "F");
  legSqrD->AddEntry(hSqrDVec, "Vector", "F");
  legSqrD->AddEntry(hSqrDFQ3, "Q_{3}", "F");
  legSqrD->AddEntry(hSqrDFu3, "u_{3}", "F");
  legSqrD->AddEntry(hSqrDFd3, "d_{3}", "F");
  legSqrD->AddEntry(hSqrDFxy, "Others", "F");

  THStack* hRelD = new THStack();
  hRelD->SetTitle("Relative D factor;Overlap size;D factor fraction");
  hRelD->Add(hRelDFxy);
  hRelD->Add(hRelDFd3);
  hRelD->Add(hRelDFu3);
  hRelD->Add(hRelDFQ3);
  hRelD->Add(hRelDVec);
  hRelD->Add(hRelDSca);
  TLegend* legRelD = new TLegend(0.75, 0.65, 0.95, 0.90);
  legRelD->SetFillStyle(0);
  legRelD->SetBorderSize(0);
  legRelD->AddEntry(hRelDSca, "Higgs", "F");
  legRelD->AddEntry(hRelDVec, "Vector", "F");
  legRelD->AddEntry(hRelDFQ3, "Q_{3}", "F");
  legRelD->AddEntry(hRelDFu3, "u_{3}", "F");
  legRelD->AddEntry(hRelDFd3, "d_{3}", "F");
  legRelD->AddEntry(hRelDFxy, "Others", "F");

  TCanvas* c = new TCanvas("cProfile", "cProfile", 500, 500);
  //c->SetLogy();

  fSca->SetLineColor(kBlack);
  fVec->SetLineColor(kOrange+1);

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

  TLegend* leg = new TLegend(0.17, 0.75, 0.70, 0.92);
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

  leg->AddEntry(fVec, "Vector", "l");

  TH1F* hFrame = new TH1F("hFrame", ";Coordinate in 5th dimension;Profile amplitude #chi(y)", 100, minL, maxL);
  hFrame->SetMinimum(0);
  hFrame->SetMaximum(5);
  hFrame->Draw();

  fVec->Draw("same");
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

  TCanvas* cProfIR = new TCanvas("cProfIR", "cProfIR", 500, 500);
  hProfIR->Draw();

  TCanvas* cSqrD = new TCanvas("cSqrD", "cSqrD", 500, 500);
  hSqrD->Draw();
  hSqrD->GetYaxis()->SetRangeUser(0,1);
  hSqrD->GetXaxis()->SetNdivisions(505);
  legSqrD->Draw();

  TCanvas* cRelD = new TCanvas("cRelD", "cRelD", 500, 500);
  hRelD->Draw();
  hRelD->GetYaxis()->SetRangeUser(0,1);
  hRelD->GetXaxis()->SetNdivisions(505);
  legRelD->Draw();

*/
}

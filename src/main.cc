#include "include/ConfigReader.h"
#include "include/ADDModel.h"
#include "include/RSModel.h"
#include "include/Utility.h"

#include <boost/format.hpp>

#include <iostream>
#include <vector>
#include <string>

#ifdef DEBUGROOT
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
TFile* f = new TFile("debug.root", "recreate");
TH2F* _hMJLoss = new TH2F("hMJLoss", "hMJLoss", 100, 0., 1., 100, 0., 1.);
TH1F* _hNDecay = new TH1F("hNDecay", "hNDecay", 30, 0, 30);
TH1F* _hEDecay = new TH1F("hEDecay", "hEDecay", 100, 0, 1500);
TGraph* _grpFlux[3];
TGraph* _grpTemVsPeakPos[3];
TGraph* _grpTemVsTotalFlux[3];
std::vector<TGraph*> _grpMBHHistory, _grpJBHHistory;
#endif

using namespace std;

void printUsageAndExit();

int main(int argc, char* argv[])
{
  if ( argc < 2 ) printUsageAndExit();
  for ( int i=0; i<argc; ++i )
  {
    if ( strcmp(argv[i], "-h") == 0 or strcmp(argv[i], "--help") == 0 )
    {
      printUsageAndExit();
    }
  }

  ConfigReader cfg(argv[1], argc-2, argv+2);

  const int nEvent = cfg.get<int>("nEvent", 0, 100000);

  AbsModel* model = 0;
  const string modelName = cfg.get<string>("model");
  if ( modelName == "ADD" ) model = new ADDModel(cfg);
  else if ( modelName == "RS" ) model = new RSModel(cfg);
  else
  {
    cout << "!!" << argv[0] << ": model " << modelName << " not supported" << endl;
    return 2;
  }

  cout << "\n+ Calculating cross section...\n\n";
  cout << "******************************************************\n";
  model->calculateCrossSection();
  const double xsec = model->getCrossSection();
  const double xsecErr = model->getCrossSectionError();
  printCrossSection(xsec, xsecErr);
  cout << boost::format("** Maximum weight    = %-21.5g        **\n") % model->getWeightMax();
  cout << "******************************************************\n";

#ifdef DEBUGROOT
  TDirectory* dirMBHHistory = f->mkdir("MBHHistory");
  TDirectory* dirJBHHistory = f->mkdir("JBHHistory");
  for ( int i=0; i<3; ++i )
  {
    _grpFlux[i] = new TGraph();
    _grpFlux[i]->SetName(Form("grpFlux_%d", i));
    _grpFlux[i]->SetTitle(Form("flux curve s2 = %d", i));

    _grpTemVsTotalFlux[i] = new TGraph();
    _grpTemVsTotalFlux[i]->SetName(Form("grpTemVsTotalFlux_%d", i));
    _grpTemVsTotalFlux[i]->SetTitle(Form("Temperature vs total flux s2 = %d", i));

    _grpTemVsPeakPos[i] = new TGraph();
    _grpTemVsPeakPos[i]->SetName(Form("grpTemVsPeakPos_%d", i));
    _grpTemVsPeakPos[i]->SetTitle(Form("Temperature vs peak position s2 = %d", i));
  }
#endif
  cout << "\n+ Starting to generate events...\n\n";
  model->beginJob();
  for ( int i=0; i<nEvent; ++i )
  {
#ifdef DEBUGROOT
    _grpMBHHistory.push_back(new TGraph());
    _grpJBHHistory.push_back(new TGraph());
#endif
    printEventNumber(i, nEvent);
    model->event();
#ifdef DEBUGROOT
    dirMBHHistory->cd();
    _grpMBHHistory.back()->SetName(Form("grpMBHHistory_%d", i));
    _grpMBHHistory.back()->SetTitle(Form("BH history %d;Iteration;Mass/Initial mass", i));
    _grpMBHHistory.back()->Write();
    dirJBHHistory->cd();
    _grpJBHHistory.back()->SetName(Form("grpJBHHistory_%d", i));
    _grpJBHHistory.back()->SetTitle(Form("BH history %d;Iteration;J/Initial J", i));
    _grpJBHHistory.back()->Write();
    f->cd();
#endif
  }
  model->endJob();

#ifdef DEBUGROOT
  for ( int i=0; i<3; ++i )
  {
    _grpFlux[i]->Write();
    _grpTemVsTotalFlux[i]->Write();
    _grpTemVsPeakPos[i]->Write();
  }
  f->Write();
#endif

  return 0;
}

void printUsageAndExit()
{
  cout << "bex, A MC generator for Blackhole in EXtra dimension models\n\n";
  cout << "Usage : bex [-h] configFile [item1=value1 [item2=value2 ...] ] \n\n";
  cout << "  -h --help : print help message\n\n";
  cout << "You can override configurations by adding ITEM=VALUE format\n";

  exit(1);
}

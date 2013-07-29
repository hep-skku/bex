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
TFile* f = new TFile("debug.root", "recreate");
TH2F* _hMJLoss = new TH2F("hMJLoss", "hMJLoss", 100, 0., 1., 100, 0., 1.);
TH2F* _hMJLossNoWeight = new TH2F("hMJLossNoWeight", "hMJLossNoWeight", 100, 0., 1., 100, 0., 1.);
TH2F* _hMJLossUnWeight = new TH2F("hMJLossUnWeight", "hMJLossUnWeight", 100, 0., 1., 100, 0., 1.);
#endif

using namespace std;

int main(int argc, char* argv[])
{
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
  cout << "********************************************\n";
  model->calculateCrossSection();
  const double xsec = model->getCrossSection();
  const double xsecErr = model->getCrossSectionError();
  printCrossSection(xsec, xsecErr);
  cout << boost::format("** Maximum weight = %-21.5g **\n") % model->getWeightMax();
  cout << "********************************************\n";

  cout << "\n+ Starting to generate events...\n\n";
  model->beginJob();
  for ( int i=0; i<nEvent; ++i )
  {
    printEventNumber(i, nEvent);
    model->event();
  }
  model->endJob();

#ifdef DEBUGROOT
  f->Write();
#endif

  return 0;
}


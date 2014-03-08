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
TGraph* _grpFlux;
#endif

using namespace std;

void printUsageAndExit();

int main(int argc, char* argv[])
{
#ifdef DEBUGROOT
_grpFlux = new TGraph();
#endif
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
  _grpFlux->Write();
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

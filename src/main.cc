#include "include/ConfigReader.h"
#include "include/ADDModel.h"
#include "include/RSModel.h"
#include "include/Utility.h"
#include "include/Blackhole.h"

#include <boost/format.hpp>

#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
  ConfigReader cfg(argv[1], argc-2, argv+2);
  cfg.print();

  const int nEvent = cfg.get<int>("nEvent");

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
  cout << "##############################################\n";
  model->calculateCrossSection();
  const double xsec = model->getCrossSection();
  const double xsecErr = model->getCrossSectionError();
  printCrossSection(xsec, xsecErr);
  cout << boost::format("## Maximum weight = %-23.5g ##\n") % model->getWeightMax();
  cout << "##############################################\n";

  cout << "\n+ Starting to generate events...\n\n";
  model->beginJob();
  for ( int i=0; i<nEvent; ++i )
  {
    printEventNumber(i, nEvent);
    model->event();
  }
  model->endJob();

  return 0;
}


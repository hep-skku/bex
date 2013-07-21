#include "include/ConfigReader.h"
#include "include/ADDModel.h"
#include "include/RSModel.h"
#include "include/Utility.h"
#include "include/Blackhole.h"

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

  model->calculateCrossSection();
  const double xsec = model->getCrossSection();
  const double xsecErr = model->getCrossSectionError();
  Utility::printCrossSection(xsec, xsecErr);

  for ( int i=0; i<nEvent; ++i )
  {
    model->produce();
  }

  return 0;
}


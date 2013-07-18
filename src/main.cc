#include "include/ConfigReader.h"
#include "include/ADDModel.h"
#include "include/RSModel.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
  ConfigReader cfg(argv[1], argc-2, argv+2);
  cfg.print();

  AbsModel* model = 0;
  const string modelName = cfg.get<string>("model");
  if ( modelName == "ADD" ) model = new ADDModel(cfg);
  else if ( modelName == "RS" ) model = new RSModel(cfg);

  return 0;
}


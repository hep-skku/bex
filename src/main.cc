#include "include/ConfigReader.h"

#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[])
{
  ConfigReader cfg(argv[1], argc-2, argv+2);
  cfg.print();

  return 0;
}


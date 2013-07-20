#include "include/Utility.h"
#include <iostream>

using namespace Utility;
using namespace std;

void Utility::printCrossSection(const double xsec, const double xsecErr)
{
  cout << "Cross section = " << xsec << " +- " << xsecErr << " (pb)\n";
}

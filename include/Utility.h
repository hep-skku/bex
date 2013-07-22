#ifndef Utility_H
#define Utility_H

#include <iostream>
#include <vector>

void printCrossSection(const double xsec, const double xsecErr);
std::istream& operator>>(std::istream& in, std::vector<std::pair<double, double> >& data);

#endif


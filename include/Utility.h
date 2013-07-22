#ifndef Utility_H
#define Utility_H

#include <iostream>
#include <vector>

namespace physics
{

extern const double GevToPbarn; // = 3.894e8; pbarn*GeV2;
extern const double Pi; // Usual pi
extern const double OmegaDs[]; // Surface area of unit sphere in D-1 dimension

}

void printCrossSection(const double xsec, const double xsecErr);
std::istream& operator>>(std::istream& in, std::vector<std::pair<double, double> >& data);

#endif


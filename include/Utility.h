#ifndef Utility_H
#define Utility_H

#include <iostream>
#include <vector>

namespace physics
{

extern const double GevToPbarn; // = 3.894e8; pbarn*GeV2;
extern const double Pi; // Usual pi
extern const double OmegaDs[]; // Surface area of unit sphere in D-1 dimension
extern const double kn[]; // overall common factors k_n in BH xsec calculation

double getMassByPdgId(const int pdgId);
double r0ToRs(const int nDim, const double r0);

}

void printCrossSection(const double xsec, const double xsecErr);
void printEventNumber(const int eventNumber, const int nEvent);
std::istream& operator>>(std::istream& in, std::vector<std::pair<double, double> >& data);
std::istream& operator>>(std::istream& in, std::vector<double>& data);

#endif


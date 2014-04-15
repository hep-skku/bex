#ifndef Utility_H
#define Utility_H

#include <iostream>
#include <utility>
#include <vector>

namespace physics
{

extern const double GevToPbarn; // = 3.894e8; pbarn*GeV2;
extern const double Pi; // Usual pi
extern const double OmegaDs[]; // Surface area of unit sphere in D-1 dimension
extern const double kn[]; // overall common factors k_n in BH xsec calculation

int get3ChargeByPdgId(const int pdgId);
double getMassByPdgId(const int pdgId);
double r0ToRs(const int nDim, const double r0);

void rotate(const double phi, double& x, double& y);
void boost(const double b[], double p[]);
void boost(const double bz, double& e, double& p);
}

double interpolate(const std::vector<std::pair<double, double> >& data, const double x);
void printCrossSection(const double xsec, const double xsecErr);
void printEventNumber(const int eventNumber, const int nEvent);
std::istream& operator>>(std::istream& in, std::vector<std::pair<double, double> >& data);
std::istream& operator>>(std::istream& in, std::vector<std::vector<double> >& data);
std::istream& operator>>(std::istream& in, std::vector<double>& data);
void readValues(const char* in, std::vector<double>& data);

#endif


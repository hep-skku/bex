#include "include/Utility.h"

#include <iostream>
#include <sstream>
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

using namespace std;

namespace physics
{

const double GevToPbarn = 3.894e8; // pbarn*GeV2;
const double Pi = 3.1415926L;
const double OmegaDs[] = {
  2,
  2*Pi, 4*Pi,
  2*Pi*Pi, 8./3*Pi*Pi,
  pow(Pi, 3), 16./15.*pow(Pi, 3),
  pow(Pi, 4)/3., pow(Pi, 4)/105.*32,
  pow(Pi, 5)/12.
};

double getMassByPdgId(const int pdgId)
{
  const int absPdgId = std::abs(pdgId);
  switch ( absPdgId )
  {
    // values from http://pdg.lbl.gov/2013/mcdata/mass_width_2013.mcd
    case 2212: return  9.38272046E-01; // proton
    //case   21: return               0; // gluon
    //case   22: return               0; // photon
    case   24: return      8.0385E+01; // W
    case   23: return     9.11876E+01; // Z
    case   35: return       1.259E+02; // Higgs
    case   11: return  5.10998928E-04; // e
    case   13: return 1.056583715E-01; // mu
    case   15: return     1.77682E+00; // tau
    case    1: return         4.8E-03; // d
    case    2: return         2.3E-03; // u
    case    3: return         9.5E-02; // s
    case    4: return       1.275E+00; // c
    case    5: return       4.180E+00; // b
    case    6: return      1.7307E+02; // t
  }
  return 0;
}


}

void printCrossSection(const double xsec, const double xsecErr)
{
  boost::format fmt("## Cross section = %-.3f +- %-.3f ##");
  if ( xsec > 1e3 ) cout << fmt % (xsec/1e3) % (xsecErr/1e3) << " (nb) ##\n";
  else if ( xsec > 1e-1 ) cout << fmt % xsec % xsecErr << " (pb) ##\n";
  else cout << fmt % (xsec*1e3) % (xsecErr*1e3) << " (fb) ##\n";
}

void printEventNumber(int eventNumber, const int nEvent)
{
  const static int fw = int(log10(nEvent+1)+1);
  const static std::string fmt = (boost::format("+ Producing event: %%%dd") % fw).str();
  ++eventNumber;
  if ( eventNumber <= nEvent )
  {
    const int x = int(pow(10, floor(log10(eventNumber))));
    if ( eventNumber % x != 0 ) return;
  }
  cout << boost::format(fmt.c_str()) % eventNumber << endl;
}

istream& operator>>(istream& in, std::vector<std::pair<double, double> >& data)
{
  string line;
  while ( getline(in, line) )
  {
    const size_t commentPos = line.find('#');
    if ( commentPos != string::npos ) line.erase(commentPos);
    boost::algorithm::trim(line);
    if ( line.empty() or line[0] == '#' ) continue;
    std::replace(line.begin(), line.end(), ',', ' ');

    stringstream ss(line);
    double x, y;
    ss >> x >> y;
    data.push_back(make_pair(x, y));
  }
  return in;
}


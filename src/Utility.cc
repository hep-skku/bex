#include "include/Utility.h"

#include <iostream>
#include <sstream>
#include <cmath>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

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
const double kn[] = {
  1./4./Pi                 , sqrt(2./3./Pi)            ,
  pow(3./4., 1./3.)        , 2.*pow(5., -1./4.)        ,
  pow(5.*Pi, 1./5.)        , 2.*pow(3./7.*Pi, 1./6.)   ,
  pow(105./2.*Pi*Pi, 1./7.), 2.*pow(4.*Pi*Pi/3., 1./8.)
};

int get3ChargeByPdgId(const int pdgId)
{
  const int absPdgId = std::abs(pdgId);
  if ( absPdgId == 0 /* Undefined particle */
    or (absPdgId >= 21 and absPdgId <= 23 ) /* Neutral gauge bosons */
    or absPdgId == 12 or absPdgId == 14 or absPdgId == 16 /* Neutrinos */ ) return 0;
  const int sign = pdgId/absPdgId;
  if ( absPdgId == 11 or absPdgId == 13 or absPdgId == 15 ) return -3*sign; // Leptons
  else if ( absPdgId <= 6 )
  {
    if ( absPdgId % 2 == 1 ) return -sign; // Down type quarks : -1/3 (and +1/3 for antiquarks)
    else return sign*2; // Up type quarks : +2/3 (and -2/3 for antiquarks)
  }
  if ( absPdgId == 35 ) return 0; // Higgs

  return 0;
}

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

double r0ToRs(const int nDim, const double r0)
{
  return r0*pow( (nDim-2.)*OmegaDs[nDim-2]/4./OmegaDs[nDim-3], 1./(nDim-3) );
}

void rotate(const double phi, double& x, double& y)
{
  const double newX =  cos(phi)*x + sin(phi)*y;
  const double newY = -sin(phi)*x + cos(phi)*y;
  x = newX;
  y = newY;
}

void boost(const double b[], double p4[])
{
  const double beta2 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  const double gamma = 1./sqrt(1-beta2);
  const double sFact = (gamma-1)/beta2;

  double newP[4] = {0,};
  newP[0] = gamma*(p4[0]-b[0]*p4[1]-b[1]*p4[2]-b[2]*p4[3]);
  for ( int i=1; i<=3; ++i )
  {
    newP[i] = -gamma*b[i-1]*p4[0];
    newP[i] += p4[i];
    for ( int j=1; j<=3; ++j )
    {
      newP[i] += sFact*b[i-1]*b[j-1]*p4[j];
    }
  }
  for ( int i=0; i<=3; ++i )
  {
    p4[i] = newP[i];
  }
}

void boost(const double b, double& e, double& p)
{
  const double gamma = 1./sqrt(1-b*b);

  const double newE = gamma*(e - b*p);
  const double newP = gamma*(-b*e + p);

  p = newP;
  e = newE;
}

}
// namespace physics

// Linear interpolation. Input data have to be sorted
double interpolate(const std::vector<std::pair<double, double> >& data, double x)
{
  // Binary search on 1st item of pair
  unsigned int lo = 0, hi = data.size()-1;
  if ( hi < 1 ) throw runtime_error("interpolate: need at least 2 data points");
  //Check boundedness
  double xLo = data[lo].first;
  double xHi = data[hi].first;
  if ( x == xLo ) return data[0].second;
  else if ( x == xHi ) return data.back().second;
  else if ( x < xLo ) throw std::underflow_error((boost::format("interpolate: %1% < LB(%2%)") % x % xLo).str());
  else if ( x > xHi ) throw std::overflow_error( (boost::format("interpolate: %1% > UB(%2%)") % x % xHi).str());

  while ( true )
  {
    const unsigned int curr = (hi+lo)/2;
    const double currX = data[curr].first;
    if ( x < currX ) hi = curr;
    else if ( currX <= x ) lo = curr;
    if ( hi - lo <= 1 ) break;
  }
  xLo = data[lo].first;
  xHi = data[hi].first;
  // Now xLo <= x < xHi

  // Do linear interpolation
  const double yLo = data[lo].second;
  const double yHi = data[hi].second;
  const double frac = (x-xLo)/(xHi-xLo);
  const double y = frac*(yHi-yLo) + yLo;

  return y;
}

void printCrossSection(const double xsec, const double xsecErr)
{
  boost::format fmt("** Cross section     = %-.3f +- %-.3f");
  if      ( xsec > 1e3  ) cout << fmt % (xsec/1e3) % (xsecErr/1e3) << " (nb)          **\n";
  else if ( xsec > 1e-1 ) cout << fmt % xsec % xsecErr << " (pb)          **\n";
  else if ( xsec > 1e-6 ) cout << fmt % (xsec*1e3) % (xsecErr*1e3) << " (fb)          **\n";
  else                    cout << fmt % (xsec*1e6) % (xsecErr*1e6) << " (ab)          **\n";
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

istream& operator>>(istream& in, std::vector<std::vector<double> >& data)
{
  const size_t nCol = data.size();
  if ( nCol == 0 )
  {
    std::cerr << "Data column is empty." << std::endl;
    return in;
  }

  string line;
  while ( getline(in, line) )
  {
    const size_t commentPos = line.find('#');
    if ( commentPos != string::npos ) line.erase(commentPos);
    boost::algorithm::trim(line);
    if ( line.empty() or line[0] == '#' ) continue;
    std::replace(line.begin(), line.end(), ',', ' ');

    stringstream ss(line);
    for ( size_t i=0; i<nCol; ++i )
    {
      double x;
      ss >> x;
      data[i].push_back(x);
    }
  }
  return in;
}

istream& operator>>(istream& in, std::vector<double>& data)
{
  string content;
  string line;
  while ( getline(in, line) )
  {
    const size_t commentPos = line.find('#');
    if ( commentPos != string::npos ) line.erase(commentPos);
    boost::algorithm::trim(line);
    if ( line.empty() or line[0] == '#' ) continue;
    std::replace(line.begin(), line.end(), ',', ' ');
    content += line + " ";
  }

  stringstream ss(content);
  double x;
  while ( ss >> x )
  {
    data.push_back(x);
  }
  return in;
}

void readValues(const char* in, std::vector<double>& data)
{
  std::vector<string> fields;
  boost::split(fields, in, boost::is_any_of(" ,\t\n"));
  BOOST_FOREACH(const string& s, fields)
  {
    if ( s.empty() ) continue;
    const double x = boost::lexical_cast<double>(s);
    data.push_back(x);
  }
}

int findNearest(const double x, const std::vector<double>& v)
{
  // Do binary search
  size_t lo = 0, hi=v.size()-1;
  // Special case when hitting upper bound, x == v[hi]
  // This case can appear depending on implementation of random number algorithm
  if ( x == v[hi] ) return hi;
  while ( true )
  {
    const size_t curr = (hi+lo)/2;
    const double currX = v[curr];
    if ( x < currX ) hi = curr;
    else if ( currX <= x ) lo = curr;
    if ( hi - lo <= 1 ) break;
  }

  return lo;
}

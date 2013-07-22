#include "include/Utility.h"
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

using namespace std;

void printCrossSection(const double xsec, const double xsecErr)
{
  cout << "Cross section = " << xsec << " +- " << xsecErr << " (pb)\n";
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


#include "include/ConfigReader.h"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace boost::algorithm;

ConfigReader::ConfigReader(const char* fileName, int argc, char* argv[])
{
  ifstream fin(fileName);
  if ( fin )
  {
    string line;
    while ( getline(fin, line) )
    {
      const size_t commentPos = line.find('#');
      if ( commentPos != string::npos) line.erase(commentPos);
      trim(line);
      if ( line.empty() or line[0] == '#' ) continue;
      processInputCommand(line);
    }
  }
  else
  {
    cerr << "!!ConfigReader: cannot open file " << fileName << endl;
    return;
    // FIXME: raise exception
  }

  // Read command line options. We don't allow whitespaces here
  for ( int i=0; i<argc; ++i )
  {
    processInputCommand(argv[i]);
  }
}

void ConfigReader::processInputCommand(const string line)
{
  const size_t assignPos = line.find('=');
  if ( assignPos == string::npos )
  {
    cerr << "!!ConfigReader: invalid input command \"" << line << "\". we skip this input\n";
    return;
  }

  string name = line.substr(0, assignPos);
  string value = line.substr(assignPos+1);
  trim(name);
  trim(value);

  data_[name] = value;
}

bool ConfigReader::hasOption(const std::string name) const
{
  return data_.find(name) != data_.end();
}

void ConfigReader::print() const
{
  // Print out configurations
  for ( std::map<string, string>::const_iterator key = data_.begin();
        key != data_.end(); ++key )
  {
    cout << key->first << " = " << key->second << endl;
  }
}

template<>
std::vector<double> ConfigReader::get(const std::string name) const
{
  std::vector<double> l;
  if ( !hasOption(name) ) return l;

  std::stringstream ss(data_.find(name)->second);
  double val;
  while ( ss>>val ) l.push_back(val);

  return l;
}



#include "include/ConfigReader.h"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

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

void ConfigReader::print() const
{
  // Print out configurations
  for ( std::map<string, string>::const_iterator key = data_.begin();
        key != data_.end(); ++key )
  {
    cout << key->first << " = " << key->second << endl;
  }
}



#ifndef ConfigReader_H
#define ConfigReader_H

#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <exception>

class ConfigReader
{
public:
  ConfigReader(const char* fileName, int argc, char* argv[]);

  void processInputCommand(const std::string line);

  void print() const;
  bool hasOption(const std::string name) const;
  template<typename T>
  T get(const std::string name) const
  {
    //if ( !hasOption(name) ) throw ;
    const std::string valueStr = data_.at(name);

    T value;
    std::stringstream ss(valueStr);//data_.find(name)->second);
    ss >> value;

    return value;
  }

private:
  // Case insensitive map, from http://stackoverflow.com/questions/1801892/making-mapfind-operation-case-insensitive
  struct ci_less : std::binary_function<std::string, std::string, bool>
  {
    // case-independent (ci) compare_less binary function
    struct nocase_compare : public std::binary_function<unsigned char,unsigned char,bool>
    {
      bool operator() (const unsigned char& c1, const unsigned char& c2) const
      {
          return tolower(c1) < tolower(c2);
      }
    };
    bool operator() (const std::string & s1, const std::string & s2) const
    {
      return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(), s2.end(), nocase_compare());  // comparison
    }
  };

  typedef std::map<std::string, std::string, ci_less> DataMap;
  DataMap data_;
};

template<> std::vector<double> ConfigReader::get(const std::string name) const;


#endif


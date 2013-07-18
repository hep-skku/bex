#include <string>
#include <map>

class ConfigReader
{
public:
  ConfigReader(const char* fielName, int argc, char* argv[]);
  template<typename T>

  T get(const std::string name);
  void processInputCommand(const std::string line);
  void print();

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

  std::map<std::string, std::string, ci_less> data_;
};


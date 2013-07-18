#include <string>
#include <map>

class ConfigReader
{
public:
  ConfigReader(const char* fielName, int argc, char* argv[]);
  template<typename T>
  T get(const std::string name);
  void print();

private:
  std::map<std::string, std::string> data_;
};


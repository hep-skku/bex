#include "boost/random.hpp"
#include "boost/random/mersenne_twister.hpp"

class Random
{
public:
  Random(const unsigned int seed);
  double uniform(const double min, const double max);

private:
  boost::mt19937 rnd_;

  double shift(const double x, const double min, const double max);
};


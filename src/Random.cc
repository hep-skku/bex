#include "include/Random.h"
#include <cmath>

using namespace std;

Random::Random(const unsigned int seed)
{
  rnd_.seed(seed);
}

double Random::shift(const double x, const double min, const double max) const
{
  return (x-min)/(max-min) + min;
}

double Random::uniform(const double min, const double max)
{
  return shift(rnd_(), min, max);
}

double Random::ramp(const double min, const double max)
{
  return shift(sqrt(rnd_()), min, max);
}


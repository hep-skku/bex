#include "include/Random.h"

Random::Random(const unsigned int seed)
{
  rnd_.seed(seed);
}

double Random::shift(const double x, const double min, const double max)
{
  return (x-min)/(max-min) + min;
}

double Random::uniform(const double min, const double max)
{
  return shift(rnd_(), min, max);
}


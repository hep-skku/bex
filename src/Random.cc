#include "include/Random.h"

#include <cmath>

using namespace std;

Random::Random(const unsigned int seed)
{
  rnd_.seed(seed);
  rnd_min_ = rnd_.min();
  rnd_width_ = rnd_.max()-rnd_min_;
}

double Random::shift(const double x, const double min, const double max) const
{
  return x*(max-min) + min;
}

double Random::uniform(const double min, const double max)
{
  return shift(rand(), min, max);
}

double Random::ramp(const double min, const double max)
{
  return shift(sqrt(rand()), min, max);
}

void Random::sphere(const double r, double& x, double& y, double& z)
{
  double rr = 1e9;
  while ( rr <= 1 and rr > 0.01 )
  {
    x = rand();
    y = rand();
    z = rand();
    rr = x*x + y*y + z*z;
  }
  const double scale = r/sqrt(rr);
  x *= scale;
  y *= scale;
  z *= scale;
}

int Random::pick(const int min, const int max)
{
  return int(rand()*(max+1-min));
}


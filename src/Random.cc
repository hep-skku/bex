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

int Random::pickFromCDF(const std::vector<double>& v)
{
  // Check validity of CDF : is it monolothic array?
  for ( int i=0, n=v.size()-1; i<n; ++i )
  {
    if ( v[i] > v[i+1] )
    {
      std::cerr << "Invalid CDF, array is not monolothic\n";
      return -1;
    }
  }

  const double x = uniform(0, v.back());

  // Do binary search
  unsigned int lo = 0, hi=v.size()-1;
  // Special case when hitting upper bound, x == v[hi]
  // This case can appear depending on implementation of generator
  if ( x == v[hi] ) return hi;
  while ( true )
  {
    const unsigned int curr = (hi+lo)/2;
    const double currX = v[curr];
    if ( x < currX ) hi = curr;
    else if ( currX <= x ) lo = curr;
    if ( hi - lo <= 1 ) break;
  }

  return lo;
}

int Random::pickFromHist(const std::vector<double>& v)
{
  std::vector<double> cdf(v.size()+1);
  cdf[0] = 0;
  // Make it CDF
  for ( int i=0, n=v.size(); i<n; ++i )
  {
    cdf[i+1] = cdf[i]+v[i];
  }
  return pickFromCDF(cdf);
}


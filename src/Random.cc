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
  while ( rr > 1 or rr < 0.01 )
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

int Random::pickFromCHist(const std::vector<double>& v)
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

  return find(x, v);
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
  return pickFromCHist(cdf);
}

double Random::curve(const std::vector<std::pair<double, double> >& points, const double xmin, const double xmax)
{
  // Make CDF
  // FIXME : xmin is not applied. to be implemented
  // FIXME : xmax is not true xmax. to be updated
  std::vector<double> cdf(1);
  cdf.reserve(points.size());
  cdf[0] = 0;
  for ( int i=1, n=points.size(); i<n; ++i )
  {
    const double x0 = points[i-1].first;
    const double y0 = max(0., points[i-1].second);
    const double x1 = points[i].first;
    const double y1 = max(0., points[i].second);
    const double area = (x1-x0)*(y1+y0)/2;

    cdf.push_back(cdf.back()+area);

    if ( x1 > xmax ) break;
  }

  // Generate by inverse method
  const double y = uniform(0, cdf.back());
  const size_t index = find(y, cdf);

  const double x0 = points[index].first;
  const double x1 = points[index+1].first;
  const double y0 = cdf[index];
  const double y1 = cdf[index+1];
  const double dy = y1-y0;

  // Special case if zero prob. in this range
  if ( std::abs(dy) < 1e-12  ) return x0;
  const double invSlope = (x1-x0)/dy;

  return invSlope*(y-y0) + x0;

}

size_t Random::find(const double x, const std::vector<double>& v) const
{
  // Do binary search
  size_t lo = 0, hi=v.size()-1;
  // Special case when hitting upper bound, x == v[hi]
  // This case can appear depending on implementation of random number algorithm
  if ( x == v[hi] ) return hi;
  while ( true )
  {
    const size_t curr = (hi+lo)/2;
    const double currX = v[curr];
    if ( x < currX ) hi = curr;
    else if ( currX <= x ) lo = curr;
    if ( hi - lo <= 1 ) break;
  }

  return lo;
}

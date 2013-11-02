#ifndef Random_H
#define Random_H

#include "boost/random.hpp"
#include "boost/random/mersenne_twister.hpp"
#include <map>
#include <vector>

class Random
{
public:
  Random(const unsigned int seed);

  // Uniform distribution
  double uniform(const double min, const double max);
  // Linear ramp distribution
  double ramp(const double min, const double max);
  // Sphere with radius r
  void sphere(const double r, double& x, double& y, double& z);
  // Pick an integer from range [min, max]
  int pick(const int min, const int max);
  // Pick one from list
  template<typename IteratorType>
  IteratorType pick(IteratorType begin, IteratorType end)
  {
    const size_t shift(rand()*(end-begin));
    return begin+shift;
  }
  template<typename VectorType>
  int pickFromCDF(const VectorType& v)
  {
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

private:
  boost::mt19937 rnd_;
  double rnd_min_, rnd_width_;

  inline double rand()
  {
    return (rnd_()-rnd_min_)/rnd_width_;
  }
  double shift(const double x, const double min, const double max) const;
};

#endif


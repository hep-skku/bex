#ifndef Random_H
#define Random_H

#include "boost/random.hpp"
#include "boost/random/mersenne_twister.hpp"

class Random
{
public:
  Random(const unsigned int seed);

  // Uniform distribution
  double uniform(const double min, const double max);
  // Linear ramp distribution
  double ramp(const double min, const double max);
  // Pick an integer from range [min, max]
  int pick(const int min, const int max);
  // Pick one from list
  template<typename IteratorType>
  IteratorType pick(IteratorType begin, IteratorType end)
  {
    const size_t shift(rand()*(end-begin));
    return begin+shift;
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


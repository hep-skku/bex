#include "include/Blackhole.h"

NVector::NVector()
{
  for ( int i=0; i<nDim; ++i )
  {
    data_[i] = 0.;
  }
}

Blackhole::Blackhole()
{
}

NVector& Blackhole::nMomentum()
{
  return nMomentum_;
}

NVector& Blackhole::nPosition()
{
  return nPosition_;
}


#include "include/NVector.h"

NVector::NVector()
{
  for ( int i=0; i<nDim; ++i )
  {
    data_[i] = 0.;
  }
}


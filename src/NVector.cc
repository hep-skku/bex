#include "include/NVector.h"
#include <cmath>

using namespace std;

NVector::NVector()
{
  for ( int i=0; i<nDim; ++i )
  {
    data_[i] = 0.;
  }
}

NVector::NVector(const NVector& op)
{
  *this = op;
}

double NVector::mD2()
{
  return data_[0]*data_[0] - rD2();
}

double NVector::rD2()
{
  double retVal = 0;
  for ( int i=1; i<nDim; ++i )
  {
    retVal += data_[i]*data_[i];
  }
  return retVal;
}

NVector& NVector::operator=(const NVector& op)
{
  for ( int i=0; i<nDim; ++i )
  {
    data_[i] = op.data_[i];
  }
  return *this;
}

NVector operator+(NVector op1, const NVector& op2)
{
  op1 += op2;
  return op1;
}

NVector operator-(NVector op1, const NVector& op2)
{
  op1 -= op2;
  return op1;
}

NVector& operator+=(NVector& op1, const NVector& op2)
{
  for ( int i=0; i<op2.nDim; ++i )
  {
    op1.data_[i] += op2.data_[i];
  }
  return op1;
}

NVector& operator-=(NVector& op1, const NVector& op2)
{
  for ( int i=0; i<op2.nDim; ++i )
  {
    op1.data_[i] -= op2.data_[i];
  }
  return op1;
}

NVector& operator*=(NVector& op, const double scale)
{
  for ( int i=0; i<op.nDim; ++i )
  {
    op.data_[i] *= scale;
  }
  return op;
}

NVector operator*(const double scale, NVector op)
{
  op *= scale;
  return op;
}

NVector operator*(NVector op, const double scale)
{
  op *= scale;
  return op;
}


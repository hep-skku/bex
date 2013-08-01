#include "include/NVector.h"
#include <cmath>

using namespace std;

NVector::NVector()
{
  for ( unsigned int i=0; i<nDim; ++i )
  {
    data_[i] = 0.;
  }
}

NVector::NVector(const NVector& op)
{
  *this = op;
}

void NVector::set(const double t, const double x, const double y, const double z)
{
  data_[0] = t;
  data_[1] = x;
  data_[2] = y;
  data_[3] = z;
}

void NVector::set(const int i, const double p)
{
  data_[i] = p;
}

double NVector::mD2() const
{
  return data_[0]*data_[0] - rD2();
}

double NVector::rD2() const
{
  double retVal = 0;
  for ( unsigned int i=1; i<nDim; ++i )
  {
    retVal += data_[i]*data_[i];
  }
  return retVal;
}

NVector& NVector::operator=(const NVector& op)
{
  for ( unsigned int i=0; i<nDim; ++i )
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
  for ( unsigned int i=0; i<op2.nDim; ++i )
  {
    op1.data_[i] += op2.data_[i];
  }
  return op1;
}

NVector& operator-=(NVector& op1, const NVector& op2)
{
  for ( unsigned int i=0; i<op2.nDim; ++i )
  {
    op1.data_[i] -= op2.data_[i];
  }
  return op1;
}

NVector& operator*=(NVector& op, const double scale)
{
  for ( unsigned int i=0; i<op.nDim; ++i )
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


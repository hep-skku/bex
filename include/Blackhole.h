#ifndef Blackhole_H
#define Blackhole_H

class NVector
{
public:
  NVector();
  double mass();
  const static unsigned int nDim = 11;

private:
  double data_[nDim];

};

class Blackhole
{
public:
  Blackhole();
  NVector& nMomentum();
  NVector& nPosition();

private:
  NVector nMomentum_, nPosition_;

};

#endif


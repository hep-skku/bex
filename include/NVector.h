#ifndef NVector_H
#define NVector_H

class NVector
{
public:
  NVector();
  double mass();
  const static unsigned int nDim = 11;

private:
  double data_[nDim];

};

#endif


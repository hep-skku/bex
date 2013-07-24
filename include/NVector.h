#ifndef NVector_H
#define NVector_H

#include <cmath>

class NVector
{
public:
  NVector();
  NVector(const NVector& op);
  NVector& operator=(const NVector& op);

  const static unsigned int nDim = 11;

  // Generic D dimensional variables
  double mD2(); // D dimensional mass-squared
  double mD() { return std::sqrt(mD2()); } // D dimensional mass
  double rD2(); // D dimensional momentum/distance squared
  double rD() { return std::sqrt(rD2()); } // D dimensional momentum/distance
  double p(const int i) { return data_[i]; } // Accessor for n'th coordinate

  // Variables in 3+1 D projection
  double mass() { return m(); }
  double m() { return std::sqrt(m2()); }
  double m2() { return data_[0]*data_[0] - p2(); }
  double p2() { return x()*x() + y()*y() + z()*z(); }
  double pt2() { return x()*x() + y()*y(); }
  double pt() { return sqrt(pt2()); }

  double t() { return data_[0]; }
  double x() { return data_[1]; }
  double y() { return data_[2]; }
  double z() { return data_[3]; }

  double e()  { return e(); }
  double px() { return x(); }
  double py() { return y(); }
  double pz() { return z(); }

  friend NVector& operator+=(NVector& op1, const NVector& op2);
  friend NVector& operator-=(NVector& op1, const NVector& op2);
  friend NVector  operator+(NVector op1, const NVector& op2);
  friend NVector  operator-(NVector op1, const NVector& op2);
  friend NVector& operator*=(NVector& op, const double scale);
  friend NVector  operator*(const double scale, NVector op);
  friend NVector  operator*(NVector op, const double scale);

private:
  double data_[nDim];

};

NVector& operator+=(NVector& op1, const NVector& op2);
NVector& operator-=(NVector& op1, const NVector& op2);
NVector  operator+(NVector op1, const NVector& op2);
NVector  operator-(NVector op1, const NVector& op2);
NVector& operator*=(NVector& op, const double scale);
NVector  operator*(const double scale, NVector op);
NVector  operator*(NVector op, const double scale);

#endif


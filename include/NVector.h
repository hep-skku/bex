#ifndef NVector_H
#define NVector_H

#include <cmath>
#include <ostream>

class NVector
{
public:
  NVector();
  NVector(const NVector& op);
  NVector& operator=(const NVector& op);

  const static unsigned int nDim = 11;

  // Setters
  void set(const double t = 0, const double x = 0, const double y = 0, const double z = 0);
  void set(const int i, const double p);

  // Generic D dimensional variables
  double mD2() const; // D dimensional mass-squared
  double mD()  const { return std::sqrt(mD2()); } // D dimensional mass
  double rD2() const; // D dimensional momentum/distance squared
  double rD()  const { return std::sqrt(rD2()); } // D dimensional momentum/distance
  double p (const int i) const { return data_[i]; } // Accessor for n'th coordinate

  // Variables in 3+1 D projection
  double mass() const { return m(); }
  double m()    const { return std::sqrt(m2()); }
  double m2()   const { return data_[0]*data_[0] - p2(); }
  double p2()   const { return x()*x() + y()*y() + z()*z(); }
  double pt2()  const { return x()*x() + y()*y(); }
  double pt()   const { return sqrt(pt2()); }

  double t() const { return data_[0]; }
  double x() const { return data_[1]; }
  double y() const { return data_[2]; }
  double z() const { return data_[3]; }

  double e()  const { return e(); }
  double px() const { return x(); }
  double py() const { return y(); }
  double pz() const { return z(); }

  friend NVector& operator+=(NVector& op1, const NVector& op2);
  friend NVector& operator-=(NVector& op1, const NVector& op2);
  friend NVector  operator+(NVector op1, const NVector& op2);
  friend NVector  operator-(NVector op1, const NVector& op2);
  friend NVector& operator*=(NVector& op, const double scale);
  friend NVector  operator*(const double scale, NVector op);
  friend NVector  operator*(NVector op, const double scale);

  friend std::ostream& operator<<(std::ostream& out, const NVector& op);

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

std::ostream& operator<<(std::ostream& out, const NVector& op);

#endif


#ifndef RSModel_H
#define RSModel_H

#include "include/AbsModel.h"
#include <vector>
#include <map>

class RSModel : public AbsModel
{
public:
  RSModel(const ConfigReader& cfg);
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2);
  virtual void selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2);
  virtual bool selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           Particle& daughter);

private:
  double extraDimSize_;
  std::vector<double> nuQ_, nuU_, nuD_;
  std::map<int, double> cFactors_;
  std::map<int, double> dFactors_;
};

#endif


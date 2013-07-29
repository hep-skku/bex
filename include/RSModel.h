#ifndef RSModel_H
#define RSModel_H

#include "include/AbsModel.h"

class RSModel : public AbsModel
{
public:
  RSModel(const ConfigReader& cfg);
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2);
  virtual void selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2);
  virtual bool selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           int& dau_id, double& dau_energy);

private:
  std::vector<double> prodWeights_;
  double rs_wGG_, rs_wBG_, rs_wBB_;
};

#endif


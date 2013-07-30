#ifndef ADDModel_H
#define ADDModel_H

#include "include/AbsModel.h"

class ADDModel : public AbsModel
{
public:
  ADDModel(const ConfigReader& cfg);
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2);
  virtual bool selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           Particle& daughter);

private:

};

#endif


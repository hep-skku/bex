#ifndef RSModel_H
#define RSModel_H

#include "include/AbsModel.h"

class RSModel : public AbsModel
{
public:
  RSModel(const ConfigReader& cfg);
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2);

private:
  double kn_, kn2_;

  std::vector<double> prodWeights_;
  double rs_wGG_, rs_wBG_, rs_wBB_;
};

#endif


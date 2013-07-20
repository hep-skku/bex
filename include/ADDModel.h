#ifndef ADDModel_H
#define ADDModel_H

#include "include/AbsModel.h"

class ADDModel : public AbsModel
{
public:
  ADDModel(const ConfigReader& cfg);
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2);

private:
  int nDim_;

  double kn_, kn2_;
};

#endif


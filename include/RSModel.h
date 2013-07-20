#ifndef RSModel_H
#define RSModel_H

#include "include/AbsModel.h"

class RSModel : public AbsModel
{
public:
  RSModel(const ConfigReader& cfg);
  virtual double calculatePartonCrossSection();
};

#endif


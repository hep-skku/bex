#ifndef ADDModel_H
#define ADDModel_H

#include "include/AbsModel.h"

class ADDModel : public AbsModel
{
public:
  ADDModel(const ConfigReader& cfg);

private:
  int nDim_;
};

#endif


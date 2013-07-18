#ifndef AbsModel_H
#define AbsModel_H

#include "include/ConfigReader.h"
#include "include/PDFInterface.h"
#include "include/Random.h"

class AbsModel
{
public:
  AbsModel(const ConfigReader& cfg);
  virtual ~AbsModel();

protected:
  Random* rnd_;
  PDFInterface* pdf_;

  double beamEnergy_;
  double massMin_, massMax_;
  double mD_;
};

#endif


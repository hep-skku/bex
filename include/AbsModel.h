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

  void calculateCrossSection();
  virtual double calculatePartonCrossSection() = 0;

protected:
  Random* rnd_;
  PDFInterface* pdf_;

  double beamEnergy_;
  double massMin_, massMax_;
  double mD_;

  double xsec_, xsecErr_;

  const static int nXsecIter_ = 100000;
};

#endif


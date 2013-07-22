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

  struct FormFactorType { enum X { YOSHINO, FIOP, END }; };
  struct MassLossType { enum X { YOSHINO, CONST, END }; };

  void calculateCrossSection();
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2) = 0;
  double getCrossSection();
  double getCrossSectionError();

  void produce();

protected:
  bool isValid_;
  Random* rnd_;
  PDFInterface* pdf_;

  int nDim_;
  double beamEnergy_;
  double massMin_, massMax_;
  double mD_;

  double weightMax_;
  double xsec_, xsecErr_;

  const static double gevToPbarn_ = 3.894e8; //pbarn*GeV2
  const static int nXsecIter_ = 100000;
  const static double pi_ = 3.141592L;

  // Cached variables for convenience
  int formFactorType_, mLossType_;
  double s_;
  double formFactor_;
  std::vector<std::pair<double, double> > mLossTab_;
};

#endif


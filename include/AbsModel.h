#ifndef AbsModel_H
#define AbsModel_H

#include "include/ConfigReader.h"
#include "include/PDFInterface.h"
#include "include/Random.h"

#include <fstream>

struct Particle
{
  Particle(const int id, const int status,
                 const int mother1, const int mother2,
                 const double px, const double py, const double pz);

  int id_, status_;
  int mother1_, mother2_, color1_, color2_;
  double px_, py_, pz_, e_, m_;
  double vt_, spin_;
};

class AbsModel
{
public:
  AbsModel(const ConfigReader& cfg);
  virtual ~AbsModel();

  struct FormFactorType { enum X { YOSHINO, FIOP, PiR2, END }; };
  struct MassLossType { enum X { YOSHINO, LINEAR, UNIFORM, END }; };

  void calculateCrossSection();
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2) = 0;
  double getCrossSection();
  double getCrossSectionError();
  double getWeightMax() { return weightMax_; };

  virtual void selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2);

  virtual void beginJob();
  virtual void endJob();
  virtual void event();

protected:
  void loadYoshinoDataTable();
  double computeRs(const double m0);
  double computeMirr(const double m0, const double mFrac, const double jFrac);

protected:
  bool isValid_;
  std::string name_;

  std::ofstream fout_;

  ConfigReader cfg_;
  Random* rnd_;
  PDFInterface* pdf_;

  int beamIds_[2];
  int nDim_;
  double cmEnergy_;
  double massMin_, massMax_;
  double mD_;

  double weightMax_;
  double xsec_, xsecErr_;

  const static int nXsecIter_ = 100000;

  // Cached variables for convenience
  int formFactorType_, mLossType_;
  double s_;
  double bMax_, formFactor_;
  std::vector<std::pair<double, double> > mLossTab_;
  double kn_, kn2_;
};

#endif


#ifndef AbsModel_H
#define AbsModel_H

#include "include/ConfigReader.h"
#include "include/PDFInterface.h"
#include "include/Random.h"
#include "include/NVector.h"

#include <fstream>

struct Particle
{
  Particle(const int id, const int status,
           const int mother1, const int mother2,
           const double px, const double py, const double pz);
  Particle(const int id, const int status,
           const int mother1, const int mother2,
           const double energy);

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
  struct MassLossType { enum X { YOSHINO, LINEAR, UNIFORM, CONST, NONE, END }; };

  void calculateCrossSection();
  virtual double calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2) = 0;
  double getCrossSection();
  double getCrossSectionError();
  double getWeightMax() { return weightMax_; };

  // Select incident parton pair
  virtual void selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2);
  // Select decay particle via Hawking radiation.
  // Daughter particles varialbes are calculated at the BH rest frame.
  virtual bool selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           Particle& daughter);

  virtual void beginJob();
  virtual void endJob();
  virtual void event();

  typedef std::vector<std::pair<double, double> > Pairs;

protected:
  void loadYoshinoDataTable();
  double computeRs(const double m0) const;
  double computeMirr(const double m0, const double mFrac, const double jFrac) const;

protected:
  bool isValid_;
  std::string name_;

  std::ofstream fout_;

  ConfigReader cfg_;
  Random* rnd_;
  PDFInterface* pdf_;

  int beamId1_, beamId2_;
  int nDim_;
  double cmEnergy_;
  double massMin_, massMax_;
  double mD_;
  int formFactorType_, mLossType_, jLossType_;
  double jLossFactor_;

  double weightMax_;
  double xsec_, xsecErr_;

  const static int nXsecIter_ = 100000;

  // Cached variables for convenience
  double s_;
  double bMax_, formFactor_;
  Pairs mLossTab_;
  double kn_, kn2_;
};

#endif


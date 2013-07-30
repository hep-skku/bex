#include "include/ADDModel.h"
#include "include/Utility.h"

#include <cmath>

using namespace std;
using physics::Pi;

ADDModel::ADDModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  name_ += ":ADDModel";
}

double ADDModel::calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2)
{
  const double rssqr = pow(m/mD_, 2./(nDim_-3.))/mD_/mD_;
  const double u = m*m/s_;
  const double weightParton = rssqr*(1./massMax_-1./massMin_)*m*u*log(u);

  return weightParton*pdf1.getSumPDF()*pdf2.getSumPDF();
}

bool ADDModel::selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           Particle& daughter)
{
  return false;
}

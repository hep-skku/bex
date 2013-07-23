#include "include/ADDModel.h"
#include "include/Utility.h"

#include <cmath>

using namespace std;
using physics::Pi;

ADDModel::ADDModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  name_ += ":ADDModel";

  kn_ = pow(pow(2., nDim_-4.)*pow(Pi, (nDim_-7.)/2.)*tgamma((nDim_-1.)/2)/(nDim_-2.), 1./(nDim_-3));
  kn2_ = pow(pow(2., nDim_-4.)*pow(Pi, (nDim_-7.)/2.)*tgamma((nDim_-1.)/2)/(nDim_-2.), 2./(nDim_-3));
  formFactor_ = kn2_*Pi;
  if ( formFactorType_ == FormFactorType::YOSHINO )
  {
    const double bmax = mLossTab_.back().first;
    formFactor_ *= bmax*bmax;
  }
  else if ( formFactorType_ == FormFactorType::FIOP )
  {
    const double fiop = 4.*pow(1.+(nDim_-2.)*(nDim_-2.)/4., -2./(nDim_-3.));
    formFactor_ *= fiop;
  }
}

double ADDModel::calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2)
{
  const double rssqr = pow(m/mD_, 2./(nDim_-3.))/mD_/mD_;
  const double u = m*m/s_;
  const double weightParton = rssqr*(1./massMax_-1./massMin_)*m*u*log(u);

  return weightParton*pdf1.getSumPDF()*pdf2.getSumPDF();
}


#include "include/ADDModel.h"
#include <cmath>

using namespace std;

ADDModel::ADDModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  nDim_ = cfg.get<int>("dimension");

  kn_ = pow(pow(2., nDim_-4.)*pow(pi_, (nDim_-7.)/2.)*tgamma((nDim_-1.)/2)/(nDim_-2.), 1./(nDim_-3));
  kn2_ = pow(pow(2., nDim_-4.)*pow(pi_, (nDim_-7.)/2.)*tgamma((nDim_-1.)/2)/(nDim_-2.), 2./(nDim_-3));
  formFactor_ = kn2_*pi_;
  if ( formFactorType_ == FormFactorType::YOSHINO )
  {
    formFactor_ = 1; // FIXME : Implement Yoshino-Rychkov
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


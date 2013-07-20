#include "include/RSModel.h"
#include <iostream>
using namespace std;

RSModel::RSModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  kn2_ = 2./3*pi_;
  formFactor_ = kn2_*pi_;
  if ( formFactorName_ == "YOSHINO" )
  {
    formFactor_ = 1; // FIXME: Implement  Yoshino-Rychkov
  }
  else if ( formFactorName_ == "FIOP" )
  {
    const double fiop = 16./13;
    formFactor_ *= fiop;
  }

  prodWeights_ = cfg.get<std::vector<double> >("prodWeights");
  // wGG = C_gg
  rs_wGG_ = prodWeights_[0];
  // wBG = C_gQ + C_gd
  rs_wBG_ = prodWeights_[2] + prodWeights_[8];
  // wBB = C_QQ + C_Qd + C_dd
  rs_wBB_ = prodWeights_[1] + prodWeights_[7] + prodWeights_[6];
}

double RSModel::calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2)
{
  const double rssqr = m/mD_/mD_/mD_; // We will do 5D only
  const double u = m*m/s_;
  const double weightParton = rssqr*(1./massMax_-1./massMin_)*m*u*log(u);

  // Consider suppression factor due to profile in extra dimension
  // Assume partons are not polarized; 50% left and 50% right handed
  double weightProfile = 0;
  weightProfile += rs_wGG_*pdf1(21)*pdf2(21); // gluon+gluon
  weightProfile += 2*rs_wBG_*(pdf1(21)*pdf2(5) + pdf1(5)*pdf2(21)); // gluon+bottom (and bbar)
  weightProfile += 4*rs_wBB_*pdf1(5)*pdf2(5); // bottom+bottom (and bbar)

  return weightParton*weightProfile;
}


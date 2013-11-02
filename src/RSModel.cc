#include "include/RSModel.h"
#include "include/Utility.h"

#include <iostream>

using namespace std;
using physics::Pi;

RSModel::RSModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  name_ += ":RSModel";

  if ( nDim_ != 5 )
  {
    cerr << "!! RSModel: We do not support RS Blackhole model at " <<  nDim_ << " dimension\n";
    cerr << "!!          Changing dimension to 5D\n";
    nDim_ = 5;
    cfg_.set("dimension", nDim_);

    kn_ = physics::kn[nDim_-4];
    kn2_ = kn_*kn_;

    if ( mLossType_ == MassLossType::YOSHINO )
    {
      loadYoshinoDataTable(); // Reload mass loss data table since nDim_ is altered
    }

    if ( formFactorType_ == FormFactorType::YOSHINO )
    {
      // Same code in AbsModel.cc
      std::vector<double> maxBValues;
      std::ifstream fin("data/yoshino/max_b.data");
      fin >> maxBValues;
      bMax_ = physics::r0ToRs(nDim_, maxBValues[nDim_-4]);
    }
    else if ( formFactorType_ == FormFactorType::FIOP )
    {
      bMax_ = 2.*pow(1.+(nDim_-2.)*(nDim_-2.)/4., -1./(nDim_-3.));
    }
    kn_ = physics::kn[nDim_-4];
    kn2_ = kn_*kn_;
    formFactor_ = kn2_*physics::Pi*bMax_*bMax_;
  }
}

double RSModel::calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2)
{
  const double rssqr = m/mD_/mD_/mD_; // We will do 5D only
  const double u = m*m/s_;
  const double weightParton = rssqr*(1./massMax_-1./massMin_)*m*u*log(u);

  // Consider suppression factor due to profile in extra dimension
  // Assume partons are not polarized; 50% left and 50% right handed
  double weightProfile = 0.0140845*pdf1(21)*pdf2(21);
                       + 0.0021186*(pdf1(21)*pdf2(1) + pdf1(1)*pdf2(21))
                       + 0.0021403*(pdf1(21)*pdf2(2) + pdf1(2)*pdf2(21))
                       + 0.0103402*(pdf1(21)*pdf2(3) + pdf1(3)*pdf2(21))
                       + 0.0139810*(pdf1(21)*pdf2(4) + pdf1(4)*pdf2(21))
                       + 0.0279787*(pdf1(21)*pdf2(5) + pdf1(5)*pdf2(21))
                       + 0.0003192*pdf1(1)*pdf2(1) + 0.0003258*pdf1(2)*pdf2(2)
                       + 0.0075920*pdf1(3)*pdf2(3) + 0.0138812*pdf1(4)*pdf2(4)
                       + 0.0555812*pdf1(5)*pdf2(5)
                       + 0.0003225*(pdf1(1)*pdf2(2)+pdf1(2)*pdf2(1))
                       + 0.0015560*(pdf1(1)*pdf2(3)+pdf1(3)*pdf2(1))
                       + 0.0021043*(pdf1(1)*pdf2(4)+pdf1(4)*pdf2(1))
                       + 0.0042075*(pdf1(1)*pdf2(5)+pdf1(5)*pdf2(1))
                       + 0.0015719*(pdf1(2)*pdf2(3)+pdf1(3)*pdf2(2))
                       + 0.0021258*(pdf1(2)*pdf2(4)+pdf1(4)*pdf2(2))
                       + 0.0042506*(pdf1(2)*pdf2(5)+pdf1(5)*pdf2(2))
                       + 0.0102657*(pdf1(3)*pdf2(4)+pdf1(4)*pdf2(3))
                       + 0.0205396*(pdf1(3)*pdf2(5)+pdf1(5)*pdf2(3))
                       + 0.0277709*(pdf1(4)*pdf2(5)+pdf1(5)*pdf2(4));

  return weightParton*weightProfile;
}

void RSModel::selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2)
{
/*
  std::vector<double> weights;
  const double pdf_gg = pdf1(21)*pdf2(21);
  weights.push_back(prodWeights_[0]*pdf_gg); // BH production via Gluon+Gluon

  const double pdf_gb = pdf1(21)*pdf2(5);
  const double pdf_bg = pdf1(5)*pdf2(21);
  weights.push_back(prodWeights_[2]*pdf_gb); // g + B_R
  weights.push_back(prodWeights_[8]*pdf_gb); // g + b_L
  weights.push_back(prodWeights_[2]*pdf_bg); // B_R + g
  weights.push_back(prodWeights_[8]*pdf_bg); // b_L + g

  const double pdf_bb = pdf1(5)*pdf2(5);
  weights.push_back(prodWeights_[1]*pdf_bb); // B_R + B_R
  weights.push_back(prodWeights_[7]*pdf_bb); // b_L + B_R
  weights.push_back(prodWeights_[7]*pdf_bb); // B_R + b_L
  weights.push_back(prodWeights_[6]*pdf_bb); // b_L + b_L

  // Make weights be cumulative
  for ( int i=1, n=weights.size(); i<n; ++i )
  {
    weights[i] += weights[i-1];
  }

  // Select one case from this weights CDF
  const double LEFT = -1, RIGHT = 1;
  const int index = rnd_->pickFromCDF(weights);
  int id1, id2;
  double spin1, spin2;
  if ( index == 0 )
  {
    id1 = id2 = 21;
    spin1 = spin2 = 9;
  }
  else if ( index <= 2 )
  {
    id1 = 21;
    id2 = 5;
    spin1 = 9;
    spin2 = index % 2 ? RIGHT : LEFT;
  }
  else if ( index <= 4 )
  {
    id1 = 5;
    id2 = 21;
    spin1 = index % 2 ? RIGHT : LEFT;
    spin2 = 9;
  }
  else
  {
    id1 = id2 = 5;
    spin1 = index % 2 ? RIGHT : LEFT;
    spin2 = ((index+1)/2) % 2 ? RIGHT : LEFT;
  }

  parton1 = Particle(id1, -1, 1, 1, 0., 0., parton1.pz_);
  parton2 = Particle(id2, -1, 2, 2, 0., 0., parton2.pz_);
  parton1.spin_ = spin1;
  parton2.spin_ = spin2;
*/
}

bool RSModel::selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                          const int bh_charge, const double bh_spin,
                          Particle& daughter)
{
  return false;
}

#include "include/RSModel.h"
#include "include/Utility.h"

#include <iostream>

using namespace std;
using physics::Pi;

RSModel::RSModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  name_ += ":RSModel";
  extraDimSize_ = cfg.get<double>("extraDimSize");
  // nu-factors to fit Fermion Yukawa structure
  nuQ_ = cfg.get<std::vector<double> >("nuQ");
  nuU_ = cfg.get<std::vector<double> >("nuU");
  nuD_ = cfg.get<std::vector<double> >("nuD");

  // Number of dimension should be 5 in RS model. Reset all parameters to 5D
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

  // Calculate C-factors and D-factors from nuFermion parameters
  const double kL = 11.3*physics::Pi;
  cFactors_[2121] = (1-exp(-2*kL))/(2*kL); // Gluon+Gluon = 2121
  for ( int i=0; i<3; ++i )
  {
    const int idD1 = i==0 ? 2 : 2*i+1;
    const int idU1 = i==0 ? 1 : 2*i+2;

    const double nQ1 = 1+2*nuQ_[i];
    const double nU1 = 1+2*nuU_[i];
    const double nD1 = 1+2*nuD_[i];

    // C_gq factors
    const double c_gQ = sqrt(nQ1/kL/(exp(kL*nQ1)-1)) * (exp(kL*nQ1/2)-exp(-2*kL)) * 2 / (4+nQ1);
    const double c_gU = sqrt(nU1/kL/(exp(kL*nU1)-1)) * (exp(kL*nU1/2)-exp(-2*kL)) * 2 / (4+nU1);
    const double c_gD = sqrt(nD1/kL/(exp(kL*nD1)-1)) * (exp(kL*nD1/2)-exp(-2*kL)) * 2 / (4+nD1);

    cFactors_[idU1*100+21] = cFactors_[21*100+idU1] = (c_gQ + c_gU)/2;
    cFactors_[idD1*100+21] = cFactors_[21*100+idD1] = (c_gQ + c_gD)/2;

    for ( int j=0; j<3; ++j )
    {
      const int idD2 = j==0 ? 2 : 2*j+1;
      const int idU2 = j==0 ? 1 : 2*j+2;

      const double nQ2 = 1+2*nuQ_[j];
      const double nU2 = 1+2*nuU_[j];
      const double nD2 = 1+2*nuD_[j];

      // C_qq factors
      const double c_QQ = sqrt( nQ1*nQ2/(exp(kL*nQ1)-1)/(exp(kL*nQ2)-1) ) * 2 * (exp(kL/2*(nQ1+nQ2)) - exp(-2*kL)) / (4+nQ1+nQ2);
      const double c_QU = sqrt( nQ1*nU2/(exp(kL*nQ1)-1)/(exp(kL*nU2)-1) ) * 2 * (exp(kL/2*(nQ1+nU2)) - exp(-2*kL)) / (4+nQ1+nU2);
      const double c_QD = sqrt( nQ1*nD2/(exp(kL*nQ1)-1)/(exp(kL*nD2)-1) ) * 2 * (exp(kL/2*(nQ1+nD2)) - exp(-2*kL)) / (4+nQ1+nD2);
      const double c_UU = sqrt( nU1*nU2/(exp(kL*nU1)-1)/(exp(kL*nU2)-1) ) * 2 * (exp(kL/2*(nU1+nU2)) - exp(-2*kL)) / (4+nU1+nU2);
      const double c_UD = sqrt( nU1*nD2/(exp(kL*nU1)-1)/(exp(kL*nD2)-1) ) * 2 * (exp(kL/2*(nU1+nD2)) - exp(-2*kL)) / (4+nU1+nD2);
      const double c_DD = sqrt( nD1*nD2/(exp(kL*nD1)-1)/(exp(kL*nD2)-1) ) * 2 * (exp(kL/2*(nD1+nD2)) - exp(-2*kL)) / (4+nD1+nD2);

      // c_UU = (uL+uL) + (uL+uR) + (uR+uL) + (uR+uR)
      cFactors_[idU1*100 + idU2] = (c_QQ + 2*c_QU + c_UU)/4;
      // c_UD = (uL+dL) + (uL+dR) + (uR+dL) + (uR+dR)
      cFactors_[idU1*100 + idD2] = cFactors_[idU1 + 100*idD2] = (c_QQ + c_QD + c_QU + c_UD)/4;
      // c_DD = (dL+dL) + (dL+dR) + (dR+dL) + (dR+dR)
      cFactors_[idD1*100 + idD2] = cFactors_[idD1 + 100*idD2] = (c_QQ + 2*c_QD + c_DD)/4;
    }
  }
}

double RSModel::calculatePartonWeight(const double m, const PDF& pdf1, const PDF& pdf2)
{
  const double rssqr = m/mD_/mD_/mD_; // We will do 5D only
  const double u = m*m/s_;
  const double weightBH = rssqr*(1./massMax_-1./massMin_)*m*u*log(u);

  // Consider suppression factor due to profile in extra dimension
  // Assume partons are not polarized; 50% left and 50% right handed
  double weightParton = 0;
  for ( std::map<int, double>::const_iterator x = cFactors_.begin(); x != cFactors_.end(); ++x )
  {
    const int idPair = x->first;
    const double cFactor = x->second;

    const int id1 = idPair/100;
    const int id2 = idPair%100;

    double pdf = pdf1(id1)*pdf2(id2);
    if ( id1 != 21 ) pdf += pdf1(-id1)*pdf2(id2);
    if ( id2 != 21 ) pdf += pdf1(id1)*pdf2(-id2);
    if ( id1 != 21 and id2 != 21 ) pdf += pdf1(-id1)*pdf2(-id2);
    weightParton += cFactor * pdf;
  }

  return weightBH*weightParton;
}

void RSModel::selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2)
{
  std::vector<int> idPairs;
  // Change id definition little bit for full CDF, include sign
  std::vector<double> cFactorCDF;
  cFactorCDF.push_back(0);
  for ( std::map<int, double>::const_iterator x = cFactors_.begin(); x != cFactors_.end(); ++x )
  {
    const int idPair = x->first;
    const int id1 = idPair/100;
    const int id2 = idPair%100;
    const double cFactor = x->second;

    idPairs.push_back(id1*1000+id2);
    cFactorCDF.push_back(cFactorCDF.back()+cFactor*pdf1(id1)*pdf2(id2));

    if ( id1 != 21 )
    {
      idPairs.push_back((id1+1000)*1000+id2);
      cFactorCDF.push_back(cFactorCDF.back()+cFactor*pdf1(-id1)*pdf2(id2));
    }
    if ( id2 != 21 )
    {
      idPairs.push_back(id1*1000+(id2+1000));
      cFactorCDF.push_back(cFactorCDF.back()+cFactor*pdf1(id1)*pdf2(-id2));
    }
    if ( id1 != 21 and id2 != 21 )
    {
      // qbar + qbar
      idPairs.push_back((id1+1000)*1000+(id2+1000));
      cFactorCDF.push_back(cFactorCDF.back()+cFactor*pdf1(-id1)*pdf2(-id2));
    }
  }
  const int index = rnd_->pickFromCDF(cFactorCDF);
  const int idPair = idPairs[index];
  const int id1 = (idPair/1000)%100;
  const int id2 = (idPair%1000)%100;
  const int sign1 = (id1 != 21 and idPair/1000 >= 1000) ? -1 : 1;
  const int sign2 = (id2 != 21 and idPair%1000 >= 1000) ? -1 : 1;

  const double spin1 = 9., spin2 = 9.;

  parton1 = Particle(id1*sign1, -1, 1, 1, 0., 0., parton1.pz_);
  parton2 = Particle(id2*sign2, -1, 2, 2, 0., 0., parton2.pz_);
  parton1.spin_ = spin1;
  parton2.spin_ = spin2;

}

bool RSModel::selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                          const int bh_charge, const double bh_spin,
                          Particle& daughter)
{
  return false;
}

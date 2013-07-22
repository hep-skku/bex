#include "include/AbsModel.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include "include/Utility.h"

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg)
{
  isValid_ = false;

  beamEnergy_ = cfg.get<double>("beamEnergy", 0, 1e9);
  massMin_ = cfg.get<double>("massMin", 0., beamEnergy_);
  massMax_ = cfg.get<double>("massMax", massMin_, beamEnergy_);
  massMin_ = cfg.get<double>("massMin", 0., massMax_); // Check once again to check reversed input
  mD_   = cfg.get<double>("mD", 0., 1e9);
  nDim_ = cfg.get<int>("dimension", 4, 11);

  ConfigReader::MenuType formFactorMenu;
  formFactorMenu["Yoshino"] = FormFactorType::YOSHINO;
  formFactorMenu["FIOP"] = FormFactorType::FIOP;
  formFactorType_ = cfg.get("formFactor", formFactorMenu);

  if ( formFactorType_ == FormFactorType::YOSHINO and !cfg.hasOption("mLossType") )
  {
    mLossType_ = MassLossType::YOSHINO;
  }
  else
  {
    ConfigReader::MenuType mLossMenu;
    mLossMenu["Yoshino"] = MassLossType::YOSHINO;
    mLossMenu["Const"] = MassLossType::CONST;
    mLossMenu["Fixed"] = MassLossType::CONST;
    mLossType_ = cfg.get("mLossType", mLossMenu);
  }
  // Build Impact parameter vs mass loss factor tables
  if ( mLossType_ == MassLossType::CONST )
  {
    const double mLossFactor = cfg.get<double>("mLossFactor", 0, 1);
    mLossTab_.push_back(std::make_pair(0., mLossFactor));
    mLossTab_.push_back(std::make_pair(1., mLossFactor));
  }
  else if ( mLossType_ == MassLossType::YOSHINO )
  {
    loadYoshinoDataTable();
  }

  rnd_ = new Random(cfg.get<int>("seed"));
  pdf_ = new PDFInterface(cfg.get<int>("PDFSet"));

  // Calculate constants for speed up
  xsec_ = xsecErr_ = -1;
  s_ = beamEnergy_*beamEnergy_;

  isValid_ = true;
}

void AbsModel::loadYoshinoDataTable()
{
  const std::string mLossFileName = (boost::format("data/yoshino/MLB_N%1%.data") % (nDim_-4)).str();
  ifstream mLossFile(mLossFileName.c_str());
  mLossFile >> mLossTab_;

  // This Yoshino parameters are given in R0 unit.
  // Convert to Rs unit for convenience
  using physics::OmegaDs;
  const double r0ToRs = pow( (nDim_-2.)*OmegaDs[nDim_-2]/4./OmegaDs[nDim_-3], 1./(nDim_-3) );
  for ( int i=0, n=mLossTab_.size(); i<n; ++i )
  {
    const double x = mLossTab_[i].first * r0ToRs;
    const double y = mLossTab_[i].second;
    mLossTab_[i] = std::make_pair(x, y);
  }
}

AbsModel::~AbsModel()
{
  delete rnd_;
  delete pdf_;
}

double AbsModel::getCrossSection()
{
  if ( xsecErr_ < 0 ) calculateCrossSection();
  return xsec_;
}

double AbsModel::getCrossSectionError()
{
  if ( xsecErr_ < 0 ) calculateCrossSection();
  return xsecErr_;
}

void AbsModel::calculateCrossSection()
{
  if ( !isValid_ ) return;

  double sumW = 0, sumW2 = 0;
  PDF pdf1, pdf2;
  for ( int i=0; i<nXsecIter_; ++i )
  {
    // Generate BH mass
    const double m0 = 1/rnd_->uniform(1/massMax_, 1/massMin_);
    // Calculate parton level weight (overall numerical factors will be multiplied later
    const double u = m0*m0/s_;

    const double x1 = exp(rnd_->uniform(0, log(u)));
    const double x2 = u/x1;

    // Convoution with Parton distribution function
    pdf_->loadPDF(x1, m0, pdf1);
    pdf_->loadPDF(x2, m0, pdf2);

    const double weight = calculatePartonWeight(m0, pdf1, pdf2);

    if ( weightMax_ < weight ) weightMax_ = weight;
    sumW += weight;
    sumW2 += weight*weight;
  }

  // Now calculate total cross section and its error
  // Convert xsec to pico-barn unit by multiplying overall constants
  const double scale = 2*physics::GevToPbarn*formFactor_;
  xsec_ = scale*sumW/nXsecIter_;
  xsecErr_ = scale*sqrt(sumW2 - sumW*sumW/nXsecIter_)/nXsecIter_;
}

void AbsModel::produce()
{
  if ( !isValid_ ) return;

  PDF pdf1, pdf2;
  while ( true )
  {
    // Same routine in the xsec calculation
    const double m0 = 1/rnd_->uniform(1/massMax_, 1/massMin_);
    const double u = m0*m0/s_;

    const double x1 = exp(rnd_->uniform(0, log(u)));
    const double x2 = u/x1;

    pdf_->loadPDF(x1, m0, pdf1);
    pdf_->loadPDF(x2, m0, pdf2);

    const double weight = calculatePartonWeight(m0, pdf1, pdf2);
    if ( weight > rnd_->uniform(0, weightMax_) ) continue;

    // Apply mass loss
    break;
  }
}


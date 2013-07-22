#include "include/AbsModel.h"
#include <cmath>
#include <iostream>

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg)
{
  isValid_ = false;

  beamEnergy_ = cfg.get<double>("beamEnergy");
  massMin_ = cfg.get<double>("massMin");
  massMax_ = cfg.get<double>("massMax");
  mD_   = cfg.get<double>("mD");

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

  rnd_ = new Random(cfg.get<int>("seed"));
  pdf_ = new PDFInterface(cfg.get<int>("PDFSet"));

  // Calculate constants for speed up
  xsec_ = xsecErr_ = -1;
  s_ = beamEnergy_*beamEnergy_;

  isValid_ = true;
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
  const double scale = 2*gevToPbarn_*formFactor_;
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


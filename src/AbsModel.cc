#include "include/AbsModel.h"
#include <cmath>

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg)
{
  beamEnergy_ = cfg.get<double>("beamEnergy");
  massMin_ = cfg.get<double>("massMin");
  massMax_ = cfg.get<double>("massMax");
  mD_   = cfg.get<double>("mD");

  rnd_ = new Random(cfg.get<int>("seed"));
  pdf_ = new PDFInterface(cfg.get<int>("pdfId"));
}

AbsModel::~AbsModel()
{
  delete rnd_;
  delete pdf_;
}

void AbsModel::calculateCrossSection()
{
  double sumw = 0, sumw2 = 0;
  for ( int i=0; i<nXsecIter_; ++i )
  {
    double weight = calculatePartonCrossSection();

    sumw += weight;
    sumw2 += weight*weight;
  }
  xsec_ = sumw/nXsecIter_;
  xsecErr_ = sqrt(sumw2/nXsecIter_-xsec_*xsec_)/nXsecIter_;
}


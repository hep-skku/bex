#include "include/AbsModel.h"

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

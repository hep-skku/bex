#include "include/ADDModel.h"

ADDModel::ADDModel(const ConfigReader& cfg):
  AbsModel(cfg)
{
  nDim_ = cfg.get<int>("Dimension");
}

double ADDModel::calculatePartonCrossSection()
{
  return 0;
}


#include "include/PDFInterface.h"
#include "LHAPDF/LHAPDF.h"

PDFInterface::PDFInterface(const ConfigReader& cfg)
{
  pdfId_ = cfg.get<int>("id");
  LHAPDF::initPDF(pdfId_);
}

std::vector<double> PDFInterface::loadPDF(const double x, const double q) const
{
  return LHAPDF::xfx(x, q);
}


#ifndef PDFInterface_H
#define PDFInterface_H

#include "include/ConfigReader.h"
#include <vector>

class PDFInterface
{
public:
  PDFInterface(const ConfigReader& cfg);
  std::vector<double> loadPDF(const double x, const double q) const;

private:
  int pdfId_;

};

#endif


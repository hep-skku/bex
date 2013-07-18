#ifndef PDFInterface_H
#define PDFInterface_H

#include <vector>

class PDFInterface
{
public:
  PDFInterface(const int pdfId);
  std::vector<double> loadPDF(const double x, const double q) const;

private:
  int pdfId_;

};

#endif


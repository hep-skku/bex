#ifndef PDFInterface_H
#define PDFInterface_H

#include <vector>

class PDFInterface;

class PDF
{
public:
  PDF();
  double operator()(const int pdgId) const;
  double getSumPDF() const { return sumPDF_; }
  const static int nParton = 13;

private:
  void update(); // Re-calculate cached values

  double pdfValues_[nParton];
  double sumPDF_;

  friend class PDFInterface;
};

class PDFInterface
{
public:
  PDFInterface(const int pdfSet);
  void loadPDF(const double x, const double q, PDF& pdf) const;

private:
  int pdfSet_;

};

#endif


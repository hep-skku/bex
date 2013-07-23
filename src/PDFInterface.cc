#include "include/PDFInterface.h"
#include "LHAPDF/LHAPDF.h"

#include <iostream>
#include <algorithm>
#include <stdexcept>

using namespace std;

PDFInterface::PDFInterface(const int pdfSet)
{
  pdfSet_ = pdfSet;
  LHAPDF::initPDFSet(pdfSet, 0);
}

void PDFInterface::loadPDF(const double x, const double q, PDF& pdf) const
{
  LHAPDF::xfx(x, q, pdf.pdfValues_);
  for ( int i=0; i<PDF::nParton; ++i )
  {
    pdf.pdfValues_[i] /= x; 
  }
  pdf.update();
}

PDF::PDF()
{
  for ( int i=0; i<nParton; ++i )
  {
    pdfValues_[i] = 0;
  }
}

void PDF::update()
{
  sumPDF_ = 0;
  for ( int i=0; i<nParton; ++i )
  {
    if ( pdfValues_[i] < 0 ) pdfValues_[i] = 0;
    sumPDF_ += pdfValues_[i];
  }
}

double PDF::operator()(const int pdgId) const
{
  // Comment from LHAPDF/LHAPDF.h
  /// Nucleon PDF: fills primitive 13 element array pointed at by @a results with
  /// \f$ x f(x, Q) \f$ with index \f$ 0 < i < 12 \f$.
  /// @arg 0..5 = \f$ \bar{t} \f$, ..., \f$ \bar{u} \f$, \f$ \bar{d} \f$;
  /// @arg 6 = \f$ g \f$;
  /// @arg 7..12 = \f$ d \f$, \f$ u \f$, ..., \f$ t \f$.
  if ( pdgId == 21 ) return pdfValues_[6];
  if ( std::abs(pdgId) == 0 or std::abs(pdgId) > 6 ) throw out_of_range("Invalid pdgId");

  return pdfValues_[pdgId+6];
}

void PDF::getStackPDF(std::vector<double>& retVal) const
{
  retVal.resize(nParton);
  retVal[0] = pdfValues_[0];
  for ( int i=1; i<nParton; ++i )
  {
    retVal[i] = retVal[i-1]+pdfValues_[i];
  }
}

int PDF::indexToPdgId(const int index)
{
  if ( index < 0 or index > 12 ) throw out_of_range("PDF::indexToPdgId() : Invalid parton index");
  if ( index == 6 ) return 21;
  return index-6;
}

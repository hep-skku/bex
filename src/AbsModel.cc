#include "include/AbsModel.h"
#include "include/Utility.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/date_time/local_time/local_time.hpp>

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg):cfg_(cfg)
{
  isValid_ = false;
  name_ = "bex";

  beamIds_[0] = beamIds_[1] = 2212;

  cmEnergy_ = cfg.get<double>("cmEnergy", 0, 1e9);
  massMin_ = cfg.get<double>("massMin", 0., cmEnergy_);
  massMax_ = cfg.get<double>("massMax", massMin_, cmEnergy_);
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

  weightMax_ = 0;

  // Calculate constants for speed up
  xsec_ = xsecErr_ = -1;
  s_ = cmEnergy_*cmEnergy_;

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

void AbsModel::beginJob()
{
  using namespace boost::posix_time;

  // get cross section
  const double xsec = getCrossSection();
  const double xsecErr = getCrossSectionError();
  const double weightMax = getWeightMax();
  const double beamEnergy = cmEnergy_/2;

  // Open output file
  fout_.open(cfg_.get<string>("lheFile").c_str());

  fout_ << "<LesHouchesEvents version=\"1.0\">" << endl;
  fout_ << "<!--//\n";
  fout_ << boost::format("## Started %1% on %2%\n\n") % name_ % second_clock::local_time();
  cfg_.print(fout_);
  fout_ << "## Cross section = " << xsec << " +- " << xsecErr << endl;
  fout_ << "## Maximum weight during xsec calculation = " << weightMax << endl;
  fout_ << "//-->\n";
  fout_ << "<init>\n";
  // Add <init> header in LHE
  // HEPRUP run common block
  // Line 1 : IDBMUP(2) EBMUP(2) PDFGUP(2)=0 PDFSUP(2) IDWTUP=3 NPRUP=1
  //   IDBMUP(2) : Incident beam1 and beam2's PDG ID
  //   EBMUP(2)  : Incident beam1 and beam2's energy
  //   PDFGUP(2) : PDF group IDs, for backward compatibility. set to zero.
  //   PDFSUP(2) : PDF set IDs for beam1 and beam2
  //   IDWTUP    : How to set event weight. +3 for unweighted, accept all events
  //   NPRUP     : Number of user subprocess. We will consider only one subprocess
  fout_ << boost::format(" %5d %5d %12.5e %12.5e %5d %5d %5d %5d %5d %5d\n")
         % beamIds_[0] % beamIds_[1] % beamEnergy % beamEnergy
         % 0 % 0 % pdf_->getPDFSet() % pdf_->getPDFSet() % 3 % 1;
  // Line 2+ : XSECUP(NPRUP) XERRUP(NPRUP) XMAXUP(NPRUP) LPRUP(NPRUP)
  //   XSECUP(NPRUP) : Cross section of this process
  //   XERRUP(NPRUP) : Cross section error of this process
  //   XMAXUP(NPRUP) : Maximum weight of this process
  //   LPRUP(NPRUP)  : Listing of all user process IDs (set to 1 here)
  fout_ << boost::format(" %12.5e %12.5e %12.5e %5d\n")
         % xsec % xsecErr % weightMax % 1;
  fout_ << "</init>" << endl;

}

void AbsModel::endJob()
{
  // Finalize output file
  fout_ << "</LesHouchesEvents>" << endl;
  fout_.close();
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

void AbsModel::event()
{
  if ( !isValid_ ) return;

  // vectors to store decay product information
  std::vector<Particle> decays;

  // Set incoming particles
  const double beamEnergy = cmEnergy_/2;
  decays.push_back(Particle(beamIds_[0], -9, 0, 0, 0., 0., +beamEnergy));
  decays.push_back(Particle(beamIds_[1], -9, 0, 0, 0., 0., -beamEnergy));

  // Default values of Blackhole property
  //NVector bh_position, bh_momentum;
  int bh_charge = 0; // Blackhole charge
  double q2 = 0; // Initial CM energy before mass loss, Q^2

  // Start BH production
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

    // Now BH satisfies Hoop conjecture condition,
    // thus we can do initial calculations for this energy chunk
    q2 = m0*m0;
    // pick parton flavors (interface considering extension to RS model)
    Particle parton1(0, -1, 1, 1, 0., 0., +beamEnergy*x1);
    Particle parton2(0, -1, 2, 2, 0., 0., -beamEnergy*x2);
    selectParton(pdf1, pdf2, parton1, parton2);
    decays.push_back(parton1);
    decays.push_back(parton2);

    // Apply mass loss
    double mFrac = 1.0, jFrac = 1.0;
    while ( true )
    {
      // FIXME : Implement M/J loss
      break;
    }

    // Initial mass loss by 2 graviton radiation,
    // thus BH 3-momentum will be conserved

    break;
  }

  // Start evaporation

  // Remant decay

  // Put BH to event record
  fout_ << "<event>\n";
  // Header of user process common block
  // NUP IDPRUP=1 XWGTUP=1 SCALUP AQEDUP=-1 AQCDUP=-1
  fout_ << boost::format(" %5d %5d %15.10e %15.10e %15.10e\n") 
         % decays.size() % 1 % q2 % -1 % -1;
  for ( int i=0, n=decays.size(); i<n; ++i )
  {
    // User process ID, particles list
    // IDUP ISTUP MOTHUP(2) ICOLUP(2) PUP(5) VTIMUP SPINUP
    fout_ << boost::format(" %5d %5d %5d %5d %5d %5d %+15.10e %+15.10e %+15.10e %+15.10e %+15.10e %.1f %.1f\n")
          % decays[i].id_ % decays[i].status_
          % decays[i].mother1_ % decays[i].mother2_ % decays[i].color1_ % decays[i].color2_
          % decays[i].px_ % decays[i].py_ % decays[i].pz_ % decays[i].e_ % decays[i].m_
          % decays[i].vt_ % decays[i].spin_;
  }
  fout_ << "</event>" << endl;
}

void AbsModel::selectParton(const PDF& pdf1, const PDF& pdf2, Particle& parton1, Particle& parton2)
{
  std::vector<double> stackPDF1(PDF::nParton);
  std::vector<double> stackPDF2(PDF::nParton);
  pdf1.getStackPDF(stackPDF1);
  pdf2.getStackPDF(stackPDF2);
  const int id1 = PDF::indexToPdgId(rnd_->pickFromCDF(stackPDF1));
  const int id2 = PDF::indexToPdgId(rnd_->pickFromCDF(stackPDF2));
  parton1 = Particle(id1, -1, 1, 1, 0., 0., parton1.pz_);
  parton2 = Particle(id2, -1, 2, 2, 0., 0., parton2.pz_);
}

Particle::Particle(const int id, const int status,
                   const int mother1, const int mother2,
                   const double px, const double py, const double pz)
{
  id_ = id;
  status_ = status;
  mother1_ = mother1;
  mother2_ = mother2;
  color1_ = 0; // Colorless as default value
  color2_ = 0;
  px_ = px;
  py_ = py;
  pz_ = pz;
  m_ = physics::getMassByPdgId(id);
  e_ = sqrt(px*px+py*py+pz*pz + m_*m_);
  vt_ = 0;
  spin_ = 9;
}

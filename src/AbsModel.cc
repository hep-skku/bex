#include "include/AbsModel.h"
#include "include/NVector.h"
#include "include/Utility.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/date_time/local_time/local_time.hpp>
#include <boost/assign/std/vector.hpp>

#ifdef DEBUGROOT
#include "TH2F.h"
#include "TGraph.h"
extern TH2F* _hMJLoss;
extern TH1F* _hNDecay, * _hEDecay;
extern TGraph* _grpFlux[], * _grpTemVsTotalFlux[], * _grpTemVsPeakPos[];
extern std::vector<TGraph*> _grpMBHHistory;
#endif

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg):name_("bex"),cfg_(cfg)
{
  isValid_ = false;

  beamId1_ = beamId2_ = 2212;

  cmEnergy_ = cfg_.get<double>("cmEnergy", 0, 1e9);
  massMin_ = cfg_.get<double>("massMin", 0., cmEnergy_);
  massMax_ = cmEnergy_; // Let maximum value to CM energy
  mD_   = cfg_.get<double>("mD", 0., 1e9);
  nDim_ = cfg_.get<int>("dimension", 4, 11);

  // Routines for FormFactor definition
  ConfigReader::MenuType formFactorMenu;
  formFactorMenu["Yoshino"] = FormFactorType::YOSHINO;
  formFactorMenu["FIOP"] = FormFactorType::FIOP;
  formFactorMenu["PiR2"] = FormFactorType::PiR2;
  formFactorType_ = cfg_.get("formFactor", formFactorMenu);

  if ( formFactorType_ == FormFactorType::YOSHINO )
  {
    std::vector<double> maxBValues;
    std::ifstream fin("data/yoshino/max_b.data");
    fin >> maxBValues;
    bMax_ = physics::r0ToRs(nDim_, maxBValues[nDim_-4]);
  }
  else if ( formFactorType_ == FormFactorType::FIOP )
  {
    bMax_ = 2.*pow(1.+(nDim_-2.)*(nDim_-2.)/4., -1./(nDim_-3.));
  }
  else // simple PiR^2
  {
    bMax_ = 1.;
  }

  // Mass loss types
  ConfigReader::MenuType mLossMenu;
  mLossMenu["Yoshino"] = MJLossType::YOSHINO;
  mLossMenu["Uniform"] = MJLossType::UNIFORM;
  mLossMenu["Linear" ] = MJLossType::LINEAR ;
  mLossMenu["Const"  ] = MJLossType::CONST  ;
  mLossMenu["None"   ] = MJLossType::NONE   ;
  mLossType_ = cfg_.get("mLossType", mLossMenu);

  // Build Impact parameter vs mass loss factor tables
  if ( mLossType_ == MJLossType::YOSHINO )
  {
    loadYoshinoDataTable();
  }
  else
  {
    if ( mLossType_ == MJLossType::NONE ) cfg_.set("mLossFactor", 1);
    const double mLossFactor = cfg_.get<double>("mLossFactor", 0, 1);
    mLossTab_.push_back(std::make_pair(0., mLossFactor));
    mLossTab_.push_back(std::make_pair(bMax_, mLossFactor));
  }

  // Angular momentum loss
  if ( mLossType_ == MJLossType::YOSHINO and !cfg_.hasOption("jLossType") )
  {
    jLossType_ = MJLossType::YOSHINO;
    jLossFactor_ = 1.0;
  }
  else
  {
    ConfigReader::MenuType jLossMenu;
    jLossMenu["Yoshino"] = MJLossType::YOSHINO;
    jLossMenu["Uniform"] = MJLossType::UNIFORM;
    jLossMenu["Linear" ] = MJLossType::LINEAR ;
    jLossMenu["Const"  ] = MJLossType::CONST  ;
    jLossMenu["None"   ] = MJLossType::NONE   ;
    jLossType_ = cfg_.get("jLossType", jLossMenu);
    if ( jLossType_ == MJLossType::NONE ) jLossFactor_ = 1.0;
    else if ( jLossType_ == MJLossType::YOSHINO ) jLossFactor_ = 1.0;
    else jLossFactor_ = cfg_.get<double>("jLossFactor");
  }

  rnd_ = new Random(cfg_.get<int>("seed"));
  pdf_ = new PDFInterface(cfg_.get<int>("PDFSet"));

  weightMax_ = 0;

  // Calculate constants for speed up
  xsec_ = xsecErr_ = -1;
  s_ = cmEnergy_*cmEnergy_;
  kn_ = physics::kn[nDim_-4];
  kn2_ = kn_*kn_;
  formFactor_ = kn2_*physics::Pi*bMax_*bMax_;

  // Set Number of Degree of Freedom for each particles and spins
  using namespace boost::assign;
  // Bosons
  nDoF_[0] = 1;
  decayPdgIds_ += 35;
  decayNDoFs_  +=  1;
  // Spinors : charged leptons + neutrinos + quarks
  //  quarks : 3 generations with up/down types and 3 colors
  //  charged leptons : 3 generations
  //  neutrinos : 3 generations but no right handed neutrinos
  nDoF_[1] = 2*( 3*2*3 + 3 + 3./2 );
  decayPdgIds_ += 1, 2, 3, 4, 5, 6, 11, 13, 15, 12, 14, 16;
  decayNDoFs_  += 6, 6, 6, 6, 6, 6,  2,  2,  2,  1,  1,  1;
  // Vector bosons : gluons + photons + W + Z
  //  gluons = 8 : 8 combintions of bi-colors
  //  photons = 2/3. : 2 spin out of 3
  //  Z = 1 : 1 Z boson
  //  W = 2 : W+ and W-
  nDoF_[2] = 8+2/3.+1+2;
  decayPdgIds_ += 21,   22, 23, 24;
  decayNDoFs_  +=  8, 2/3.,  1,  2;

  isValid_ = true;

}

void AbsModel::loadYoshinoDataTable()
{
  const std::string mLossFileName = (boost::format("data/yoshino/MLB_N%1%.data") % (nDim_-4)).str();
  ifstream mLossFile(mLossFileName.c_str());
  mLossFile >> mLossTab_;

  // This Yoshino parameters are given in R0 unit.
  // Convert to Rs unit for convenience
  for ( int i=0, n=mLossTab_.size(); i<n; ++i )
  {
    const double x = physics::r0ToRs(nDim_, mLossTab_[i].first);
    const double y = mLossTab_[i].second;
    mLossTab_[i] = std::make_pair(x, y);
  }
}

void AbsModel::loadFluxDataTable()
{
  // Load flux data. data is stored in the data/flux/*D_*flux.dat
  const std::string nFluxFileName = (boost::format("data/flux/%1%D_Nflux.dat") % nDim_).str();
  ifstream nFluxFile(nFluxFileName.c_str());

  // Data is 3 columned data, uniform step
  std::vector<std::vector<double> > nFluxData(3);
  nFluxFile >> nFluxData;

  // Make cumulative distribution
  nFluxTabS0_.push_back(make_pair(0., nFluxData[0][0]));
  nFluxTabS1_.push_back(make_pair(0., nFluxData[1][0]));
  nFluxTabS2_.push_back(make_pair(0., nFluxData[2][0]));
  for ( int i=0, n=nFluxData[0].size(); i<n; ++i )
  {
    const double omega = 0.2*(i+1);

    nFluxTabS0_.push_back(make_pair(omega, nFluxTabS0_.back().second + nFluxData[0][i]));
    nFluxTabS1_.push_back(make_pair(omega, nFluxTabS1_.back().second + nFluxData[1][i]));
    nFluxTabS2_.push_back(make_pair(omega, nFluxTabS2_.back().second + nFluxData[2][i]));
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
  fout_ << boost::format("** Started %1% on %2%\n\n") % name_ % second_clock::local_time();
  cfg_.print(fout_);
  fout_ << "** Cross section = " << xsec << " +- " << xsecErr << endl;
  fout_ << "** Maximum weight during xsec calculation = " << weightMax << endl;
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
         % beamId1_ % beamId2_ % beamEnergy % beamEnergy
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

  cfg_.print();

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

    double weight = calculatePartonWeight(m0, pdf1, pdf2);

    // Suppression by mass and angular momentum loss
    double mFrac = 1.0;
    const double b0 = rnd_->ramp(0, bMax_);
    const double mFracMin = interpolate(mLossTab_, b0);

    // Generate mass fraction after balding phase
    if ( mLossType_ == MJLossType::UNIFORM )
    {
      mFrac = rnd_->uniform(mFracMin, 1);
    }
    else if ( mLossType_ == MJLossType::LINEAR or mLossType_ == MJLossType::YOSHINO )
    {
      mFrac = rnd_->ramp(mFracMin, 1);
    }
    else if ( mLossType_ == MJLossType::CONST or mLossType_ == MJLossType::NONE )
    {
      mFrac = mFracMin;
    }
    if ( !checkBHState(mFrac*m0) ) weight = 0;

    // Does cross section supressed by angular momentum loss? Maybe not.
    // From the hoop conjecture, we assume a BH forms only if mass of chunk exceeds the minimum mass criteria
    // Angular momentum loss especially due to the upper limit of J in low dimensions
    // are calculated AFTER BH is formed.
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
  decays.push_back(Particle(beamId1_, -9, 0, 0, 0., 0., +beamEnergy));
  decays.push_back(Particle(beamId2_, -9, 0, 0, 0., 0., -beamEnergy));

  // Default values of Blackhole property
  double bh_mass = 0, bh_spin = 0;
  NVector bh_momentum, bh_position;
  int bh_charge = 0; // Blackhole charge
  double qsqr = 0; // Initial CM energy before mass loss, Q^2

#ifdef DEBUGROOT
  TGraph* grpMBHHistory = _grpMBHHistory.back();
#endif

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
    qsqr = m0*m0;
    //const double rs0 = computeRs(m0);
    // pick parton flavors (interface considering extension to RS model)
    Particle parton1(0, -1, 1, 1, 0., 0., +beamEnergy*x1);
    Particle parton2(0, -1, 2, 2, 0., 0., -beamEnergy*x2);
    selectParton(pdf1, pdf2, parton1, parton2);
    bh_charge = physics::get3ChargeByPdgId(parton1.id_) + physics::get3ChargeByPdgId(parton2.id_);

    // Apply mass loss
    double mFrac = 1.0, jFrac = 1.0;
    const double b0 = rnd_->ramp(0, bMax_);
    const double mFracMin = interpolate(mLossTab_, b0);

    // Generate mass fraction after balding phase
    if ( mLossType_ == MJLossType::UNIFORM )
    {
      mFrac = rnd_->uniform(mFracMin, 1);
    }
    else if ( mLossType_ == MJLossType::LINEAR or mLossType_ == MJLossType::YOSHINO )
    {
      mFrac = rnd_->ramp(mFracMin, 1);
    }
    else if ( mLossType_ == MJLossType::CONST or mLossType_ == MJLossType::NONE )
    {
      mFrac = mFracMin;
    }

    // Retry if final mass is below minimum mass range
    if ( !checkBHState(mFrac*m0) ) continue;

    while ( true )
    {
      // There's upper bound of angular momentum for low dimensional cases
      // Adjust jFracMax for low dimensional cases
      double jFracMax = jLossFactor_;
      if ( jLossType_ == MJLossType::LINEAR or jLossType_ == MJLossType::UNIFORM )
      {
        if ( nDim_ == 4 ) jFracMax = min(mFrac*mFrac/b0, jLossFactor_);
        else if ( nDim_ == 5 ) jFracMax = min(pow(mFrac, 3./2.)/b0, jLossFactor_);
      }

      // Generate angular momentum fraction after balding phase
      if ( jLossType_ == MJLossType::UNIFORM )
      {
        jFrac = rnd_->uniform(0, jFracMax);
      }
      else if ( jLossType_ == MJLossType::LINEAR or jLossType_ == MJLossType::YOSHINO )
      {
        jFrac = rnd_->ramp(0, jFracMax);
      }
      else if ( jLossType_ == MJLossType::CONST or jLossType_ == MJLossType::NONE )
      {
        jFrac = jFracMax;
      }

      // BH should obey area theorem and cosmic censorship condition
      // according to the Yoshino-Rychkov, irreducible mass have to be
      // greater than minimum mass
      if ( jLossType_ == MJLossType::YOSHINO )
      {
        const double mIrr = computeMirr(mFrac, jFrac, b0);
        if ( mIrr < mFracMin ) continue;
      }
      break;
    }

    // Build initial BH
    bh_mass = mFrac*m0;
    const double bh_pz = parton1.pz_+parton2.pz_;
    //bh_momentum.set(std::sqrt(bh_mass*bh_mass+bh_pz*bh_pz), 0, 0, bh_pz);
    bh_momentum.set(std::sqrt(m0*m0+bh_pz*bh_pz), 0, 0, bh_pz);

#ifdef DEBUGROOT
    _hMJLoss->Fill(mFrac, jFrac);
#endif

    // Initial mass loss by 2 graviton radiation along +- z direction,
    // thus BH 3-momentum will be conserved
    // NOTE : This assumes isotropic gravitational radiation
    const double graviton_e = (1-mFrac)*m0/2;
    double gravPz1 =  graviton_e, gravE1 = graviton_e;
    double gravPz2 = -graviton_e, gravE2 = graviton_e;
    physics::boost(-bh_pz/bh_momentum.p(0), gravE1, gravPz1);
    physics::boost(-bh_pz/bh_momentum.p(0), gravE2, gravPz2);
    Particle graviton1(39, 1, 3, 4, 0., 0., gravPz1);
    Particle graviton2(39, 1, 3, 4, 0., 0., gravPz2);
    bh_momentum -= graviton1.p4();
    bh_momentum -= graviton2.p4();

    decays.push_back(parton1);
    decays.push_back(parton2);
    decays.push_back(graviton1);
    decays.push_back(graviton2);

    break;
  }

#ifdef DEBUGROOT
  grpMBHHistory->SetPoint(grpMBHHistory->GetN(), grpMBHHistory->GetN(), 1);
#endif
  // Start evaporation by hawking radiation
  double b[3], p4[4];
  while ( true )
  {
    // Select daughter particle
    Particle daughter(0, 1, 3, 4, 0., 0., 0.); // A dummy particle
    if ( !selectDecay(bh_momentum, bh_position, bh_charge, bh_spin, daughter) ) break;

    // Particle is selected.
    // Boost along BH momentum direction
    b[0] = -bh_momentum.p(1)/bh_momentum.p(0);
    b[1] = -bh_momentum.p(2)/bh_momentum.p(0);
    b[2] = -bh_momentum.p(3)/bh_momentum.p(0);
    p4[0] = daughter.e_; p4[1] = daughter.px_; p4[2] = daughter.py_; p4[3] = daughter.pz_;
    physics::boost(b, p4);
    daughter = Particle(daughter.id_, 1, 0, 0, p4[1], p4[2], p4[3]);

    decays.push_back(daughter);
    bh_momentum -= daughter.p4();
#ifdef DEBUGROOT
    _hEDecay->Fill(daughter.e_);
    grpMBHHistory->SetPoint(grpMBHHistory->GetN(), grpMBHHistory->GetN(), bh_momentum.mass()/bh_mass);
#endif
  }
#ifdef DEBUGROOT
  _hNDecay->Fill(decays.size()-6);
#endif

  // Remant decay

  // Apply overall phi rotations for all particles in the event
  // since we assumed BH rotation axis is parallel to x axis
  double totalpx = 0, totalpy = 0, totalpz = 0;
  const double phi = rnd_->uniform(0, 2*physics::Pi);
  for ( int i=0, n=decays.size(); i<n; ++i )
  {
    physics::rotate(phi, decays[i].px_, decays[i].py_);
    if ( decays[i].status_ == 1 )
    {
      totalpx += decays[i].px_;
      totalpy += decays[i].py_;
      totalpz += decays[i].pz_;
    }
  }

  // Put BH to event record
  fout_ << "<event>\n";
  // Header of user process common block
  // NUP IDPRUP=1 XWGTUP=1 SCALUP AQEDUP=-1 AQCDUP=-1
  fout_ << boost::format(" %5d %5d %15.10e %15.10e %15.10e\n")
         % decays.size() % 1 % qsqr % -1 % -1;
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
  std::vector<double> stackPDF1 = pdf1.getStackPDF();
  std::vector<double> stackPDF2 = pdf2.getStackPDF();
  const int id1 = PDF::indexToPdgId(rnd_->pickFromCHist(stackPDF1));
  const int id2 = PDF::indexToPdgId(rnd_->pickFromCHist(stackPDF2));
  parton1 = Particle(id1, -1, 1, 1, 0., 0., parton1.pz_);
  parton2 = Particle(id2, -1, 2, 2, 0., 0., parton2.pz_);
}

bool AbsModel::checkBHState(const double bh_mass, const int bh_charge, const double bh_spin) const
{
  if ( bh_mass < massMin_ ) return false;
  const double rh = computeRh(bh_mass, bh_spin);
  if ( rh < 0 ) return false;

  return true;
}

bool AbsModel::selectDecay(const NVector& bh_momentum, const NVector& bh_position,
                           const int bh_charge, const double bh_spin,
                           Particle& daughter)
{
  const double bh_mass2 = bh_momentum.m2();
  if ( bh_mass2 < 0 ) return false;
  const double bh_mass = sqrt(bh_mass2);
  // Stop decay if phase space is too small.
  // INFO : We applied slightly tight max energy requirement
  //        the true maximum is (bh_mass^2 - massMin^2 + particle_mass^2)/bh_mass/2.
  const double maxE = (bh_mass - massMin_*massMin_/bh_mass)/2;
  if ( maxE < 1e-5 ) return false;
  // Check BH state for safety
  if ( !checkBHState(bh_mass, bh_charge, bh_spin) ) return false;

  //const double bh_pos5 = bh_position.p(5); // Position in 5th dimension - not used in ADD model
  //const double rs = computeRs(bh_mass);
  const double rh = computeRh(bh_mass, bh_spin);
  const double astar = (nDim_-2.)/2*bh_spin/bh_mass/rh;

  // Choose particle type and its energy
  // Step 1 : pick particle spin for a given M and J, considering DoF
  //    <- we need values of g \int_0^\infty d\omega dN/d\omega
  std::vector<double> fluxes(3);
  fluxes[0] = getIntegratedFlux(0, rh, astar)*nDoF_[0]; // Scalar integrated flux * nDoF
  fluxes[1] = getIntegratedFlux(1, rh, astar)*nDoF_[1]; // Spinor integrated flux * nDoF
  fluxes[2] = getIntegratedFlux(2, rh, astar)*nDoF_[2]; // Vector integrated flux * nDoF

  while ( true )
  {
    const int selectedSpin2 = rnd_->pickFromHist(fluxes);

    // Step 2 : pick particle energy from cumulative dN/dw distribution
    //    <- assume we already have full energy flux curve.
    Pairs fluxCurve = getFluxCurve(selectedSpin2, rh, astar);
    const double energy = rnd_->curve(fluxCurve, 0, maxE);

    // Step 3 : pick specific particle
    const int id = decayPdgIds_[rnd_->pickFromHist(decayNDoFs_)];
    // Check the particle has enough energy to create daughter particle
    // If particle does not hold on shell condition, retry from the particle spin selection
    const double m = physics::getMassByPdgId(id);
    if ( energy < m ) continue;
    const double p = sqrt(energy*energy - m*m);

    // Choose particle
    // Spherical symmetric as the first implementation
    double px, py, pz;
    rnd_->sphere(p, px, py, pz);
    daughter = Particle(id, -1, 0, 0, px, py, pz);

    return true;
  }

  return false;
}

double AbsModel::computeMirr(const double mFrac, const double jFrac, const double b0) const
{
  // Implementation from CHMJLSPA() in charybdis2-1.0.3.F

  // Find horizon radius rh with satisfying the condition,
  //     rh^2 + ( (D-2)J/2/M )^2 = rh^(5-D) * rs^(D-3)
  // or, rh^(D-3) + k*rh^(D-5) - rs = 0
  // Here rh = (horizon radius of Kerr BH), rs = (schwartzschild radius)
  //
  // If we choose rs(m0)=1 unit, we can simplify equation even more,
  //     rs(M)^(D-3) = rs(mFrac*m0)^(D-3) = mFrac*rs(m0) = mFrac
  //     sqrt(k) = (D-2)J/2/M = (D-2)jFrac*(m0*b0/2)/2/(mFrac*m0) = (D-2)jFrac/4/mFrac*b0
  // Then problem becomes finding root of the following equation,
  //     x^(D-3) + k*x^(D-5) - mFrac = 0
  const double sqrtk = (nDim_-2.)/4.*jFrac/mFrac*b0;
  const double k = sqrtk*sqrtk;

  // Solve x^(D-3) + k*x^(D-5) - mFrac = 0
  double x;
  if ( nDim_ == 5 ) // 5D case : We can directly solve equation
  {
    const double det = mFrac - k;
    if ( det < 0 ) return -1;

    x = sqrt(det);
  }
  else if ( nDim_ >= 6 )
  {
    // In general, solution have to be found numerically
    // We will use Newton-Raphson method with 100 iterations
    x = 1;
    double dx = 0;
    for ( int i=0; i<100; ++i )
    {
      const double f      = pow(x, nDim_-3.) + k*pow(x, nDim_-5.) - mFrac;
      const double fprime = (nDim_-3.)*pow(x, nDim_-4.) + k*(nDim_-5.)*pow(x, nDim_-6.);
      dx = f/fprime;
      x -= dx;
    }
    // Check convergence
    if ( std::abs(dx) > 1e-4*std::abs(x) ) return -1;

    // Horizon radius must be positive definite
    if ( x < 0 ) x = 0;
  }
  else // 4D case : We can directly solve equation
  {
    // Equation is : x + k/x - mFrac = 0
    // Then the solution is: x = mFrac/2 +- sqrt( mFrac^2/4 - k )
    // where x > 0, 0 < mFrac <= 1, k >= 0
    const double det = mFrac*mFrac/4. - k;
    if ( det < 0 ) return -1;

    x = mFrac/2. + sqrt(det);
  }

  // Finallly, apply inverse function of rs^(D-2)[Mirr] = rs^(D-3)rh
  const double mirr = pow(x*mFrac, (nDim_-3.)/(nDim_-2.));
  return mirr;
}

double AbsModel::computeRs(const double m0) const
{
  return kn_*pow(m0/mD_, 1./(nDim_-3.))/mD_;
}

double AbsModel::computeRh(const double m0, const double j0) const
{
  const double rs = computeRs(m0);
  const double a = (nDim_-2)/2*j0/m0;

  // Horizon radius can be solved in closed form for 4D and 5D
  if ( nDim_ == 4 )
  {
    const double det = rs*rs - 4*a*a;
    if ( det < 0 ) return -1;

    return (rs + sqrt(det))/2;
  }
  else if ( nDim_ == 5 )
  {
    const double det = rs*rs - a*a;
    if ( det < 0 ) return -1;

    return sqrt(det);
  }
  else
  {
    // Otherwise, solve it numerically
    // We will use Newton-Raphson method with 100 iterations
    double rh = rs;
    double drh = 1e9;
    const double p = 1./(nDim_-3);
    while ( std::abs(drh) > 1e-3*std::abs(rh) )
    {
      const double det = rh*rh+a*a;
      const double f      = rh*pow(det, p) - rs*pow(rh, 2*p);
      const double fprime = pow(det, p) + 2*p*rh*rh*pow(det, p-1) - 2*p*rs*pow(rh, 2*p-1);
      drh = f/fprime;
      rh -= drh;
    }

    // Horizon radius must be positive definite
    if ( rh < 0 ) rh = 0;

    return rh;
  }

  return -1;

}

double AbsModel::getIntegratedFlux(const int spin2, const double rh, const double astar) const
{
  const double astar2 = astar*astar;
  //const double bh_omega = astar/(1+astar2)/rh;

  // Integrate fluxes from the data table and put them
  Pairs curve = getFluxCurve(spin2, rh, astar);
  double cFlux = 0;
  for ( int i=1, n=curve.size(); i<n; ++i )
  {
    const double x1 = curve[i-1].first;
    const double y1 = curve[i-1].second;
    const double x2 = curve[i+0].first;
    const double y2 = curve[i+0].second;
    const double area = (y2+y1)/2*(x2-x1);
    cFlux += area;
  }
#ifdef DEBUGROOT
  const double bh_tem = ((nDim_-3) + (nDim_-5)*astar2)/4/physics::Pi/(1+astar2)/rh;
  TGraph* grp = _grpTemVsTotalFlux[spin2];
  grp->SetPoint(grp->GetN(), bh_tem, cFlux);
#endif

  return cFlux;
}

AbsModel::Pairs AbsModel::getFluxCurve(const int spin2, const double rh, const double astar) const
{
  // A black body radiation curve
  const int signFactor = (spin2 % 2 == 0) ? -1 : 1;

  const double astar2 = astar*astar;
  const double bh_tem = ((nDim_-3) + (nDim_-5)*astar2)/4/physics::Pi/(1+astar2)/rh;
  //const double bh_Omega = astar/(1+astar2)/rh;
  const double bh_mFactor = 4*physics::Pi*astar/((nDim_-3)+(nDim_-5)*astar2); // factor in exponent : Omega/T

  Pairs fluxCurve;
  const double xmax = 2500;
  const int nPoint = 1000;
  fluxCurve.push_back(std::make_pair(0., 0.));
  for ( int i=1; i<=nPoint; ++i )
  {
    const double x = xmax/nPoint*i;
    double y = 0;
    for ( int l2=0; l2<=spin2; l2+=2 )
    {
      for ( int m2=-l2; m2<=l2; m2+=2 )
      {
        y += std::max(0., x*x/(exp(x/bh_tem - m2/2.*bh_mFactor)+signFactor));
      }
    }
    fluxCurve.push_back(std::make_pair(x, y));
  }

#ifdef DEBUGROOT
  double peakX = -1, peakY = -1e9;
  TGraph* grpFlux = _grpFlux[spin2];
  const bool doDraw = grpFlux->GetN() == 0;
  for ( int i=0; i<fluxCurve.size(); ++i )
  {
    const double x = fluxCurve[i].first;
    const double y = fluxCurve[i].second;
    if ( y > peakY )
    {
      peakX = x;
      peakY = y;
    }
    if ( doDraw ) grpFlux->SetPoint(i, x, y);
  }

  _grpTemVsPeakPos[spin2]->SetPoint(_grpTemVsPeakPos[spin2]->GetN(), bh_tem, peakX);
#endif

  return fluxCurve;
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

Particle::Particle(const int id, const int status,
                  const int mother1, const int mother2,
                  const double energy)
{
  id_ = id;
  status_ = status;
  mother1_ = mother1;
  mother2_ = mother2;
  color1_ = 0; // Colorless as default value
  color2_ = 0;
  px_ = 0;
  py_ = 0;
  m_ = physics::getMassByPdgId(id);
  e_ = energy;
  pz_ = max(0., sqrt(energy*energy - m_*m_));
  vt_ = 0;
  spin_ = 9;

}

NVector Particle::p4()
{
  NVector v;
  v.set(e_, px_, py_, pz_);
  return v;
}

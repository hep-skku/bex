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

#ifdef DEBUGROOT
#include "TH2F.h"
extern TH2F* _hMJLoss, * _hMJLossNoWeight, * _hMJLossUnWeight;
#endif

using namespace std;

AbsModel::AbsModel(const ConfigReader& cfg):name_("bex"),cfg_(cfg)
{
  isValid_ = false;

  beamId1_ = beamId2_ = 2212;

  cmEnergy_ = cfg_.get<double>("cmEnergy", 0, 1e9);
  massMin_ = cfg_.get<double>("massMin", 0., cmEnergy_);
  massMax_ = cfg_.get<double>("massMax", massMin_, cmEnergy_);
  massMin_ = cfg_.get<double>("massMin", 0., massMax_); // Check once again to check reversed input
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
  mLossMenu["Yoshino"] = MassLossType::YOSHINO;
  mLossMenu["Uniform"] = MassLossType::UNIFORM;
  mLossMenu["Linear" ] = MassLossType::LINEAR ;
  mLossMenu["Const"  ] = MassLossType::CONST  ;
  mLossMenu["None"   ] = MassLossType::NONE   ;
  mLossType_ = cfg_.get("mLossType", mLossMenu);

  // Build Impact parameter vs mass loss factor tables
  if ( mLossType_ == MassLossType::YOSHINO )
  {
    loadYoshinoDataTable();
  }
  else
  {
    if ( mLossType_ == MassLossType::NONE ) cfg_.set("mLossFactor", 1);
    const double mLossFactor = cfg_.get<double>("mLossFactor", 0, 1);
    mLossTab_.push_back(std::make_pair(0., mLossFactor));
    mLossTab_.push_back(std::make_pair(bMax_, mLossFactor));
  }

  // Angular momentum loss
  if ( mLossType_ == MassLossType::YOSHINO and !cfg_.hasOption("jLossType") )
  {
    jLossType_ = MassLossType::YOSHINO;
    jLossFactor_ = 1.0;
  }
  else
  {
    ConfigReader::MenuType jLossMenu;
    jLossMenu["Yoshino"] = MassLossType::YOSHINO;
    jLossMenu["Uniform"] = MassLossType::UNIFORM;
    jLossMenu["Linear" ] = MassLossType::LINEAR ;
    jLossMenu["Const"  ] = MassLossType::CONST  ;
    jLossMenu["None"   ] = MassLossType::NONE   ;
    jLossType_ = cfg_.get("jLossType", jLossMenu);
    if ( jLossType_ == MassLossType::NONE ) jLossFactor_ = 1.0;
    else if ( jLossType_ == MassLossType::YOSHINO ) jLossFactor_ = 1.0;
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

#ifdef DEBUGROOT
  double jFracV = 1.0;
#endif
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
    while ( true )
    {
      const double b0 = rnd_->ramp(0, bMax_);
      const double mFracMin = interpolate(mLossTab_, b0);

      // Generate mass fraction after balding phase
      if ( mLossType_ == MassLossType::UNIFORM )
      {
        mFrac = rnd_->uniform(mFracMin, 1);
      }
      else if ( mLossType_ == MassLossType::LINEAR or mLossType_ == MassLossType::YOSHINO )
      {
        mFrac = rnd_->ramp(mFracMin, 1);
      }
      else if ( mLossType_ == MassLossType::CONST or mLossType_ == MassLossType::NONE )
      {
        mFrac = mFracMin;
      }

      if ( jLossType_ == MassLossType::YOSHINO )
      {
        // Generate angular momentum fraction after balding phase
        const double jFrac = rnd_->ramp(0, 1);
        // BH should obey area theorem and cosmic censorship condition
        // according to the Yoshino-Rychkov, irreducible mass have to be
        // greater than minimum mass
        const double mIrr = computeMirr(mFrac, jFrac, b0);
        if ( mIrr < mFracMin ) continue;
#ifdef DEBUGROOT
        jFracV = jFrac;
#endif
      }

      break;
    }
    if ( mFrac*m0 < massMin_ ) weight = 0;

#ifdef DEBUGROOT
    _hMJLoss->Fill(mFrac, jFracV, weight);
    _hMJLossNoWeight->Fill(mFrac, jFracV);
#endif

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
    decays.push_back(parton1);
    decays.push_back(parton2);
    bh_charge = physics::get3ChargeByPdgId(parton1.id_) + physics::get3ChargeByPdgId(parton2.id_);

    // Apply mass loss
    double mFrac = 1.0, jFrac = 1.0;
    const double b0 = rnd_->ramp(0, bMax_);
    const double mFracMin = interpolate(mLossTab_, b0);
    while ( true )
    {
      // Generate mass fraction after balding phase
      if ( mLossType_ == MassLossType::UNIFORM )
      {
        mFrac = rnd_->uniform(mFracMin, 1);
      }
      else if ( mLossType_ == MassLossType::LINEAR or mLossType_ == MassLossType::YOSHINO )
      {
        mFrac = rnd_->ramp(mFracMin, 1);
      }
      else if ( mLossType_ == MassLossType::CONST or mLossType_ == MassLossType::NONE )
      {
        mFrac = mFracMin;
      }

      // Retry if final mass is below minimum mass range
      if ( mFrac*m0 < massMin_ ) continue;

      // There's upper bound of angular momentum for low dimensional cases
      // Adjust jFracMax for low dimensional cases
      double jFracMax = jLossFactor_;
      if ( jLossType_ == MassLossType::LINEAR or jLossType_ == MassLossType::UNIFORM )
      {
        if ( nDim_ == 4 ) jFracMax = min(mFrac*mFrac/b0, jLossFactor_);
        else if ( nDim_ == 5 ) jFracMax = min(pow(mFrac, 3./2.)/b0, jLossFactor_);
      }

      // Generate angular momentum fraction after balding phase
      if ( jLossType_ == MassLossType::UNIFORM )
      {
        jFrac = rnd_->uniform(0, jFracMax);
      }
      else if ( jLossType_ == MassLossType::LINEAR or jLossType_ == MassLossType::YOSHINO )
      {
        jFrac = rnd_->ramp(0, jFracMax);
      }
      else if ( jLossType_ == MassLossType::CONST or jLossType_ == MassLossType::NONE )
      {
        jFrac = jFracMax;
      }

      // BH should obey area theorem and cosmic censorship condition
      // according to the Yoshino-Rychkov, irreducible mass have to be
      // greater than minimum mass
      if ( jLossType_ == MassLossType::YOSHINO )
      {
        const double mIrr = computeMirr(mFrac, jFrac, b0);
        if ( mIrr < mFracMin ) continue;
      }

      break;
    }
#ifdef DEBUGROOT
    _hMJLossUnWeight->Fill(mFrac, jFrac);
#endif

    // Initial mass loss by 2 graviton radiation along +- z direction,
    // thus BH 3-momentum will be conserved
    // NOTE : This assumes isotropic gravitational radiation
    const double graviton_e = (1-mFrac)*m0/2;
    decays.push_back(Particle(39, 1, 3, 4, 0., 0., +graviton_e));
    decays.push_back(Particle(39, 1, 3, 4, 0., 0., -graviton_e));

    // Build initial BH
    bh_mass = mFrac*m0;
    const double bh_pz = parton1.pz_+parton2.pz_;
    bh_momentum.set(std::sqrt(bh_mass*bh_mass+bh_pz*bh_pz), 0, 0, bh_pz);

    break;
  }

  // Start evaporation by hawking radiation
  while ( true )
  {
    int dau_id;
    double dau_energy;
    if ( !selectDecay(bh_momentum, bh_position, bh_charge, bh_spin, dau_id, dau_energy) ) break;
  }

  // Remant decay

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
  const int id1 = PDF::indexToPdgId(rnd_->pickFromCDF(stackPDF1));
  const int id2 = PDF::indexToPdgId(rnd_->pickFromCDF(stackPDF2));
  parton1 = Particle(id1, -1, 1, 1, 0., 0., parton1.pz_);
  parton2 = Particle(id2, -1, 2, 2, 0., 0., parton2.pz_);
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

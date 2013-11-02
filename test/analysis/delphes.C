/*
root -l examples/Example2.C\(\"delphes_output.root\"\)
*/

#include "TH1.h"
#include "TSystem.h"

class ExRootResult;
class ExRootTreeReader;

void delphes(const char *inputFile, const char* outputFileName)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();
  TH1* hNJet = result->AddHist1D("hNJets", "Jet multiplicity", "Multiplicity", "Events", 20, 0., 20.);
  TH1* hNLeptons = result->AddHist1D("hNLeptons", "Lepton multiplicity", "Multiplicity", "Events", 20, 0., 20.);
  TH1* hJetPt = result->AddHist1D("hJetPt", "Jet pT", "Transverse momentum p_{T} (GeV/c)", "Entries per 10GeV/c", 50, 0., 500.);
  TH1* hNBjets = result->AddHist1D("hNBjets", "B Jet multiplicity", "B jet multiplicity", "Events", 20, 0., 20.);
  TH1* hBjetPt = result->AddHist1D("hBjetPt", "B jet pT", "Transverse momentum p_{T} (GeV/c)", "Entries per 10GeV/c", 50, 0., 500.);
  TH1* hST = result->AddHist1D("hST", "ST", "Scalar sum S_{T} (GeV/c)", "Events per 500GeV/c", 40, 0., 20000.);

  // Analyze
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  int nEventsTotal = 0, nEvents4J = 0, nEvents2B = 0, nEvents3B = 0, nEvents4B = 0;
  int nEvents4J2B = 0, nEvents6J2B = 0, nEvents6J3B = 0;

  // Loop over all events
  for(Long64_t entry = 0, allEntries = treeReader->GetEntries(); entry < allEntries; ++entry)
  {
    double sT = 0;
    int nLeptons = 0, nJets = 0, nBjets = 0;

    treeReader->ReadEntry(entry);

    for ( int i=0, n=branchMuon->GetEntriesFast(); i<n; ++i )
    {
      const Muon* muon = (Muon*)branchMuon->At(i);
      const double muonPt = muon->PT;
      if ( muonPt < 50 || fabs(muon->Eta) > 2.5 ) continue;
      sT += muonPt;
      nLeptons += 1;
    }

    for ( int i=0, n=branchElectron->GetEntriesFast(); i<n; ++i )
    {
      const Electron* electron = (Electron*)branchElectron->At(i);
      const double electronPt = electron->PT;
      if ( electronPt < 50 || fabs(electron->Eta) > 2.5 ) continue;
      sT += electronPt;
      nLeptons += 1;
    }

    for ( int i=0, n=branchJet->GetEntriesFast(); i<n; ++i )
    {
      const Jet* jet = (Jet*)branchJet->At(i);
      const double jetPt = jet->PT;
      if ( jetPt < 50 ) continue;
      if ( fabs(jet->Eta) > 2.5 ) continue;
      ++nJets;
      hJetPt->Fill(jetPt);
      sT += jetPt;
      if ( !jet->BTag ) continue;
      ++nBjets;
      hBjetPt->Fill(jetPt);
    }
    hNLeptons->Fill(nLeptons);
    hNJets->Fill(nJets);
    hNBjets->Fill(nBjets);

    // Analyse missing ET
    if( branchMissingET->GetEntriesFast() > 0 )
    {
      const MissingET* met = (MissingET*)branchMissingET->At(0);
      sT += met->MET;
    }

    hST->Fill(sT);

    ++nEventsTotal;
    if ( nJets >= 4 ) ++nEvents4J;
    if ( nBjets >= 2 ) ++nEvents2B;
    if ( nBjets >= 3 ) ++nEvents3B;
    if ( nBjets >= 4 ) ++nEvents4B;
    if ( nJets >= 4 && nBjets >= 2 ) ++nEvents4J2B;
    if ( nJets >= 6 && nBjets >= 2 ) ++nEvents6J2B;
    if ( nJets >= 6 && nBjets >= 3 ) ++nEvents6J3B;
  }

  cout << inputFile << endl;
  cout << "nJet4 " << 1.*nEvents4J/nEventsTotal << endl;
  cout << "nBjet2 " << 1.*nEvents2B/nEventsTotal << endl;
  cout << "nBjet3 " << 1.*nEvents3B/nEventsTotal << endl;
  cout << "nBjet4 " << 1.*nEvents4B/nEventsTotal << endl;
  cout << "nJet4 nBjet2 " << 1.*nEvents4J2B/nEventsTotal << endl;
  cout << "nJet6 nBjet2 " << 1.*nEvents6J2B/nEventsTotal << endl;
  cout << "nJet6 nBjet3 " << 1.*nEvents6J3B/nEventsTotal << endl;
  //result->Print("png");
  result->Write(outputFileName);

  delete result;
  delete treeReader;
  delete chain;
}


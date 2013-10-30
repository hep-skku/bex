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
  TH1* hNJet = result->AddHist1D("hNJets", "Jet multiplicity", "Jet multiplicity", "Events", 20, 0., 20.);
  TH1* hJetPt = result->AddHist1D("hJetPt", "Jet pT", "Transverse momentum p_{T} (GeV/c)", "Entries per 10GeV/c", 50, 0., 500.);
  TH1* hNBjets = result->AddHist1D("hNBjets", "B Jet multiplicity", "B jet multiplicity", "Events", 20, 0., 20.);
  TH1* hBjetPt = result->AddHist1D("hBjetPt", "B jet pT", "Transverse momentum p_{T} (GeV/c)", "Entries per 10GeV/c", 50, 0., 500.);
  TH1* hST = result->AddHist1D("hST", "ST", "Scalar sum of transverse momentum S_{T} (GeV/c)", "Events per 10GeV/c", 50, 0., 10000.);

  // Analyze
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  // Loop over all events
  for(Long64_t entry = 0, allEntries = treeReader->GetEntries(); entry < allEntries; ++entry)
  {
    double sT = 0;
    int nJets = 0, nBjets = 0;

    treeReader->ReadEntry(entry);

    for ( int i=0, n=branchMuon->GetEntriesFast(); i<n; ++i )
    {
      const Muon* muon = (Muon*)branchMuon->At(i);
      const double muonPt = muon->PT;
      if ( muonPt < 50 || fabs(muon->Eta) > 2.5 ) continue;
      sT += muonPt;
    }

    for ( int i=0, n=branchElectron->GetEntriesFast(); i<n; ++i )
    {
      const Electron* electron = (Electron*)branchElectron->At(i);
      const double electronPt = electron->PT;
      if ( electronPt < 50 || fabs(electron->Eta) > 2.5 ) continue;
      sT += electronPt;
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
    hNJets->Fill(nJets);
    hNBjets->Fill(nBjets);

    // Analyse missing ET
    if( branchMissingET->GetEntriesFast() > 0 )
    {
      const MissingET* met = (MissingET*)branchMissingET->At(0);
      sT += met->MET;
    }

    hST->Fill(sT);
  }

  //result->Print("png");
  result->Write(outputFileName);

  delete result;
  delete treeReader;
  delete chain;
}


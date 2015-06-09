#define dielectronMass_cxx
#include "dielectronMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"

using namespace std;

double dielectronMass::dR(double eta1, double phi1, double eta2, double phi2){
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

int dielectronMass::MCTruthMatch(int jele){

  int eleInd = -1;
  for(int imc = 0; imc < nMC; ++imc){

    if(fabs((*mcPID)[imc]) != 11) continue;
    if((*mcPt)[imc] < 10) continue;

    bool match_gen = dR((*mcEta)[imc], (*mcPhi)[imc], (*eleSCEta)[jele], (*elePhi)[jele]) < 0.05;
    if(match_gen && eleInd < 0) eleInd = imc;
  }
  if(eleInd >= 0){
    if(((*mcParentage)[eleInd]& 4)==0) return 1;
    else
      return 2;
  } else {
    return 3;
  }
}

bool dielectronMass::eleCut(int i){

  bool pass = true;
  if (fabs((*eleSCEta)[i]) < 1.4442) {
    if((*eleSigmaIEtaIEta)[i] > 0.011) pass = false;
  } else {
    if((*eleSigmaIEtaIEta)[i] > 0.033) pass = false;
  }
  return pass;
}

void dielectronMass::Loop()
{
  if (fChain == 0) return;
  TH1F* hMass = new TH1F("hMass", "", 100, 0, 200);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  //  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<200000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry%10000==0) cout << "Processed " << jentry
			     << "events" << endl;

    // LOOP OVER FOR ELECTRONS
    vector <int> ielectrons;
    for(int iele = 0; iele < nEle; ++iele){
      if((*elePt)[iele] < 15) continue;

      if(fabs((*eleSCEta)[iele])>2.5) continue;
      if(fabs((*eleSCEta)[iele])>1.4442 && fabs((*eleSCEta)[iele]) < 1.566) continue;

      int mc_truth = MCTruthMatch(iele);

      if(eleCut(iele)) ielectrons.push_back(iele);
    }
    if(ielectrons.size() < 2) continue;

    TLorentzVector e1, e2;
    int i = ielectrons[0];
    e1.SetPtEtaPhiM((*elePt)[i],(*eleEta)[i],(*elePhi)[i],0);
    int j = ielectrons[1];
    e2.SetPtEtaPhiM((*elePt)[j],(*eleEta)[j],(*elePhi)[j],0);

    double mass = (e1 + e2).M();

    hMass->Fill(mass);
  }
  TFile* f1 = new TFile("hist_mass.root", "recreate");
  hMass->Write();
  f1->Write();
  f1->Close();
}

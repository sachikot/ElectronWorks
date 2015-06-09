#define Wmass_cxx
#include "Wmass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <vector>
#include "TLorentzVector.h"

using namespace std;

double Wmass::dR(double eta1, double phi1, double eta2, double phi2){
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

void Wmass::Loop()
{
  if (fChain == 0) return;
  TH1F* W_mass = new TH1F("W_mass", "", 100, 0, 200);
  TH1F* MET = new TH1F("MET", "", 100, 0, 200);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry%1000==0) cout << "Processed " << jentry
			    << " events" << endl;

    if (pfMET < 40.) continue;

    vector <int> ielectrons;
    for(int iele = 0; iele < nEle; ++iele){

      if (fabs((*eleSCEta)[iele]) > 2.5) continue;
      if ( fabs((*eleSCEta)[iele]) > 1.4442 && fabs((*eleSCEta)[iele]) < 1.566 ) continue;
      if ((*elePt)[iele] < 20) continue;

      ielectrons.push_back(iele);
    }
    if(ielectrons.size()==0) continue;

    TLorentzVector ele;
    int i = ielectrons[0];

    vector <int> iphotons;
    for(int ipho = 0; ipho < nPho; ++ipho){

      if((*phoEt)[ipho] < 15.) continue;
      if(fabs((*phoEta)[ipho]) > 2.5) continue;

      iphotons.push_back(ipho);
    }
    if (iphotons.size() == 0) continue;

    double mass_w = sqrt( 2 * (*elePt)[i] * pfMET * (1 - cos(dR(0.0, (*elePhi)[i], 0.0, pfMETPhi))) );

    W_mass->Fill(mass_w);
    MET->Fill(pfMET);

  }
  TFile* f1 = new TFile("hist_13TeV_data.root", "recreate");
  W_mass->Write();
  MET->Write();

  f1->Write();
  f1->Close();
}

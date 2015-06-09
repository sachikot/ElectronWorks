#define FillHistos_cxx
#include "FillHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <vector>

using namespace std;

double FillHistos::dR(double eta1, double phi1, double eta2, double phi2){
  double dphi = acos(cos(phi1 - phi2));
  double deta = eta1 - eta2;
  return sqrt(dphi*dphi + deta*deta);
}

int FillHistos::MCTruthMatch(int jele){

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


void FillHistos::Loop()
{
  if (fChain == 0) return;
  TH1F* hPt = new TH1F("hPt", "electron pT", 100, 0, 200);
  TH1F* hEta = new TH1F("hEta", "electron eta", 100,-3, 3);
  TH1F* hSigmaIetaIeta_barrel = new TH1F("hSigmaIetaIeta_barrel", "", 100, 0, 0.05);
  TH1F* hSigmaIetaIeta_endcap = new TH1F("hSigmaIetaIeta_endcap", "", 100, 0, 0.05);

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  //  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<500000;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(jentry%10000==0) cout << "Processed " << jentry
			     << "events" << endl;

    // LOOP OVER FOR ELECTRONS
    for(int iele = 0; iele < nEle; ++iele){
      if((*elePt)[iele] < 15) continue;

      int mc_truth = MCTruthMatch(iele);

      if(mc_truth==2 || mc_truth==3){
	// Fill histograms
	hPt->Fill( (*elePt)[iele] );
	hEta->Fill( (*eleSCEta)[iele] );

	bool isBarrel = fabs( (*eleSCEta)[iele] ) < 1.479 ? true : false;
	if( isBarrel ) {
	  hSigmaIetaIeta_barrel->Fill((*eleSigmaIEtaIEta)[iele]);
	}else{
	  hSigmaIetaIeta_endcap->Fill((*eleSigmaIEtaIEta)[iele]);
	}
      }
    } // end loop over for electrons
  } // loop over for events
  //  TFile* f1 = new TFile("hist_DYtoEE_matched.root", "recreate");
  TFile* f1 = new TFile("hist_TTjets_matched.root", "recreate");
  hPt->Write();
  hEta->Write();
  hSigmaIetaIeta_barrel->Write();
  hSigmaIetaIeta_endcap->Write();

  f1->Write();
  f1->Close();
}

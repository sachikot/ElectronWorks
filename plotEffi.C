/**
AUTHOR : Lovedeep Saini, KSU, 31 Oct 2014

FUNCTION: plotRoc(TString sigFileNam, 
                  TString asFnt, TString varHist)
ARGUMENTS: 
      1. sigFileNam : TString (file containing signal histogram)
      2. asFnt    : TString (name of histogram to be used)
      3. verbose    : bool (flag to shutup the stdout)

RETURNS: void

DESCRIPTION: This one draw Effi. as a function of pt.

USAGE:  
      $ root -l
      [1]  .L plotEffi.C
      [2]  plotEffi("sigFile.root","hPt");
****************************************************************/

void plotEffi(TString sigFileNam, TString asFnt,
	     bool verbose=true){
  //Open the signal file.
  TFile* sigFile = new TFile(sigFileNam,"open");
  //Fetch the signal histos
  TH1D* hSigVar = sigFile->Get(asFnt);
  TH1D* hSigFnt = sigFile->Get(asFnt+"_iso");

  hSigFnt->Divide(hSigVar);
  hSigFnt->GetXaxis()->SetTitle("p_{T}");
  hSigFnt->SetTitle("Relative Isolation < 0.1");

  TCanvas* c = new TCanvas("c", "ROC", 500, 500);
  c->cd();
  c->SetGrid();
  //gPad->SetLogy(1);
  
  hSigFnt->Draw();

  c->Print(asFnt+".pdf");

  
}

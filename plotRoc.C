/**
AUTHOR : Lovedeep Saini, KSU, 31 Oct 2014

FUNCTION: plotRoc(TString sigFileNam, 
                  TString bkgFileNam, TString varHist)
ARGUMENTS: 
      1. sigFileNam : TString (file containing signal histogram)
      2. bkgFileNam : TString (file containing bakgd. histogram)
      3. varHist    : TString (name of histogram to be used)
      4. verbose    : bool (flag to shutup the stdout)

RETURNS: void

DESCRIPTION: This one creates ROC curves for a variable by
      analysing the signal & background efficiency for
      various values of selection threshold. It is expected
      that both input histograms have same name and exactly
      same binning pattern.

USAGE:  
      $ root -l
      [1]  .L plotRoc.C
      [2]  plotRoc("sigFile.root","bkgFile.root","sigIEtaIEta");
****************************************************************/

void plotRoc(TString sigFileNam, TString bkgFileNam, TString varHist, 
	     bool verbose=true){
  //Open the signal and background files.
  TFile* sigFile = new TFile(sigFileNam,"open");
  TFile* bkgFile = new TFile(bkgFileNam,"open");
  //Fetch the signal and background histos
  TH1D* hSig = sigFile->Get(varHist);
  TH1D* hBkg = bkgFile->Get(varHist);
  //Use any of these to extract number of bins.
  const int nRocBins = hSig->GetNbinsX();

  //Define array to store signal and bkg effis and FOM.
  double sigEff[nRocBins];
  double bkgEff[nRocBins];
  double fomPts[nRocBins];
  double varPts[nRocBins];
  //Total Signal and Bkg Strength 
  float sigDen = hSig->Integral();
  float bkgDen = hBkg->Integral();

  //A temporary variable to keep "fom" less than one
  float largest = 1;    //A bad fix indeed.
  
  //Move threshold one bin ahead per iteration, 
  // and score signal, bkg efficiencies.
  for (size_t i = 1; i <= nRocBins; i++) {
    //Signal selected upto i-th bin
    float sigNum = hSig->Integral(0, i);
    //Bkg selected upto i-th bin
    float bkgNum = hBkg->Integral(0, i);
    //Efficiency = selected/total
    sigEff[i-1]    = (sigDen > 0.0) ? (sigNum / sigDen) : 0.0;
    bkgEff[i-1]    = (bkgDen > 0.0) ? (bkgNum / bkgDen) : 0.0;
    fomPts[i-1]    = (sigNum+bkgNum > 0.0) ? (sigNum/std::sqrt(sigNum+bkgNum)) : 0.0;

    //Generate varPoints from lower bin edge and binWidth
    varPts[i-1]  = hSig->GetBinLowEdge(i)+hSig->GetBinWidth(i)/2;
    if(fomPts[i-1]>largest){largest = fomPts[i-1];}   
    //Make some noise :)
    if(verbose && sigEff[i-1]>0.60 && sigEff[i-1]<0.96 && bkgEff[i-1]<0.40)
      std::cout << "Signal eff = " << sigEff[i-1] 
		<< ";   Background eff = " << bkgEff[i-1] 
		<< ";   Cut val = " << hSig->GetBinLowEdge(i + 1) 
		<< std::endl;
  }
  
  //A loop to bring 'fom' values in range (0,1)
  for(i=1; i<=nRocBins; i++)
    fomPts[i-1] /=largest;

  //Now its time to draw ROC curve

  TCanvas* c = new TCanvas("c", "ROC", 1400, 500);
  c->Divide(3,1);
  c->cd(1);
  gPad->SetGrid();
  // gPad->SetLogy(1);
  
  /*Compare Shapes*/
  TH1D* tempSig=hSig->Clone();
  TH1D* tempBkg=hBkg->Clone();
  tempSig->Scale(1./tempSig->Integral());
  tempBkg->Scale(1./tempBkg->Integral());
  tempSig->GetXaxis()->SetTitle(varHist);
  tempSig->SetTitle("SHAPES");
  //Its good to have a look at variable.
  tempSig->Draw();
  tempBkg->Draw("same");
  tempSig->SetLineColor(kRed);
  leg = new TLegend(0.6336207,0.6525424,0.8965517,0.8771186,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(tempSig,"Signal");
  leg->AddEntry(tempBkg,"Background");
  leg->Draw();



  c->cd(2);
  gPad->SetGrid();
  TGraph* rocGraph = new TGraph(nRocBins, bkgEff, sigEff);
  rocGraph->GetXaxis()->SetTitle("BKG EFF");
  rocGraph->GetYaxis()->SetTitle("SIG EFF");
  rocGraph->GetYaxis()->SetRangeUser(0.5,1.1);
  rocGraph->GetXaxis()->SetRangeUser(0.0,0.6);
  rocGraph->SetTitle("ROC");
  rocGraph->SetLineWidth(2);
  rocGraph->SetMarkerStyle(24);
  rocGraph->SetMarkerSize(0.5);
  rocGraph->Draw("ACP");
  


  //Put SigEff-vs-Var, BkgEff-vs-Var, fom-vs-Var
  //into a single multigraph

  c->cd(3);
  gPad->SetGrid();
  
  TMultiGraph *mg     =  new TMultiGraph();
  TGraph* fomGraph    =  new TGraph(nRocBins,varPts,fomPts);
  TGraph* sigEffGraph =  new TGraph(nRocBins,varPts,sigEff);
  TGraph* bkgEffGraph =  new TGraph(nRocBins,varPts,bkgEff);

  //fomGraph->GetYaxis()->SetRangeUser(0.5,1.1);
  //fomGraph->GetXaxis()->SetRangeUser(0.0,0.6);

  mg->SetTitle("FOM");
  fomGraph->SetLineColor(kGreen);
  fomGraph->SetMarkerColor(kGreen);
  fomGraph->SetMarkerSize(0.5);
  fomGraph->SetMarkerStyle(24);
  sigEffGraph->SetLineColor(kRed);
  bkgEffGraph->SetLineColor(kBlue);
  sigEffGraph->SetLineWidth(2);
  bkgEffGraph->SetLineWidth(2);
  //sigEffGraph->SetMarkerStyle(25);
  //bkgEffGraph->SetMarkerStyle(26);
  mg->Add(fomGraph);
  mg->Add(sigEffGraph);
  mg->Add(bkgEffGraph);
  mg->Draw("ACP");
  mg->GetXaxis()->SetTitle(varHist);
  mg->GetYaxis()->SetTitle("FOM");
  leg = new TLegend(0.6336207,0.6525424,0.8965517,0.8771186,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(fomGraph,"S/#sqrt(S+B)","lp");
  leg->AddEntry(sigEffGraph,"Signal","lp");
  leg->AddEntry(bkgEffGraph,"Background","lp");
  leg->Draw();

  c->Print(varHist+".pdf");

  
}

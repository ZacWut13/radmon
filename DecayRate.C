void DecayRate(int startRun,int stopRun,int subtRun)
{
  // calibrationPeak is the uncalibrated peak for 208Tl at 2615 keV
 
  TH1::SetDefaultSumw2();

  //This is where all the built runs are
  const char* dataPath = "/Users/mjcenpa/Data/C1/Processed/built_runs";
  const char* outputPath = "/Users/mjcenpa/Data/C1/Analysis/DecayRates";

  //Putting together all the 3 hour block runs from one sample
  TChain* MGTree = new TChain("MGTree", "MGTree");
  for(int run = startRun; run <= stopRun; run++) { 
    MGTree->Add(TString::Format("%s/OR_run%d.root",dataPath,run)); 
  }

  //This draws the regular uncalibrated spectrum
  TCanvas* c1 = new TCanvas("c1","c1");
  MGTree->Draw("fEnergy>>h(1000)","fEnergy<400000");

  /* It's possible that "h" (in the FindObject)
     should be changed to "htemp" if "couldn't find h for some reason"
  */
  TH1F* h;
  h  = (TH1F*) c1->FindObject("h");
  if(h == NULL) {
    cout << "couldn't find h1 for some reason." << endl;
    return;
  }
 
  //Get fStartTime for first run (so we can use it as t=0)
  MGTree->GetEntry(0);
  MGTree->GetLeaf("fStartTime")->GetValue();
  double ti = MGTree->GetLeaf("fStartTime")->GetValue();
  
  //Get fStopTime from last run for binning purposes
  MGTree->GetEntry(MGTree->GetEntries()-1);
  MGTree->GetLeaf("fStopTime")->GetValue();
  double tf = MGTree->GetLeaf("fStopTime")->GetValue();

  TString countsVsTime;

//Calibrate the data using the calibrationPeak
  int calibrationPeak;
  calibrationPeak = 298500;
  //cout << endl << "What is the uncalibrated energy for the 208Tl peak at 2615 keV? (Usually 298500)" << endl;
  //cin >> calibrationPeak;

  cout << "Initial Gaussian fit to 2615 keV peak" << endl;
  h->Fit("gaus", "", "", calibrationPeak-3000, calibrationPeak+3000);
  TF1 *linGaus = new TF1("linGaus","gaus(0) + pol1(3)", calibrationPeak-30000, calibrationPeak+30000);
  linGaus->SetNpx(1000);
  double gpar0 = h->GetFunction("gaus")->GetParameter(0);
  double gpar1 = h->GetFunction("gaus")->GetParameter(1);
  double gpar2 = h->GetFunction("gaus")->GetParameter(2);
  cout << "Linear fit to background at 2615 keV peak" << endl;
  h->Fit("pol1","","",calibrationPeak-30000,calibrationPeak+30000);
  double lpar0 = h->GetFunction("pol1")->GetParameter(0);
  double lpar1 = h->GetFunction("pol1")->GetParameter(1);   linGaus->SetParameters(gpar0,gpar1,gpar2,lpar0,lpar1);
  cout << "Linear plus Gaussian fit to 2615 keV peak" << endl;
  h->Fit("linGaus","R");
  double calibrCentroid = h->GetFunction("linGaus")->GetParameter(1);

//Get the centroid for your signal peak
  double signalPeak = 477.6;
  double peakRange = 40;
  cout << endl << "At what energy is the signal peak (double): ";
//  cin >> signalPeak;
  cout << endl << "What is the length of the  peak range (double): ";
//  cin >> peakRange; cout <<  endl;

  MGTree->Draw(TString::Format("fEnergy*2614.5/%f >> h(10000,0,10000)",calibrCentroid)); 
  h  = (TH1F*) c1->FindObject("h");
  if(h == NULL) {
    cout << "couldn't find h2 for some reason." << endl;
    return;
  }

  cout << "Initial Gaussian fit to " << signalPeak << " keV peak" << endl;
  h->Fit("gaus", "", "", signalPeak-10, signalPeak+10);
  TF1 *linGaus2 = new TF1("linGaus2","gaus(0) + pol1(3)", signalPeak-0.5*peakRange, signalPeak+0.5*peakRange);
  linGaus2->SetNpx(1000);
  double gp0 = h->GetFunction("gaus")->GetParameter(0);
  double gp1 = h->GetFunction("gaus")->GetParameter(1);
  double gp2 = h->GetFunction("gaus")->GetParameter(2);
  cout << "Linear fit to background at " << signalPeak << " keV peak" << endl;
  h->Fit("pol1","","",signalPeak-0.5*peakRange,signalPeak+0.5*peakRange);
  double lp0 = h->GetFunction("pol1")->GetParameter(0);
  double lp1 = h->GetFunction("pol1")->GetParameter(1);
  linGaus2->SetParameters(gp0,gp1,gp2,lp0,lp1);
  cout << "Linear plus Gaussian fit to " << signalPeak << " keV peak" << endl;
  h->Fit("linGaus2","R");
  double signalCentroid = h->GetFunction("linGaus2")->GetParameter(1);

//Plot time on x axis, correctly putting the runs in chronological order
//Plot counts in specified calibrated energy range of signal on y axis
  TString calibratedRange = TString::Format("fEnergy*2614.5/%f>%f && fEnergy*2614.5/%f<%f",calibrCentroid,signalCentroid-0.25*peakRange,calibrCentroid,signalCentroid+0.25*peakRange);
  countsVsTime = TString::Format("fStartTime+fTimeStamp*1e-8-%f>>hRate(1000)",ti);
  MGTree->Draw(countsVsTime,calibratedRange);
  TH1F* hRate = (TH1F*) c1->FindObject("hRate");
  hRate->SetLineColor(kViolet);

//Works great up to this point - no errors




//For when you want background subtraction
if(subtRun != 0 )  {
    const char* bgDataPath = "/Users/mjcenpa/Data/C1/Analysis/sample_analysis";
    TFile *bg = TFile::Open(TString::Format("%s/%d_spec.root", bgDataPath, subtRun));
    TCanvas *cbg = (TCanvas*) bg->Get("c1");
    TH1D *hbg = (TH1D*) cbg->FindObject("h");
    if(hbg == NULL) {
            cout << "couldn't find hbg for some reason." << endl;
    return;
    }
    cbg->Draw();
    TH1F* h3  = (TH1F*) c1->FindObject("h");
    if(h3 == NULL) {
      cout << "couldn't find hbg2 for some reason." << endl;
      return;
    }
    TString calibratedRange = TString::Format("fEnergy*2614.5/%f>%f && fEnergy*2614.5/%f<%f",calibrCentroid,signalCentroid-0.25*peakRange,calibrCentroid,signalCentroid+0.25*peakRange);
    countsVsTime = TString::Format("fStartTime+fTimeStamp*1e-8-%f>>h3(1000)",ti);
    
cout << countsVsTime << endl;
cout << calibratedRange << endl;

    bg->Draw(countsVsTime,calibratedRange);

    TH1F* h4 = (TH1F*) c1->FindObject("h");
    if(h4 == NULL) {
            cout << "couldn't find htemp for some reason." << endl;
    return;
    }
  hRate->Add(h,-1);
  }






//Plot time on x axis, and plot counts in background energy range on y axis
  TString BGcalibratedRange = TString::Format("abs(fEnergy*2614.5/%f-%f)>%f && abs(fEnergy*2614.5/%f-%f)<%f",calibrCentroid,signalCentroid,0.25*peakRange,calibrCentroid,signalCentroid,0.5*peakRange);
  countsVsTime = TString::Format("fStartTime+fTimeStamp*1e-8-%f>>hbgRate(1000)",ti);
  MGTree->Draw(countsVsTime,BGcalibratedRange,"SAME");
  TH1F* hbgRate = (TH1F*) c1->FindObject("hbgRate");
  TH1F* hPeakRate = (TH1F*) hRate->Clone("hPeakRate");
  hPeakRate->Add(hbgRate,-1);
  hPeakRate->Draw("SAME");
  hPeakRate->SetLineColor(kRed);
  
//Plot the peak rate on a new canvas and fit with exponential
  hPeakRate->Fit("expo","","",0,20000);
  double ep0 = hPeakRate->GetFunction("expo")->GetParameter(0);
  double ep1 = hPeakRate->GetFunction("expo")->GetParameter(1);
  TF1 *expExp = new TF1("expExp","[0]*exp(-log(2)*x/[1]/60)+[2]*exp(-log(2)*x/[3]/60)",0,tf-ti);
  TF1 *expConst = new TF1("expConst","[0]*exp(-log(2)*x/[1]/60)+[2]",0,tf-ti);
  cout << "Exponential fit to counts vs time for background subtracted " << signalPeak << " keV peak" << endl;
  expExp->SetParameters(0.5*exp(ep0),-log(2)/ep1/60,0.5*exp(ep0),-log(2)/ep1/6);
  expConst->SetParameters(exp(ep0),-log(2)/ep1/60,0);
  hPeakRate->Fit("expExp","RL");

//Formatting and saving
  hRate->SetTitle(TString::Format("Runs %d-%d, %f keV",startRun,stopRun,signalPeak));
  hRate->GetXaxis()->SetTitle("Time [s]");
  hRate->GetYaxis()->SetTitle("Counts");
  hRate->GetYaxis()->SetRangeUser(-abs(2*hPeakRate->GetBinContent(hPeakRate->GetMinimumBin())),hRate->GetBinContent(hRate->GetMaximumBin())+50);
  c1->Update();
  c1->Print(TString::Format("%s/%d_%f_decayRate.root", outputPath, startRun,signalPeak).Data());
  hRate->GetXaxis()->SetRangeUser(0,50000);
  //c1->SetLogy();
  c1->Update();
  c1->Print(TString::Format("%s/%d_%f_decayRate.png", outputPath, startRun,signalPeak).Data());
  return;
}

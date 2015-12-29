int slewCorr(TH2F* hin, Double_t a, Double_t b, Double_t c){

  Int_t xb, yb, zb;
  hin->GetMaximumBin(xb, yb, zb);
  Double_t xMax = hin->GetBinCenter(xb);
  Double_t yMax = hin->GetYaxis()->GetBinCenter(yb);

  TCanvas* cSlewCorr =0;

  if( cSlewCorr) delete cSlewCorr;
  TCanvas* cSlewCorr = new TCanvas("cSlewCorr", "cSlewCorr", 0, 0, 1500, 500);
  cSlewCorr->Divide(3,1);

  /*
  Int_t nBin = hin->GetXaxis()->GetNbins();
  TH1F * g1 = new TH1F("g1", "g1", hin->GetXaxis()->GetNbins() , hin->GetXaxis()->GetXmin(), hin->GetXaxis()->GetXmax());

  for(Int_t xBin = 1; xBin < nBin; xBin ++){

    hin->ProjectionY("temp", xBin, xBin);
    Double_t xMax = hin->GetXaxis()->GetBinCenter(xBin);
    Int_t yBin = temp->GetMaximumBin();
    Double_t yMax = temp->GetBinCenter(yBin);
    Double_t yErr = temp->GetStdDev();

    //printf(" (x:y) = (%f,%f)\n", xMax, yMax);
        
    g1->Fill(xMax, yMax);
    g1->SetBinError(xBin, yErr);
  }
  /**/

  TF1* fit = new TF1("fit", "[0]/TMath::Power(x-[1],0.5) +[2]", 0.2, 1);
  Double_t par[3] = {a, b, c};
  Double_t sug[3] = {a, b, c};

  fit->SetParameters(par);

  cSlewCorr->cd(1);
  hin->Draw("colz");

  /*
  g1->Fit("fit","RQN");
  fit->GetParameters(sug);

  printf("=========== suggest parameters based on fitting: \n \ta:%6.4f, \n\tb:%6.4f, \n\tc:%7.4f\n", sug[0], sug[1], sug[2]);
  /**/
  
  printf("=========== suggest offset based on input and max : c:%7.4f\n", yMax- a/TMath::Sqrt(xMax-b));


  fit->SetParameters(par);
  fit->Draw("same");

  cSlewCorr->cd(2);
  TString hinTitle = hin->GetTitle();
  Int_t len = hinTitle.Length();
  Int_t pos = hinTitle.First(":");

  TString YStr = hinTitle;
  YStr = YStr.Remove(pos, len -pos);
  TString XStr = hinTitle;
  XStr = XStr.Remove(0, pos+1);
  
  TH2F * hout = 0;
  if( hout) delete hout;
  TH2F * hout = new TH2F("hout", "after Slew Correction", 
                         h1->GetXaxis()->GetNbins() , 
                         h1->GetXaxis()->GetXmin(), 
                         h1->GetXaxis()->GetXmax(),
                         3*(h1->GetYaxis()->GetNbins()) , 
                         h1->GetYaxis()->GetXmin(), 
                         h1->GetYaxis()->GetXmax());

  TString plotStr;
  plotStr.Form("%s- %f/TMath::Power(%s-%f,0.5)+2:%s>>hout", YStr.Data(), a, XStr.Data(), b, XStr.Data()); 
  tree->Draw(plotStr, "", "colz");

  cSlewCorr->cd(3);
  TH1D* temp = 0;
  if( temp) delete temp;
  TH1D* temp = new TH1D("temp", "temp", 
                        hout->GetYaxis()->GetNbins() , 
                        hout->GetYaxis()->GetXmin(), 
                        hout->GetYaxis()->GetXmax());

  hout->ProjectionY("temp")->Draw();

  Double_t meanY = hout->GetMean(2);
  Double_t rmsY = hout->GetRMS(2);

  TF1 * fit2 = new TF1("fit2", "gaus", hout->GetYaxis()->GetXmin(), hout->GetYaxis()->GetXmax());

  temp->Fit("fit2", "R", "", meanY - 2*rmsY, meanY + 1.0*rmsY);

  printf(" chi-squared : %f \n", fit2->GetChisquare()/fit2->GetNDF());

  return 0;

}

void slewingfit(){

  gROOT->ProcessLine(".!date");
  //gROOT->ProcessLine(".x Ana.C");

  //======================================================== InPut setting

  char * rootfile = "run1035.root";
  Int_t Div[2] = {2,2};  //x,y
  Int_t size[2] = {600,400}; //x,y

  //====================================== Load root file
  TFile *f0 = new TFile (rootfile, "read"); 
  TTree *tree = (TTree*)f0->Get("tree");
  printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
  gStyle->SetOptStat(1112211);
  gStyle->SetOptFit(111);
  //======================================================== Browser or Canvas
  TCanvas * cSlewing = new TCanvas("cSlewing", "cSlewing", 0, 0 , size[0]*Div[0], size[1]*Div[1]);
  cSlewing->Divide(Div[0],Div[1]);
  cSlewing->cd(1);

  tree->Draw("sta2vTOF:sta2vTem>>h1(10, 0, 5, 10, -159, -156)","", "colz");
  Int_t nBin = h1->GetXaxis()->GetNbins();
  TH1F * g1 = new TH1F("g1", "g1", nBin , h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax());

  for(Int_t xBin = 1; xBin < nBin; xBin ++){

    h1->ProjectionY("temp", xBin, xBin);
    Double_t xMax = h1->GetXaxis()->GetBinCenter(xBin);
    Int_t yBin = temp->GetMaximumBin();
    Double_t yMax = temp->GetBinCenter(yBin);
    Double_t yErr = temp->GetStdDev();

    //printf(" (x:y) = (%f,%f)\n", xMax, yMax);

    g1->Fill(xMax, yMax);
    g1->SetBinError(xBin, yErr);

  }
  TF1* fit = new TF1("fit", "[0]/TMath::Power(x-[1],0.5) +[2]",1, 4.5);
  Double_t par[3] = {2, 0, -159};

  fit->SetParameters(par);
  // fit->SetParLimits(1,0,1.5); 
  g1->Fit("fit","R");
  printf("reduced chi-squared = %f \n", fit->GetChisquare()/fit->GetNDF());

  h1->Draw("colz");
  g1->Draw("same");

  cSlewing->cd(2);
  fit->GetParameters(par);
  printf("%f  %f  %f\n",par[0],par[1],par[2]);
  TString plotStr1;
  plotStr1.Form("sta2vTem:(sta2vTOF-%f/TMath::Sqrt(sta2vTem-%f))>>h2(20,-161,-158,20,1,5)",par[0],par[1]);
  tree->Draw(plotStr1,"","colz");

  cSlewing->cd(3);
  // tree->Draw("sta2v:sta2vTOF>>h3(200,-25,0,200,0,3)","","colz");
  tree->Draw("sta2vTOF>>h3(500,-165,-120)");
  h3->Fit("gaus","R","",-159,-156);
  cSlewing->cd(4);
  TString plotStr2;
  plotStr2.Form("(sta2vTOF-%f/TMath::Sqrt(sta2vTem-%f))>>h4(500,-165,-120)",par[0],par[1]);
  tree->Draw(plotStr2);
  h4->Fit("gaus","R","",-161,-158);


}

   /*
  Double_t para[3] = {0.3,20,0.3};
  hist_1h = new TH2F("hist_1h","hist_1h",50,-20,-15,50,0,1);
  tree->Draw("sta1h:sta1hTOF>>hist_1h","CUTG","colz");

  hist_1hpx = new TProfile("hist_1hpx","hist_1hpx", 50,-20,-15);
  hist_1h->ProfileX("hist_1hpx");
  //hist_1h->Draw("same");
  // hist_1hpx->Draw("same");

  tree->Draw("sta1h:sta1hTOF>>hist_1h(20,-19,-16, 20, 0, 1)", "CUTG", "colz");
  Int_t nBin = hist_1h->GetXaxis()->GetNbins();
  TH1F * g1 = new TH1F("g1", "g1", nBin , hist_1h->GetXaxis()->GetXmin(), hist_1h->GetXaxis()->GetXmax());
  for(Int_t xBin = 1; xBin < nBin; xBin ++)
  {
    h1->ProjectionY("temp", xBin, xBin);
    Double_t xMax = h1->GetBinCenter(xBin);
    Int_t yBin = temp->GetMaximumBin();
    Double_t yMax = temp->GetBinCenter(yBin);
    Double_t yErr = temp->GetStdDev();
    printf(" (x:y) = (%f,%f)\n", xMax, yMax);
    g1->Fill(xMax, yMax);
    g1->SetBinError(xBin, yErr);
  }
  h1->Draw("colz");
  g1->Draw("same");

  func1h = new TF1("func1h", fitSlew, -19.2, -17, 3);
  func1h->SetParameters(para);
  //func1h->Draw("same");
  //TF1 * fit = new TF1("fit", "[0]/TMath::Power([1]+x,2) + [2]", -18.5,-17,3)
  hist_1hpx->Fit("func1h", "R");
  printf("reduced chi-squared = %f \n", func1h->GetChisquare()/func1h->GetNDF());

  func1h->GetParameters(para);
  TF1* s1 = new TF1("s1", fitSlew, -19.2, -17, 3); s1->SetParameters(&para[0]); s1->SetLineColor(2); s1->Draw("same");

}

Double_t fitSlew(Double_t *x, Double_t *par){
  Double_t arg = x[0]+par[1];
  return par[0]/(arg*arg)+par[2];
}
   */

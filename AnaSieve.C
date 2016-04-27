TF1 *fit1, *fit2;

void AnaSieve(){
  //gROOT->Reset();
  //gROOT->ProcessLine(".L constant.h");
  // gROOT->ProcessLine(".L nuclei_mass.h");
 //======================================================== InPut setting
  char * rootfile = "run1104.root";
  Int_t Div[2] = {2,2};	 //x,y
  Int_t size[2] = {600,500}; //x,y
  
  Bool_t analysis = 1; // 0 = no analysis, only load root file and gate; 1 = analysis.
  
  Double_t BGscale = 1.0;
  
  //====================================== Load root file
  TFile *f0 = new TFile (rootfile, "read"); 
  TTree *tree = (TTree*)f0->Get("tree");
  printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,	 tree->GetEntries());
  //gStyle->SetOptStat(1112211);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  /*
  //======================================================== Cut/Gate	   
  //------- Graphical Cut
  TCutG *gateP1_a = new TCutG("cutp1_a",5);
   gateP1_a->SetVarX("grTOF1");
   gateP1_a->SetVarY("grdE1");
   gateP1_a->SetPoint(0,200.303,207.99);
   gateP1_a->SetPoint(1,200.411,165.464);
   gateP1_a->SetPoint(2,207.038,152.32);
   gateP1_a->SetPoint(3,207.038,194.072);
   gateP1_a->SetPoint(4,200.303,207.99);

  TCutG *gateP1_b = new TCutG("cutp1_b",5);
   gateP1_b->SetVarX("grTOF1");
   gateP1_b->SetVarY("grdE1");
   gateP1_b->SetPoint(0,200.303+59.364,207.99);
   gateP1_b->SetPoint(1,200.411+59.364,165.464);
   gateP1_b->SetPoint(2,207.038+59.364,152.32);
   gateP1_b->SetPoint(3,207.038+59.364,194.072);
   gateP1_b->SetPoint(4,200.303+59.364,207.99);

  TCutG *gateP2_a = new TCutG("cutp2_a",5);
   gateP2_a->SetVarX("grTOF2");
   gateP2_a->SetVarY("grdE2");
   gateP2_a->SetPoint(0,207.581,201.804);
   gateP2_a->SetPoint(1,207.581,156.959);
   gateP2_a->SetPoint(2,214.533,140.722);
   gateP2_a->SetPoint(3,214.533,188.66);
   gateP2_a->SetPoint(4,207.581,201.804);

  TCutG *gateP2_b = new TCutG("cutp2_b",5);
   gateP2_b->SetVarX("grTOF2");
   gateP2_b->SetVarY("grdE2");
   gateP2_b->SetPoint(0,207.581+59.364,201.804);
   gateP2_b->SetPoint(1,207.581+59.364,156.959);
   gateP2_b->SetPoint(2,214.533+59.364,140.722);
   gateP2_b->SetPoint(3,214.533+59.364,188.66);
   gateP2_b->SetPoint(4,207.581+59.364,201.804);

   TCut gateP = "((cutp1_a || cutp1_b) && (cutp2_a || cutp2_b))";

    printf("........ loaded gates\n"); 
  */
    if( analysis == 0) {
      printf("............. end of Ana.C\n");
      return;
    }   

    printf(".......... start analysis \n");
    //======================================================== Browser or Canvas
    TCanvas * cSieve = new TCanvas("cSieve", "cSieve", 100, 0 , size[0]*Div[0], size[1]*Div[1]);
    cSieve->Divide(Div[0],Div[1]);

    cSieve->cd(1);
    tree->Draw("Atg:grXC>>h1(250,-500,-300,250,-2,2)","","colz");
    h1->GetXaxis()->SetTitle("grXC [mm]");
    h1->GetYaxis()->SetTitle("A_{tg} [deg]");

    cSieve->cd(2);
    tree->Draw("Btg:Atg>>h2(200,-2,2,200,-4.5,4.5)","grXC+420>-10&&grXC+420<10","colz");
    h2->GetXaxis()->SetTitle("A_{tg} [deg]");
    h2->GetYaxis()->SetTitle("B_{tg} [deg]");

    TLine line;
    line.SetLineColor(2);
    line.SetLineWidth(1.8);
    line.DrawLine(-1, -4.5, -1, 4.5);
    line.DrawLine(-0.5, -4.5, -0.5, 4.5);
    line.DrawLine(0, -4.5, 0, 4.5);
    line.DrawLine(0.5, -4.5, 0.5, 4.5);
    line.DrawLine(1, -4.5, 1, 4.5);
    line.DrawLine(-2, -2.4, 2, -2.4);
    line.DrawLine(-2, -1.2, 2, -1.2);
    line.DrawLine(-2, 0, 2, 0);
    line.DrawLine(-2, 1.2, 2, 1.2);
    line.DrawLine(-2, 2.4, 2, 2.4);

    Double_t para1[15]={3000,-1,0.07,4000,-0.5,0.07,5000,0,0.07,4000,0.5,0.07,3000,1,0.07};
    cSieve->cd(3);
    tree->Draw("Atg>>h3(200,-2,2)","grXC+420>-10&&grXC+420<10");
    h3->GetXaxis()->SetTitle("A_{tg} [deg]");

    fit1 = new TF1("fit1",Fit_5Gauss,-2,2,15);
    fit1->SetParameters(para1);
    fit1->SetLineColor(kRed);
    fit1->SetLineStyle(1);
    fit1->SetLineWidth(2);
    h3->Fit("fit1", "R");

    Double_t para2[15]={1000,-2.4,0.4,2000,-1.2,0.4,2500,0,0.4,2000,1.2,0.4,1000,2.4,0.4};
    cSieve->cd(4);
    tree->Draw("Btg>>h4(200,-4.5,4.5)","grXC+420>-10&&grXC+420<10");
    h4->GetXaxis()->SetTitle("B_{tg} [deg]");

    fit2 = new TF1("fit2",Fit_5Gauss,-4.5,4.5,15);
    fit2->SetParameters(para2);
    fit2->SetLineColor(kRed);
    fit2->SetLineStyle(1);
    fit2->SetLineWidth(2);
    h4->Fit("fit2", "R");

  printf("............. end of Ana.C\n");

}


Double_t Fit_5Gauss(Double_t *x, Double_t *para){

  Double_t arg1 = (x[0] - para[1])/para[2];
  Double_t arg2 = (x[0] - para[4])/para[5];
  Double_t arg3 = (x[0] - para[7])/para[8];
  Double_t arg4 = (x[0] - para[10])/para[11];
  Double_t arg5 = (x[0] - para[13])/para[14];
 
  Double_t arg11 = exp(-0.5*TMath::Power(arg1,2));
  Double_t arg21 = exp(-0.5*TMath::Power(arg2,2));
  Double_t arg31 = exp(-0.5*TMath::Power(arg3,2));
  Double_t arg41 = exp(-0.5*TMath::Power(arg4,2));
  Double_t arg51 = exp(-0.5*TMath::Power(arg5,2));

  return para[0] * arg11 + para[3] * arg21 + para[6] * arg31 + para[9] * arg41 + para[12] * arg51;
}


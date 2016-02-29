TF1 * fit;
int NPAR;

void Count3He(Double_t dx = 159.24){
  //gROOT->Reset();
  gROOT->ProcessLine(".!date");

  //======================================================== InPut setting

  Int_t numGauss = 3; //dummy
  
  char * rootfile = "run1034.root";
  Int_t Div[2] = {1,1};  //x,y
  Int_t size[2] = {600,400}; //x,y

  //====================================== Load root file
  NPAR = 3*numGauss;
  TFile *f0 = new TFile (rootfile, "read"); 
  TTree *tree = (TTree*)f0->Get("tree");
  printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
  gStyle->SetOptStat(1112211);
  gStyle->SetOptFit(1);
  //======================================================== Browser or Canvas
  TCanvas * c3HeCount = new TCanvas("c3HeCount", "c3HeCount", 200, 0 , size[0]*Div[0], size[1]*Div[1]);
  c3HeCount->Divide(Div[0],Div[1]);
  c3HeCount->cd(1);
  //======================================================== Cut/Gate

  //------- Graphical Cut
  TCutG * gate3He1_a = new TCutG("cut3He1_a", 5);
  gate3He1_a->SetVarX("grTOF1");
  gate3He1_a->SetVarY("grdE1");
  gate3He1_a->SetPoint(0,186.604,257.051); //#1034
  gate3He1_a->SetPoint(1,186.604,232.051);
  gate3He1_a->SetPoint(2,200.154,211.538);
  gate3He1_a->SetPoint(3,200.154,234.615);
  gate3He1_a->SetPoint(4,186.604,257.051);

  TCutG * gate3He1_b = new TCutG("cut3He1_b", 5);
  gate3He1_b->SetVarX("grTOF1");
  gate3He1_b->SetVarY("grdE1");
  gate3He1_b->SetPoint(0,186.604+98.85,257.051); //#1034
  gate3He1_b->SetPoint(1,186.604+98.85,232.051);
  gate3He1_b->SetPoint(2,200.154+98.85,211.538);
  gate3He1_b->SetPoint(3,200.154+98.85,234.615);
  gate3He1_b->SetPoint(4,186.604+98.85,257.051);
  /*
  TCutG *gate3He2_a = new TCutG("cut3He2_a",8);
  gate3He2_a->SetVarX("grTOF2");
  gate3He2_a->SetVarY("grdE2");
  gate3He2_a->SetPoint(0,195.517,448.418);
  gate3He2_a->SetPoint(1,195.517,412.816);
  gate3He2_a->SetPoint(2,202.037,428.639);
  gate3He2_a->SetPoint(3,207.317,433.386);
  gate3He2_a->SetPoint(4,210.819,420.728);
  gate3He2_a->SetPoint(5,210.873,465.032);
  gate3He2_a->SetPoint(6,207.425,471.361);
  gate3He2_a->SetPoint(7,195.517,448.418);

  TCutG *gate3He2_b = new TCutG("cut3He2_b",8);
  gate3He2_b->SetVarX("grTOF2");
  gate3He2_b->SetVarY("grdE2");
  gate3He2_b->SetPoint(0,195.517+98.85,448.418);
  gate3He2_b->SetPoint(1,195.517+98.85,412.816);
  gate3He2_b->SetPoint(2,202.037+98.85,428.639);
  gate3He2_b->SetPoint(3,207.317+98.85,433.386);
  gate3He2_b->SetPoint(4,210.819+98.85,420.728);
  gate3He2_b->SetPoint(5,210.873+98.85,465.032);
  gate3He2_b->SetPoint(6,207.425+98.85,471.361);
  gate3He2_b->SetPoint(7,195.517+98.85,448.418);
  */

   TCutG *gate3He2_a = new TCutG("cut3He2_a",5);
   gate3He2_a->SetVarX("grTOF2");
   gate3He2_a->SetVarY("grdE2");
   gate3He2_a->SetPoint(0,188.366,467.379);
   gate3He2_a->SetPoint(1,188.366,407.893);
   gate3He2_a->SetPoint(2,201.588,360.07);
   gate3He2_a->SetPoint(3,201.588,402.061);
   gate3He2_a->SetPoint(4,188.366,467.379);

   TCutG *gate3He2_b = new TCutG("cut3He2_b",5);
   gate3He2_b->SetVarX("grTOF2");
   gate3He2_b->SetVarY("grdE2");
   gate3He2_b->SetPoint(0,188.366+98.85,467.379);
   gate3He2_b->SetPoint(1,188.366+98.85,407.893);
   gate3He2_b->SetPoint(2,201.588+98.85,360.07);
   gate3He2_b->SetPoint(3,201.588+98.85,402.061);
   gate3He2_b->SetPoint(4,188.366+98.85,467.379);

  TCut gate3He = "((cut3He1_a || cut3He1_b) && (cut3He2_a || cut3He2_b))";
  TCut gate3He_Time = " grTOF1 ";
  TCut gateEve = "eventID < 3.8e6";
  TCut gatePhi = "-0.4 < grph*TMath::RadToDeg() && grph*TMath::RadToDeg() < 0.4";
  
  //*///======================================================== analysis

  //tree->Draw("grdE1:grTOF1>>h1(500, 100, 350, 500, 0, 500)", "", "colz");
  //tree->Draw("grdE1:grTOF1>>h1g(500, 100, 350, 500, 0, 500)", "cut3He_a || cut3He_b", "colz");
  tree->Draw("grthC*TMath::RadToDeg():grx-10*grthC*TMath::RadToDeg()>>h2(1200,-600,600,600,-2,2)", gate3He+gatePhi, "colz");
  tree->Draw("grx-10*grthC*TMath::RadToDeg()>>h2px(250,-50,200)", gate3He+gatePhi, "colz");

  //Double_t para[9] = {1000, -20, 50, 80, 120, 110, 210, 8, 360}; //3He events
  Double_t para[9] = {1700, 15, 20, 296, 5.4, -0.036, 300, 130, 20}; 
  Double_t para1[9] = {1700, 15, 20, 296, 5.4, -0.036, 300, 130, 20}; 
  //Double_t para[9] = {1700, -20, 50, 5, 110, 1000, 5, -2, 1000}; //neutron events
  Int_t binWidth = h2px->GetBinWidth(1);
  
  //fit = new TF1("fit", "gaus(0)+gaus(3)+gaus(6)", -50, 200);
  fit = new TF1("fit", CusFunc3 , -50, 200, 9);
  fit->SetParameters(para);
  //fit->FixParameter(1, 0);
  //fit->SetParLimits(3,500,800);
  //fit->SetParLimits(4,-10,-1);
  //fit->SetParLimits(7,0,20);
  fit->SetLineColor(kRed);
  fit->SetLineStyle(1);
  fit->SetLineWidth(2);
  h2px->Fit("fit", "R");
  printf("reduced chi-squared = %f \n", fit->GetChisquare()/fit->GetNDF());
  printf("bin Width = %d \n", binWidth);

  fit->GetParameters(para);
  TF1* g1 = new TF1("g1", "gaus(0)", -50, 200); g1->SetParameters(&para1[0]); g1->SetLineColor(kBlue); g1->Draw("same");
  //TF1* g1 = new TF1("g1", "gaus(0)", -600, 600); g1->SetParameters(&para[0]); g1->SetLineColor(4); g1->Draw("same");
  TF1* g2 = new TF1("g2", Poly3Func, -50, 200, 3); g2->SetParameters(&para1[3]); g2->SetLineColor(kGreen); g2->Draw("same");
  TF1* g3 = new TF1("g3", "gaus(0)", -50, 200); g3->SetParameters(&para1[6]); g3->SetLineColor(kViolet); g3->Draw("same");
  TF1* g4 = new TF1("g4", CusFunc3, -50, 200, 9); g4->SetParameters(&para1[0]); g4->SetLineColor(kBlack); g4->Draw("same");

  //Double_t xmin = para[7]-(para[4]-para[7]);
  Double_t xmin = para[1]-3*para[2];
  Double_t xmax = para[4];
        
  Double_t gSide1 = g3->Integral(xmin-dx, xmin);
  Double_t gSide2 = g3->Integral(xmax, xmax + dx);

  //printf("%9s   count (%8.2f,%8.2f): %8.2f \n","g1+g2+g3",xmin,    xmax     , g1->Integral(xmin, xmax) + g2->Integral(xmin, xmax) + g3->Integral(xmin, xmax));
  //printf("%9s   count (%8.2f,%8.2f): %8.2f \n","g1      ",xmin,    xmax     , g1->Integral(xmin, xmax));
  //printf("%9s	count (%8.2f,%8.2f): %8.2f \n","g3 Cent.",xmin,	 xmax	  , g3->Integral(xmin, xmax) );
  //printf("%9s	count (%8.2f,%8.2f): %8.2f \n","g3 Side1",xmin-dx, xmin	  , gSide1);
  //printf("%9s	count (%8.2f,%8.2f): %8.2f \n","g3 Side2",xmax,	 xmax + dx, gSide2);
  //printf("%9s  %8s%8s        count : %8.2f , diff : %8.2f \n","gSide",  "",           "", gSide1+gSide2, g3->Integral(xmin, xmax) - gSide1 - gSide2);

  //TLine line;
  //line.SetLineColor(9);
  //line.DrawLine(xmin, 0, xmin, h2px->GetMaximum());
  //line.DrawLine(xmax, 0, xmax, h2px->GetMaximum());
  //line.DrawLine(xmin-dx, 0, xmin-dx, (h2px->GetMaximum())/2);
  //line.DrawLine(xmax+dx, 0, xmax+dx, (h2px->GetMaximum())/2);

  /*
	TVirtualFitter * fitter = TVirtualFitter::GetFitter();
	assert(fitter !=0); // if abort if fitter == 0
	double * covMatrix = fitter->GetCovarianceMatrix();	
	TF1 * derv_par0 = new TF1("dfdp0",df_dPar, -600, 600, 1);
	TF1 * derv_par1 = new TF1("dfdp1",df_dPar, -600, 600, 1);
	TF1 * derv_par2 = new TF1("dfdp2",df_dPar, -600, 600, 1);
	double c[3];
	double *epar = fit->GetParErrors();
	
	derv_par0->SetParameter(0, 0);
	derv_par1->SetParameter(0, 1);
	derv_par2->SetParameter(0, 2);
	c[0] = derv_par0->Integral(para[1]-3*para[2],para[1]+3*para[2]);	
	c[1] = derv_par1->Integral(para[1]-3*para[2],para[1]+3*para[2]);
	c[2] = derv_par2->Integral(para[1]-3*para[2],para[1]+3*para[2]);
	double sigma_integral = IntegralError(3, c, epar, covMatrix);
	if( TMath::IsNaN(sigma_integral) ) {
      printf("no covMarix\n");
      sigma_integral = IntegralError(3, c, epar);
    }
	//printf("c[0]... %.2f, %.2f, %.2f , :%.2f\n", c[0], c[1], c[2], sigma_integral);

    printf("g1 integration error:%.2f\n", sigma_integral);

    /*
    // Drawing on the plot
    TString textStr;
    TLatex text;
	textStr.Form("mean: %3.1f#pm%3.1f", para[1], fit->GetParError(1));
	text.DrawLatex(0.5, 0.8, textStr.Data());
	textStr.Form("#sigma: %3.1f#pm%3.1f", para[2], fit->GetParError(2));
	text.DrawLatex(0.5, 0.77, textStr.Data());
    Double_t binWidth = h2px->GetBinWidth(1);
	textStr.Form("count: %.0f(%.0f)", g1->Integral(para[1]-3*para[2],para[1]+3*para[2])/binWidth, sigma_integral/binWidth);
	text.DrawLatex(0.5, 0.74, textStr.Data());
    /**/
}


//__________________________________________________________________
Double_t CusFunc(Double_t *x, Double_t *para){
  // A, mean, sigma, x 3 
  Double_t arg1 = (x[0] - para[1])/para[2];
  Double_t arg2 = (x[0] - para[4])/para[5];
  Double_t arg3 = (x[0] - para[7])/para[8];
   
  Double_t arg11 = exp(-0.5*TMath::Power(arg1,2));
  Double_t arg21 = exp(-0.5*TMath::Power(arg2,2));
  Double_t arg31 = exp(-0.5*TMath::Power(arg3,2));
   
  return para[0] / (sqrt(2.*TMath::Pi())*para[2]) * arg11 + para[3]* arg21 + para[6] * arg31;
}

Double_t CusFunc2(Double_t *x, Double_t *para){
  // A, mean, sigma 
  Double_t arg1 = (x[0] - para[1])/para[2];
  Double_t arg11 = exp(-0.5*TMath::Power(arg1,2));
  return para[0] / (sqrt(2.*TMath::Pi())*para[2]) * arg11;
}

Double_t CusFunc3(Double_t *x, Double_t *para){
  // A, mean, sigma *2 + polynomial 
  Double_t arg1 = (x[0] - para[1])/para[2];
  Double_t arg3 = (x[0] - para[7])/para[8];
  Double_t arg11 = exp(-0.5*TMath::Power(arg1,2));
  Double_t arg12 = para[3] + para[4]*x[0] + para[5]*x[0]*x[0];
  Double_t arg13 = exp(-0.5*TMath::Power(arg3,2));
  //return para[0] / (sqrt(2.*TMath::Pi())*para[2]) * arg11 + arg12 + para[6] / (sqrt(2.*TMath::Pi())*para[8]) * arg13;
  return para[0]*arg11 + arg12 + para[6]*arg13;
}

Double_t Poly3Func(Double_t *x, Double_t *para){
return para[0] + para[1]*x[0] + para[2]*x[0]*x[0];
}

//____________________________________________________________________
double IntegralError(int npar, double * c, double * errPar,  double * covMatrix = 0) {   
// calculate the error on the integral given the parameter error and 
// the integrals of the gradient functions c[] 

   double err2 = 0; 
   for (int i = 0; i < npar; ++i) { 
      if (covMatrix == 0) { // assume error are uncorrelated
         err2 += c[i] * c[i] * errPar[i] * errPar[i]; 
      } else {
         double s = 0; 
         for (int j = 0; j < npar; ++j) {
            s += covMatrix[i*npar + j] * c[j]; 
         }
         err2 += c[i] * s; 
      }
   }

   return TMath::Sqrt(err2);
}

//____________________________________________________________________
double df_dPar(double * x, double * p) { 
   // derivative of the function w.r..t parameters
   // use calculated derivatives from TF1::GradientPar

   double* grad = new double [NPAR]; 
   // p is used to specify for which parameter the derivative is computed 
   int ipar = int(p[0] ); 
   //printf("p[0] : %f, ipar: %d \n",p[0], ipar);
   assert (ipar >=0 && ipar < NPAR );

   assert(fit);
   fit->GradientPar(x, grad);

   return grad[ipar]; 
}
/**/

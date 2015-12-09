const Int_t NPAR=9;

void Count3He(){
        gROOT->Reset();
        gROOT->ProcessLine(".!date");

        //======================================================== InPut setting

        char * rootfile = "S_run1035.root";
        Int_t Div[2] = {1,1};  //x,y
        Int_t size[2] = {600,400}; //x,y

	//====================================== Load root file
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
        Double_t BGscale = 1.0;

        TCutG * gate3He_a = new TCutG("cut3He_a", 5);
        gate3He_a->SetVarX("grTOF1");
        gate3He_a->SetVarY("grdE1");
        gate3He_a->SetPoint(0, 183.7, 305.8);
        gate3He_a->SetPoint(1, 188.4, 211.6);
        gate3He_a->SetPoint(2, 214.2, 195.4);
        gate3He_a->SetPoint(3, 225.9, 304.6);
        gate3He_a->SetPoint(4, 208.6, 343.0);
        gate3He_a->SetPoint(5, 183.7, 305.8);

	
        TCutG * gate3He_b = new TCutG("cut3He_b", 5);
        gate3He_b->SetVarX("grTOF1");
        gate3He_b->SetVarY("grdE1");
        gate3He_b->SetPoint(0, 183.7+99, 305.8);
        gate3He_b->SetPoint(1, 188.4+99, 211.6);
        gate3He_b->SetPoint(2, 214.2+99, 195.4);
        gate3He_b->SetPoint(3, 225.9+99, 304.6);
        gate3He_b->SetPoint(4, 208.6+99, 343.0);
        gate3He_b->SetPoint(5, 183.7+99, 305.8);
  
        //*///======================================================== analysis

	//tree->Draw("grdE1:grTOF1>>h1(500, 100, 350, 500, 0, 500)", "", "colz");
	//tree->Draw("grdE1:grTOF1>>h1g(500, 100, 350, 500, 0, 500)", "cut3He_a || cut3He_b", "colz");
	tree->Draw("grth*TMath::RadToDeg():grXC>>h2(600,-1000,1000,600,-1.5,1.5)", "cut3He_a || cut3He_b", "colz");
	tree->Draw("grXC>>h2px(600,-600,600)", "cut3He_a || cut3He_b", "colz");

        //h2->ProjectionX("h2px")->Draw();
        
        //TF1* fit = new TF1("fit", "gaus(0)+gaus(3)+gaus(6)", -600, 600);
        TF1* fit = new TF1("fit", CusFunc, -600, 600, 9);
        Double_t para[9] = {1100, -20, 50, 150, 110, 120, 340, -2, 340};
        fit->SetParameters(para);
        //fit->FixParameter(1, 0);
        fit->SetLineColor(2);
        fit->SetLineStyle(1);
        fit->SetLineWidth(2);
        h2px->Fit("fit", "R");
        printf("reduced chi-squared = %f \n", fit->GetChisquare()/fit->GetNDF());

        fit->GetParameters(para);
        TF1* g1 = new TF1("g1", CusFunc2, -600, 600, 3); g1->SetParameters(&para[0]); g1->SetLineColor(4); g1->Draw("same");
        TF1* g2 = new TF1("g2", "gaus(0)", -600, 600); g2->SetParameters(&para[3]); g2->SetLineColor(3); g2->Draw("same");
        TF1* g3 = new TF1("g3", "gaus(0)", -600, 600); g3->SetParameters(&para[6]); g3->SetLineColor(2); g3->Draw("same");
        
        /*
	TVirtualFitter * fitter = TVirtualFitter::GetFitter();
	assert(fitter !=0);
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
	if( TMath::IsNaN(sigma_integral) ) { printf("no covMarix\n"); sigma_integral = IntegralError(3, c, epar); }
	printf("c[0]... %.2f, %.2f, %.2f , sigma:%.2f\n", c[0], c[1], c[2], sigma_integral);
	
	textStr.Form("mean: %3.1f#pm%3.1f", para[1], fit->GetParError(1));
	text.DrawLatex(0.5, 0.8, textStr.Data());
	textStr.Form("#sigma: %3.1f#pm%3.1f", para[2], fit->GetParError(2));
	text.DrawLatex(0.5, 0.77, textStr.Data());
        Double_t binWidth = h2px->GetBinWidth();
	textStr.Form("count: %.0f(%.0f)", g1->Integral(para[1]-3*para[2],para[1]+3*para[2])/binWidth, sigma_integral/binWidth);
	text.DrawLatex(0.5, 0.74, textStr.Data());
	*/
}

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
   // A, mean, sigma, x 3 
   Double_t arg1 = (x[0] - para[1])/para[2];
   Double_t arg11 = exp(-0.5*TMath::Power(arg1,2));
   return para[0] / (sqrt(2.*TMath::Pi())*para[2]) * arg11;
}





/*
//____________________________________________________________________
double IntegralError(int npar, double * c, double * errPar, 
   double * covMatrix = 0) {   
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

   double grad[NPAR]; 
   // p is used to specify for which parameter the derivative is computed 
   int ipar = int(p[0] ); 
   //printf("p[0] : %f, ipar: %d \n",p[0], ipar);
   assert (ipar >=0 && ipar < NPAR );

   assert(total);
   total->GradientPar(x, grad);

   return grad[ipar]; 
}


TF1 * fit;
const Int_t nPara=6;

void Fit_2Gauss(TH1F* hFit, Double_t mag1, Double_t x1, Double_t sigma1, Double_t mag2, Double_t x2, Double_t sigma2){
   gStyle->SetOptFit(1111);
   //nPara = 6;
   const Double_t minX = hFit->GetXaxis()->GetXmin();
   const Double_t maxX = hFit->GetXaxis()->GetXmax();
   const Double_t binWidth = hFit->GetBinWidth(1);
   fit = new TF1("fit", fitFunction,  minX, maxX, nPara); 
   fit->SetParameters(mag1,  x1, sigma1,
                      mag2,  x2, sigma2);
//   fit->SetParLimits(1, 0,15); // x1
//   fit->SetParLimits(2, 0,50);   // sigma1
//   fit->SetParLimits(4, 20,60); // x2
//   fit->SetParLimits(5, 0,50);   // sigma2
   //fit->FixParameter(1,0);
   //fit->FixParameter(4,0);
   hFit->Fit("fit", "R");
   //Get para
   Double_t para[nPara];
   fit->GetParameters(para);
   
   //seperate function
   TF1 *Gauss1 = new TF1("Gauss1", Gauss, minX, maxX, 3);
   Gauss1->SetLineColor(kBlue);
   TF1 *Gauss2 = new TF1("Gauss2", Gauss, minX, maxX, 3);
   Gauss2->SetLineColor(kGreen);
   
   Gauss1->SetParameters(para);
   Gauss1->Draw("same");
   Double_t para2[3]={para[3],para[4],para[5]}; 
   Gauss2->SetParameters(para2);
   Gauss2->Draw("same");

   TLine line;
   line.SetLineColor(kBlack);
   line.SetLineWidth(2);
   line.DrawLine(para[1]-para[2],0,para[1]-para[2],hFit->GetMaximum());
   line.DrawLine(para[1]+para[2],0,para[1]+para[2],hFit->GetMaximum());

   Double_t cG1 = para[0]/ binWidth;
   Double_t cG2 = para[3]/ binWidth;
   printf("=========  (%8.2f,%8.2f) Count: Gauss1(%f), Gauss2(%f) \n",para[1]-para[2],para[1]+para[2],Gauss1->Integral(para[1]-para[2],para[1]+para[2]), Gauss2->Integral(para[1]-para[2],para[1]+para[2]));
   //printf("========= Count: Gauss1(%f), Gauss2(%f) \n", cG1, cG2);

   //-------------------------Integral error   
   TVirtualFitter * fitter = TVirtualFitter::GetFitter();
   assert(fitter !=0); // if abort if fitter == 0
   double * covMatrix = fitter->GetCovarianceMatrix();
   printf(" covMatrix = \n");
   for( Int_t i = 0; i < 6; i++ ){
    for(Int_t j = 0; j < 6; j++ ){
      printf("%6.2f ", covMatrix[i*6 + j]); 
    }
      printf("\n");
   }
   //get subcovarianMatrix
   double * subcovMatrix = new double[9];
   for( Int_t i = 0; i < 3; i++ ){
    for(Int_t j = 0; j < 3; j++ ){
     subcovMatrix[i*3+j]=covMatrix[(i+3)*6 +3 + j]; //gaus1: i*6; gaus2: (i+3)*6+3
     printf("%6.2f ", covMatrix[(i+3)*6 +3 + j]); 
    }
    printf("\n");
   }	
   TF1 * derv_par0 = new TF1("dfdp0",df_dPar, -600, 600, 1);
   TF1 * derv_par1 = new TF1("dfdp1",df_dPar, -600, 600, 1);
   TF1 * derv_par2 = new TF1("dfdp2",df_dPar, -600, 600, 1);
   double c[3];
   double *epar = fit->GetParErrors();
	
   derv_par0->SetParameter(0, 3);
   derv_par1->SetParameter(0, 4);
   derv_par2->SetParameter(0, 5);
   c[0] = derv_par0->Integral(para[1]-para[2],para[1]+para[2]);	
   c[1] = derv_par1->Integral(para[1]-para[2],para[1]+para[2]);
   c[2] = derv_par2->Integral(para[1]-para[2],para[1]+para[2]);
   double sigma_integral = IntegralError(3, c, epar, subcovMatrix);
   if( TMath::IsNaN(sigma_integral) ) {
     printf("no covMarix\n");
     sigma_integral = IntegralError(3, c, epar);
    }
	//printf("c[0]... %.2f, %.2f, %.2f , :%.2f\n", c[0], c[1], c[2], sigma_integral);

   printf("Gauss2 integration error:%.2f\n", sigma_integral);

}

Double_t fitFunction(Double_t *x, Double_t *para){
   const Int_t nParaGauss = 3;
   const Double_t para2[3] = {para[nParaGauss],para[nParaGauss+1], para[nParaGauss+2]};
   return Gauss(x,para) + Gauss(x,para2);
}

Double_t Gauss(Double_t *x, Double_t *para){
   // magnitude, offset, sigma
   Double_t arg;
   if(para[2]) arg = (x[0] - para[1])/para[2];
   
   //return para[0]*TMath::Gaus(arg,0,1,0);
   //return para[0]*TMath::Gaus(x[0], para[1], para[2], 1);
   return para[0]/TMath::Sqrt(2*TMath::Pi())/para[2]*TMath::Exp(-0.5*TMath::Power(arg,2));
}

//____________________________________________________________________
double IntegralError(int npar, double * c, double * errPar,  double * covMatrix = 0) {   
// calculate the error on the integral given the parameter error and 
// the integrals of the gradient functions c[] 

   double err2 = 0; 
   for (int i = 0; i < npar; ++i) { 
      if (covMatrix == 0) { // assume error are uncorrelated
         err2 += c[i] * c[i] * errPar[i+3] * errPar[i+3]; //for gaus2, i->i+3
      } else {
         double s = 0; 
         for (int j = 0; j < npar; ++j) {
            s += covMatrix[i*npar + j] * c[j]; //covMatrix must be npar x npar
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

   double* grad = new double [nPara]; 
   // p is used to specify for which parameter the derivative is computed 
   int ipar = int(p[0] ); 
   //printf("p[0] : %f, ipar: %d \n",p[0], ipar);
   assert (ipar >=0 && ipar < nPara );

   assert(fit);
   fit->GradientPar(x, grad);

   return grad[ipar]; 
}
/**/

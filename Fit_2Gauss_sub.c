TH1F* Fit_2Gauss_sub(TH1F* hFit, Double_t mag1, Double_t x1, Double_t sigma1, Double_t mag2, Double_t x2, Double_t sigma2){

   printf("...... start to fit %s with 2 Gaussians. \n", hFit->GetName());

   gStyle->SetOptFit(1111);
   const Int_t nPara = 6;
   const Double_t minX = hFit->GetXaxis()->GetXmin();
   const Double_t maxX = hFit->GetXaxis()->GetXmax();
   const Double_t binWidth = hFit->GetBinWidth(1);

   TF1 *fit = new TF1("fit", fitFunction,  minX, maxX, nPara); 
   fit->SetParameters(mag1,  x1, sigma1,
                      mag2,  x2, sigma2);
//   fit->SetParLimits(1, 0,15); // x1
//   fit->SetParLimits(2, 0,50);   // sigma1
//   fit->SetParLimits(4, 20,60); // x2
//   fit->SetParLimits(5, 0,50);   // sigma2
   //fit->FixParameter(1,0);
   //fit->FixParameter(4,0);
   hFit->Fit("fit", "R");
   printf(" reduced chi-squared = %f \n", fit->GetChisquare()/fit->GetNDF());

   //Get para
   Double_t para[nPara];
   fit->GetParameters(para);
   
   //seperate function
   TF1 *Gauss1 = new TF1("Gauss1", Gauss, minX, maxX, 3);
   Gauss1->SetLineColor(kBlue);
   TF1 *Gauss2 = new TF1("Gauss2", Gauss, minX, maxX, 3);
   Gauss2->SetLineColor(7);
   
   Gauss1->SetParameters(para);
   Gauss1->Draw("same");
   Double_t para2[3]={para[3],para[4],para[5]}; 
   Gauss2->SetParameters(para2);
   Gauss2->Draw("same");

   Double_t cG1 = para[0]/ binWidth;
   Double_t cG2 = para[3]/ binWidth;
   //printf("========= Count: Gauss1(%f), Gauss2(%f) \n",Gauss1->Integral(para[1]-3*para[2],para[1]+3*para[2]), Gauss2->Integral(para[4]-3*para[5],para[4]+3*para[5]));
   printf("=====event from Gaussian integration: \n Gauss1 : %f \n Gauss2 : %f \n", cG1, cG2);

   //substract the hFit and save into hOut
   Int_t nBin = hFit->GetNbinsX();
   TH1F * hTemp = new TH1F("hTemp", "hTemp", nBin, minX, maxX);

   for( Int_t i = 1; i <= nBin; i++){
     double x = hFit->GetBinCenter(i); 
     hTemp->Fill(x, Gauss1->Eval(x));
   }

   TH1F* hOut = new TH1F(*hFit - *hTemp);
   TString hOutName; hOutName.Form("%s_sub", hFit->GetName());
   hOut->SetName(hOutName);
   hOut->SetTitle(hOutName);

   delete fit;
   //delete Gauss1;
   //delete Gauss2;
   delete hTemp;

   //hOut->Draw();

   return hOut;
   
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
   return para[0]*TMath::Gaus(x[0], para[1], para[2], 1);
   //return para[0]/TMath::Sqrt(2*TMath::Pi())/para[2]*TMath::Exp(-0.5*TMath::Power(arg,2));
}

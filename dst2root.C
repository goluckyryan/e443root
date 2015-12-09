#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>

using namespace std;

void DataConvertor(TString openFileName){
   //gROOT->Reset();
   gROOT->ProcessLine(".!date");
   gStyle->SetOptStat(0);
   
   TString saveFileName = openFileName;
   saveFileName.Remove(0,5); //remove head
   Int_t len = saveFileName.Length();
   saveFileName.Remove(len  - 4); // remove tail
   saveFileName = saveFileName + ".root";
   printf("%s \n", saveFileName.Data());
   
   
   TFile *f1 = new TFile (saveFileName,"recreate");
   TTree * t1 = new TTree("tree","tree");
   
   //================ make historgram
   /*
   Int_t Div[2] = {3,2};
   TCanvas * cPlot = new TCanvas("cPlot", "Plot", 100,0 , 400*Div[0], 400*Div[1]);
   cPlot->Divide(Div[0],Div[1]);
   
   TH1F * hk1 = new TH1F ("hk1", "k distribution", 30,0, 300);
   TH1F * hk2 = new TH1F ("hk2", "k distribution", 30,0, 300);
   TH1F * hk3 = new TH1F ("hk3", "k distribution", 30,0, 300);
   hk1->SetMinimum(0);
   hk2->SetLineColor(2);
   hk3->SetLineColor(3);
   hk1->SetXTitle("k [MeV/c]");
   hk1->SetYTitle("xsec [a.u.]");
   
   TH1F * htheta_k1 = new TH1F ("htheta_k1", "theta_k distribution", 18, 0, 180);
   TH1F * htheta_k2 = new TH1F ("htheta_k2", "theta_k distribution", 18, 0, 180);
   TH1F * htheta_k3 = new TH1F ("htheta_k3", "theta_k distribution", 18, 0, 180);
   htheta_k1->SetMinimum(0);
   htheta_k2->SetLineColor(2);
   htheta_k3->SetLineColor(3);
   
   TH1F * hangL1 = new TH1F ("hangL1", "Theta1", 100, 70, 90);
   TH1F * hangL2 = new TH1F ("hangL2", "Theta1", 100, 70, 90);
   TH1F * hangL3 = new TH1F ("hangL3", "Theta1", 100, 70, 90);
   hangL2->SetLineColor(2);
   hangL3->SetLineColor(3);
   
   TH1F * hangR1 = new TH1F ("hangR1", "Theta2", 30, 10, 70);  
   TH1F * hangR2 = new TH1F ("hangR2", "Theta2", 30, 10, 70);  
   TH1F * hangR3 = new TH1F ("hangR3", "Theta2", 30, 10, 70);  
   hangR2->SetLineColor(2);
   hangR3->SetLineColor(3);
   
   TH1F * OffPlane1 = new TH1F ("OffPlane1", "Off-Plane Angle", 20, -20, 20);  
   TH1F * OffPlane2 = new TH1F ("OffPlane2", "Off-Plane Angle", 20, -20, 20);  
   TH1F * OffPlane3 = new TH1F ("OffPlane3", "Off-Plane Angle", 20, -20, 20);  
   OffPlane2->SetLineColor(2);
   OffPlane3->SetLineColor(3);
   
   
   */
   //================ Tree branch
   Double_t Tc, theta_c, phi_c, Td, theta_d, phi_d, offPlane;
   Double_t k, theta_k, phi_k, theta_NN;
   Double_t kp, kt;
   Double_t T1, theta1, phi1, T2, theta2, phi2; // the offPlane angle is same in all frame. 
   Double_t xsec1s1, xsec1p3, xsec1p1,xsec1d5, xsec2s1, xsec1d3;
   Double_t asym1s1, asym1p3, asym1p1,asym1d5, asym2s1, asym1d3;
   
   Double_t x1, y1, z1;
   
   t1->Branch("Tc", &Tc, "Tc/D");
   t1->Branch("theta_c", &theta_c, "theta_c/D");
   t1->Branch("phi_d", &phi_c, "phi_c/D");
   
   t1->Branch("Td", &Td, "Td/D");
   t1->Branch("theta_d", &theta_d, "theta_d/D");
   t1->Branch("phi_d", &phi_d, "phi_d/D");
   
   t1->Branch("offPlane", &offPlane, "offPlane/D");
   
   t1->Branch("k", &k, "k/D");
   t1->Branch("theta_k", &theta_k, "theta_k/D");
   t1->Branch("phi_k", &phi_k, "phi_k/D");
   t1->Branch("theta_NN", &theta_NN, "theta_NN/D");
   
   t1->Branch("kp", &kp, "kp/D");
   t1->Branch("kt", &kt, "kt/D");
   
   t1->Branch("T1", &T1, "T1/D");
   t1->Branch("theta1", &theta1, "theta1/D");
   t1->Branch("phi1", &phi1, "phi1/D");
   
   t1->Branch("T2", &T2, "T2/D");
   t1->Branch("theta2", &theta2, "theta2/D");
   t1->Branch("phi2", &phi2, "phi2/D");
   
   t1->Branch("x1", &x1, "x1/D");
   t1->Branch("y1", &y1, "y1/D");
   t1->Branch("z1", &z1, "z1/D");
   
   t1->Branch("xsec1s1", &xsec1s1, "xsec1s1/D");
   t1->Branch("xsec1p3", &xsec1p3, "xsec1p3/D");
   t1->Branch("xsec1p1", &xsec1p1, "xsec1p1/D");
   t1->Branch("xsec1d5", &xsec1d5, "xsec1d5/D");
   t1->Branch("xsec2s1", &xsec2s1, "xsec2s1/D");
   t1->Branch("xsec1d3", &xsec1d3, "xsec1d3/D");
   
   t1->Branch("asym1s1", &asym1s1, "asym1s1/D");
   t1->Branch("asym1p3", &asym1p3, "asym1p3/D");
   t1->Branch("asym1p1", &asym1p1, "asym1p1/D");
   t1->Branch("asym1d5", &asym1d5, "asym1d5/D");
   t1->Branch("asym2s1", &asym2s1, "asym2s1/D");
   t1->Branch("asym1d3", &asym1d3, "asym1d3/D");
   
   
   //=================== read file
   ifstream fp;
   fp.open(openFileName);
   
   const double deg2rad = 3.141592654/180;
   
   Int_t lineLength = 0;
   Int_t lineLength_valid = 0;
   Int_t lineNum = 0;
   Int_t DataLineNum = 0;
   string line;
   
   double DTc;
   double Dthetac;
   double Dthetad;
   double Dphic ;
   double Dphid ;
   double totalXsec_1s1 = 0;
   double totalXsec_1p3 = 0;
   double totalXsec_1p1 = 0;
   double totalXsec_1d5 = 0;
   double totalXsec_2s1 = 0;
   double totalXsec_1d3 = 0;
   
   do{
      
      lineNum ++;
      getline(fp, line); // read fp and store in line
      
      lineLength = line.length();
      
      string lineHEAD = line.substr(0,2);
      
      if ( lineHEAD == "##" ) printf("%s\n",line.c_str());
      if ( lineHEAD == "#X" ) DTc = atoi(line.substr(17,2).c_str());
      if ( lineHEAD == "#Y" ) Dthetac = atoi(line.substr(17,2).c_str())*deg2rad;
      if ( lineHEAD == "#Z" ) Dthetad = atoi(line.substr(17,2).c_str())*deg2rad;
      if ( lineHEAD == "#B" ) Dphic = atoi(line.substr(17,2).c_str())*deg2rad;
      if ( lineHEAD == "#A" ) Dphid = atoi(line.substr(17,2).c_str())*deg2rad;
      
      
      if( lineHEAD == "# ") lineLength_valid = lineLength;
      //printf("lineNum:%d, linelength = %d, valid_length = %d, lineHEAD : %s \n", lineNum, lineLength, lineLength_valid, lineHEAD);
      
      if( /*lineHEAD == "  " && */ lineNum > 10 && lineLength == lineLength_valid){
         DataLineNum ++;
         Tc          = atof(line.substr(12*(1-1)+1,12).c_str());
         theta_c     = atof(line.substr(12*(2-1)+1,12).c_str());
         phi_c       = atof(line.substr(12*(3-1)+1,12).c_str());
         Td          = atof(line.substr(12*(4-1)+1,12).c_str());
         theta_d = abs(atof(line.substr(12*(5-1)+1,12).c_str()));
         phi_d   = abs(atof(line.substr(12*(6-1)+1,12).c_str()));
         offPlane    = atof(line.substr(12*(7-1)+1,12).c_str());
         
         k        = atof(line.substr(12*(8-1)+1 ,12).c_str());
         theta_k  = atof(line.substr(12*(9-1)+1 ,12).c_str());
         phi_k    = atof(line.substr(12*(10-1)+1,12).c_str());
         theta_NN = atof(line.substr(12*(11-1)+1,12).c_str());
         
         kt = k*sin(theta_NN*deg2rad);
         kp = k*cos(theta_NN*deg2rad);
         
         T1     = atof(line.substr(12*(13-1)+1,12).c_str());
         theta1 = atof(line.substr(12*(14-1)+1,12).c_str());
         phi1   = atof(line.substr(12*(15-1)+1,12).c_str());
         
         Double_t dis = 1200/cos((60-theta1)*deg2rad);
         
         x1 = dis*sin(theta1*deg2rad)*cos(phi1*deg2rad);
         y1 = dis*sin(theta1*deg2rad)*sin(phi1*deg2rad);
         z1 = dis*cos(theta1*deg2rad);
         
         T2     = atof(line.substr(12*(16-1)+1,12).c_str());
         theta2 = atof(line.substr(12*(17-1)+1,12).c_str());
         phi2   = atof(line.substr(12*(18-1)+1,12).c_str());
         
         
         //----------- Acceptance gate
         //if ( TMath::Abs(theta1 -45) < 15 && TMath::Abs(theta2 - 45) < 15 ){ 
         
         xsec1s1 = atof(line.substr(12*(19-1)+1,12+1).c_str());
         xsec1p3 = atof(line.substr(12*(21-1)+1,12+1).c_str());
         xsec1p1 = atof(line.substr(12*(23-1)+1,12+1).c_str());
         xsec1d5 = atof(line.substr(12*(25-1)+1,12+1).c_str());
         xsec2s1 = atof(line.substr(12*(27-1)+1,12+1).c_str());
         xsec1d3 = atof(line.substr(12*(29-1)+1,12+1).c_str());
                  
         asym1s1 = atof(line.substr(12*(20-1)+1,12+1).c_str());
         asym1p3 = atof(line.substr(12*(22-1)+1,12+1).c_str());
         asym1p1 = atof(line.substr(12*(24-1)+1,12+1).c_str());
         asym1d5 = atof(line.substr(12*(26-1)+1,12+1).c_str());
         asym2s1 = atof(line.substr(12*(28-1)+1,12+1).c_str());
         asym1d3 = atof(line.substr(12*(30-1)+1,12+1).c_str());
         
         double ang = 0;
         
         ang = sin(theta_c*deg2rad)*sin(theta_d*deg2rad)*Dthetac*Dthetad*DTc*Dphic*Dphid;
         
         //printf("%4d(%3d), T_c = %7.1f, theta_c = %5.1f(%6.2f), theta_d = %5.1f(%6.2f), phi_c = %5.1f, beta_d = %5.1f , Xsec_1s1 = %12.9f, fac = %12.9f, %12.9f\n"
         //     , lineNum, lineLength, Tc, thetac/deg2rad, theta1, thetad/deg2rad, theta2, phi1, beta_d, Xsec_1s1, ang, ang*Xsec_1s1);
            
         totalXsec_1s1 += ang*xsec1s1;
         totalXsec_1p3 += ang*xsec1p3;
         totalXsec_1p1 += ang*xsec1p1;
         totalXsec_1d5 += ang*xsec1d5;
         totalXsec_2s1 += ang*xsec2s1;
         totalXsec_1d3 += ang*xsec1d3;
         
         //---------- Fill histogram
         /*
         hk1->Fill(k, xsec1s1);
         hk2->Fill(k, xsec1p1);
         hk3->Fill(k, xsec1d5);
         
         htheta_k1->Fill(TMath::Abs(theta_k), xsec1s1);
         htheta_k2->Fill(TMath::Abs(theta_k), xsec1p1);
         htheta_k3->Fill(TMath::Abs(theta_k), xsec1d5);
         
         hangL1->Fill(theta1+theta2, xsec1s1);
         hangL2->Fill(theta1+theta2, xsec1p1);
         hangL3->Fill(theta1+theta2, xsec1d5);
         
         hangR1->Fill(theta2, xsec1s1);
         hangR2->Fill(theta2, xsec1p1);
         hangR3->Fill(theta2, xsec1d5);
         
         OffPlane1->Fill(offPlane, xsec1s1);
         OffPlane2->Fill(offPlane, xsec1p1);
         OffPlane3->Fill(offPlane, xsec1d5);
         */
         //printf("k:%.2f, xsec1s1:%.2f \n", k, xsec1s1);
         //----------- Fill        
         f1->cd(); //set focus on this file
         t1->Fill(); 
         //}
      }
   
   }while(! fp.eof());
   
   
  //  }while( lineLength != 34);
  printf("=======================================\n");
  printf(" DTc %2.0f, Dangc %2.0f, Dangd %2.0f, Dphic %2.0f, Dphid %2.0f \n", DTc, Dthetac/deg2rad, Dthetad/deg2rad, Dphic/deg2rad, Dphid/deg2rad);
  printf(" total number of line = %d\n Total X-sec: \n", DataLineNum-1);
  printf("1s1 = %14.6f ub\n", totalXsec_1s1);
  printf("1p3 = %14.6f ub\n", totalXsec_1p3);
  printf("1p1 = %14.6f ub\n", totalXsec_1p1);
  printf("1d5 = %14.6f ub\n", totalXsec_1d5);
  printf("2s1 = %14.6f ub\n", totalXsec_2s1);
  printf("1d3 = %14.6f ub\n", totalXsec_1d3);
   
   f1->cd();
   t1->Write();
   fp.close();
   //f1->Close();
   
   
   
   
   /*
   cPlot->cd(1);
   hangL->Draw();
   hangL2->Draw("same");
   hangL3->Draw("same");
   
   cPlot->cd(2);
   hangR3->Draw("");
   hangR2->Draw("same");
   hangR->Draw("same");
   
   cPlot->cd(3);
   hk->Draw();
   hk2->Draw("same");
   hk3->Draw("same");
   
   cPlot->cd(4);
   htheta_k->Draw();
   htheta_k2->Draw("same");
   htheta_k3->Draw("same");
   
   cPlot->cd(5);
   OffPlane1->Draw("same");
   OffPlane2->Draw("same");
   OffPlane3->Draw("same");
   
   cPlot->cd(3);
   TLatex text;
   text.SetNDC();
   text.SetTextColor(4);   text.DrawLatex(0.7, 0.7, "1s1/2");
   text.SetTextColor(2);   text.DrawLatex(0.7, 0.6, "1p1/2");
   text.SetTextColor(3);   text.DrawLatex(0.7, 0.5, "1d5/2");
   
   */
   /*//======================================================================
   TFile *f1 = new TFile ("test.root","read");
   TTree *t1 = (TTree*)f1->Get("t1");
   
   */
   
   printf("............ done!\n");
}

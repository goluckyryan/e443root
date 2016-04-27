#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <cstring>
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "constant.h"


using namespace std;

int main(int argc, char * argv[]){

  Int_t RunName = 0;
  Int_t nEntries = 0;

  if( argc > 3 || argc == 1){
    printf("Usage: ./dst2root RunName [opt:#event]\n");
    exit(-1);
  }
  
  if( argc == 3 ){
    nEntries = atoi(argv[2]);
    RunName = atoi(argv[1]);
  }else if(argc == 2){
    nEntries = 9999999;
    RunName = atoi(argv[1]);
  }
  
  printf("=============================\n");
  gROOT->ProcessLine(".!date");

  TString openFileName;
  openFileName.Form("./dstroot/run%04d_asci.dst", RunName);
  printf("input <====== %s \n", openFileName.Data());
  
  TString saveFileName;
  saveFileName.Form("run%04d.root", RunName);
  printf("output =====> %s , nEntries:%d\n", saveFileName.Data(), nEntries);
   
  TFile *f1 = new TFile (saveFileName,"recreate");
  TTree * t1 = new TTree("tree","tree");

  //================ Tree branch
  //---------eventID
  Int_t eventID = 0;
  //---------coinReg
  Int_t coinReg = -1;
  //----------GR
  Double_t grx,gry, grth, grph,dummy;

  Double_t grXC,grthC,grphC,Atg,Btg;


  //------------make tree branch
  t1->Branch("eventID", &eventID, "eventID/I");
  t1->Branch("coinReg", &coinReg, "coinReg/I");

  t1->Branch("grx", &grx, "grx/D");
  t1->Branch("gry", &gry, "gry/D");
  t1->Branch("grth", &grth, "grth/D");
  t1->Branch("grph", &grph, "grph/D");

  t1->Branch("grXC", &grXC, "grXC/D");
  t1->Branch("grthC", &grthC, "grthC/D");
  t1->Branch("grphC", &grphC, "grphC/D");
  t1->Branch("Atg", &Atg, "Atg/D");
  t1->Branch("Btg", &Btg, "Btg/D");
  
  //=================== read file
  ifstream fp1;
  fp1.open(openFileName);

  if( !fp1.is_open() ) {
    printf("******* cannot open dst file : %s\n", openFileName.Data());
    exit(-1);
  }

  TBenchmark clock;
  Bool_t shown = 0;
   
  clock.Reset();
  clock.Start("timer");
  shown = 0;
  
  do{
    
    eventID ++;
    if( eventID > nEntries) break;

   //================================initialization custom variables
   //---------coinReg
   coinReg = -1;
   //--------- GR
   grXC = TMath::QuietNaN();
   grthC = TMath::QuietNaN();
   grphC = TMath::QuietNaN();
   Atg = TMath::QuietNaN();
   Btg = TMath::QuietNaN();
   //grXAux = TMath::QuietNaN();

   //============================================= read from file
    //------------ CoinReg
    fp1>>coinReg;
    //----------- GR
    fp1>>grx;   if(grx == -1)  grx  = TMath::QuietNaN();
    fp1>>grth;  if(grth == -1) grth = TMath::QuietNaN();
    fp1>>gry;   if(gry == -1)  gry  = TMath::QuietNaN();
    fp1>>grph;  if(grph == -1) grph = TMath::QuietNaN();    
    fp1>>dummy;

    //===================================== 2ndary data processing
    //------------- sieve slit calibration
    //grXC = grx -20 ; //same as graf_conv 
    //grthC  = -0.4191*(grth-8.770e-3-3.591e-5*grx); //same as graf_conv, should be equal to incident theta
    grthC  = -0.4191*(grth-8.770e-3-0.52e-4*grx);//rad
    //grthC = -(grth*TMath::RadToDeg()-0.003820*grx); //degree ;
    //Atg = -0.025197+57.36396371*grthC+168.031564*TMath::Power(grthC,2);// rescale to degree
    // Atg = -0.025197+1.00119*grthC*TMath::RadToDeg()+0.051185*TMath::Power(grthC*TMath::RadToDeg(),2);// rescale to degree
    Atg = grthC*TMath::RadToDeg();

    //grXC = grx+301.7-2.302*grthC+1.784*TMath::Power(grthC,2);
    //grXC = -0.26534-111.149*grthC+1.002405*grx+0.798107*grthC*grx+32920.01*TMath::Power(grthC,2)-8.33597*TMath::Power(grthC,2)*grx;
    grXC = -75.0936*grthC+1.001007*grx+0.875889*grthC*grx+31954.5*TMath::Power(grthC,2)-95.386*grph-0.24261*grph*grx-34943.5*TMath::Power(grph,2);
 
    grphC = 5.918076442*grph+1.6082e-5*grx-0.002338054*grph*grXC;
    //Btg = 339.573*grph+0.001013*grx+5.92865*grthC-0.11561*grph*grXC+2515.476*grph*grthC;
    //Btg = 5.926667*grph*TMath::RadToDeg()+0.001013*grx+0.103474*grthC*TMath::RadToDeg()-0.00202*grph*TMath::RadToDeg()*grXC+0.766258*grph*TMath::RadToDeg()*grthC*TMath::RadToDeg();
    // Btg = 5.957953889*grph*TMath::RadToDeg()+0.001046327*grx+0.157376054*grthC*TMath::RadToDeg()-0.001746444*grph*TMath::RadToDeg()*grXC+0.82117321*grph*TMath::RadToDeg()*grthC*TMath::RadToDeg()+4.66522e-6*TMath::Power(grx,2)-2.07783e-11*TMath::Power(grx,4);
    Btg = 5.986529587*grph*TMath::RadToDeg()+0.002085198*grx+0.154225791*grthC*TMath::RadToDeg()-0.001790694*grph*TMath::RadToDeg()*grXC+0.805335911*grph*TMath::RadToDeg()*grthC*TMath::RadToDeg()+3.21763e-6*TMath::Power(grXC,2)-6.58902e-9*TMath::Power(grXC,3)-1.19045e-11*TMath::Power(grXC,4);
    // grYC1 = gry-52+0.093*(grXC-116.5)+0.0092*(grXC-116.5)*(gry-52);
    //grYC2 = grYC1+33.861*(Atg+1)*grph*TMath::RadToDeg();
    //Btg = grph*TMath::RadToDeg()-0.02299*grYC2-0.02299;
    //grXC = grx -380./2.3*grthC*TMath::RadToDeg();//#1035 corrected for (d,3He)
    //grXC = grx -10.*grthC*TMath::RadToDeg();//#1034 corrected for 12C(d,3He)
    //grXAux = grXC -380./2.3*grthC*TMath::RadToDeg();
   
    //----------- Fill       
    f1->cd(); //set focus on this file
    t1->Fill(); 
      
    ///______________________________________________________________   
    clock.Stop("timer");
    Double_t time = clock.GetRealTime("timer");
    clock.Start("timer");

    if ( !shown ) {
      if (fmod(time, 10) < 1 ){
        printf( "%10d[%5.2f%%]|%6.1f| %3.0f min %2.0f sec\n", 
                eventID,
                eventID*100./nEntries,
                time,
                TMath::Floor(time/60.), time - TMath::Floor(time/60.)*60.);
        shown = 1;
      }
    }else{
      if (fmod(time, 10) > 9 ){
        shown = 0;
      }
    }

  }while(! fp1.eof());

  f1->cd(); //set focus on this file
  f1->Write(); 
  f1->Close();
  
  printf("=======================================\n");
  printf("total number of event: %d \n", eventID);
  printf("............ done!\n");
}

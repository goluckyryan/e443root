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

void dst2root(Int_t RunName, Int_t nEntries = 990000000){ //the file name should be "XXXX"
  printf("=============================\n");
  gROOT->ProcessLine(".!date");

  TString openFileName;
  openFileName.Form("../dstroot/run%04d_asci.dst", RunName);
  printf("input <====== %s \n", openFileName.Data());
  
  TString saveFileName;
  saveFileName.Form("run%04d.root", RunName);
  printf("output =====> %s , nEntries:%d\n", saveFileName.Data(), nEntries);
   
  TFile *f1 = new TFile (saveFileName,"recreate");
  TTree * t1 = new TTree("tree","tree");
  TRandom * rand = new TRandom();

  //================ Primary data
  Double_t gradc[4],grtdc[4],grrf,adc[12],adc2[14],tdc[12],tdc2[12],lrf[3],dummy;
  

  //================ Tree branch
  //---------eventID
  Int_t eventID = 0;
  //---------coinReg
  Int_t coinReg = -1;
  //--------- gate
  Int_t vetogate = 0;
  //----------GR
  Double_t grx,gry, grth, grph;
  //--------- GR plastic
  Double_t grdE1, grdE2;
  Double_t grT1avg, grT2avg;
  //--------- BAND telesope
  Double_t badEl,badEr,blo1,blo2,blo3,blo4;
  Double_t blo1Tavg,blo2Tavg,blo3Tavg,blo4Tavg;
  //--------- BAND Stack
  Double_t sta1h,sta2h,sta1v,sta2v,sta3v,sta4v;
  Double_t sta1hTavg,sta2hTavg,sta1vTavg,sta2vTavg,sta3vTavg,sta4vTavg;
  Double_t sta1hTdif,sta2hTdif,sta1vTdif,sta2vTdif,sta3vTdif,sta4vTdif;
  Int_t staMult;
  //--------- BAND liquid
  Double_t liqlf,liqld,liqrf,liqrd;
  //--------- GRRF & BANDRF
  Double_t grf,brf;
  //-------------------------------------------------Phsyics
  Double_t grTOF1, grTOF2;
  Double_t badElTOF,badErTOF,blo1TOF,blo2TOF,blo3TOF,blo4TOF;
  Double_t sta1hTOF,sta2hTOF,sta1vTOF,sta2vTOF,sta3vTOF,sta4vTOF;
  Double_t liqlTOF,liqrTOF;
  Double_t sta_odd,sta_even,sta_ratio,sta_sum;

  Double_t grXC,grthC;
  Double_t grXAux;


  //------------make tree branch
  t1->Branch("eventID", &eventID, "eventID/I");
  t1->Branch("coinReg", &coinReg, "coinReg/I");
  t1->Branch("vetogate", &vetogate, "vetogate/I");
   
  t1->Branch("grx", &grx, "grx/D");
  t1->Branch("gry", &gry, "gry/D");
  t1->Branch("grth", &grth, "grth/D");
  t1->Branch("grph", &grph, "grph/D");

  t1->Branch("grdE1", &grdE1, "grdE1/D");
  t1->Branch("grdE2", &grdE2, "grdE2/D");
  t1->Branch("grTOF1", &grTOF1, "grTOF1/D");
  t1->Branch("grTOF2", &grTOF2, "grTOF2/D");
   
  t1->Branch("grf", &grf, "grf/D");
  t1->Branch("grXC", &grXC, "grXC/D");
  t1->Branch("grthC", &grthC, "grthC/D");
  t1->Branch("grXAux", &grXAux, "grXAux/D");

  t1->Branch("badEl", &badEl, "badEl/D");
  t1->Branch("badEr", &badEr, "badEr/D");

  
  t1->Branch("blo1", &blo1, "blo1/D");
  t1->Branch("blo2", &blo2, "blo2/D");
  t1->Branch("blo3", &blo3, "blo3/D");
  t1->Branch("blo4", &blo4, "blo4/D");
  t1->Branch("badElTOF", &badElTOF, "badElTOF/D");
  t1->Branch("badErTOF", &badErTOF, "badErTOF/D");
  t1->Branch("blo1Tavg", &blo1Tavg, "blo1Tavg/D");
  t1->Branch("blo2Tavg", &blo2Tavg, "blo2Tavg/D");
  t1->Branch("blo3Tavg", &blo3Tavg, "blo3Tavg/D");
  t1->Branch("blo4Tavg", &blo4Tavg, "blo4Tavg/D");
  t1->Branch("blo1TOF", &blo1TOF, "blo1TOF/D");
  t1->Branch("blo2TOF", &blo2TOF, "blo2TOF/D");
  t1->Branch("blo3TOF", &blo3TOF, "blo3TOF/D");
  t1->Branch("blo4TOF", &blo4TOF, "blo4TOF/D");

  t1->Branch("adc", adc, "adc[12]/D");
  t1->Branch("tdc", tdc, "tdc[12]/D");
  t1->Branch("staMult", &staMult, "staMult/I");
  t1->Branch("sta1h", &sta1h, "sta1h/D");
  t1->Branch("sta2h", &sta2h, "sta2h/D");
  t1->Branch("sta1v", &sta1v, "sta1v/D");
  t1->Branch("sta2v", &sta2v, "sta2v/D");
  t1->Branch("sta3v", &sta3v, "sta3v/D");
  t1->Branch("sta4v", &sta4v, "sta4v/D");
  
  t1->Branch("sta_odd", &sta_odd, "sta_odd/D");
  t1->Branch("sta_even", &sta_even, "sta_even/D");
  t1->Branch("sta_sum", &sta_sum, "sta_sum/D");
  t1->Branch("sta_ratio", &sta_ratio, "sta_ratio/D");
  t1->Branch("sta1hTavg", &sta1hTavg, "sta1hTavg/D");
  t1->Branch("sta2hTavg", &sta2hTavg, "sta2hTavg/D");
  t1->Branch("sta1vTavg", &sta1vTavg, "sta1vTavg/D");
  t1->Branch("sta2vTavg", &sta2vTavg, "sta2vTavg/D");
  t1->Branch("sta3vTavg", &sta3vTavg, "sta3vTavg/D");
  t1->Branch("sta4vTavg", &sta4vTavg, "sta4vTavg/D");

  t1->Branch("sta1hTdif", &sta1hTdif, "sta1hTdif/D");
  t1->Branch("sta2hTdif", &sta2hTdif, "sta2hTdif/D");
  t1->Branch("sta1vTdif", &sta1vTdif, "sta1vTdif/D");
  t1->Branch("sta2vTdif", &sta2vTdif, "sta2vTdif/D");
  t1->Branch("sta3vTdif", &sta3vTdif, "sta3vTdif/D");
  t1->Branch("sta4vTdif", &sta4vTdif, "sta4vTdif/D");
  
  t1->Branch("sta1hTOF", &sta1hTOF, "sta1hTOF/D");
  t1->Branch("sta2hTOF", &sta2hTOF, "sta2hTOF/D");
  t1->Branch("sta1vTOF", &sta1vTOF, "sta1vTOF/D");
  t1->Branch("sta2vTOF", &sta2vTOF, "sta2vTOF/D");
  t1->Branch("sta3vTOF", &sta3vTOF, "sta3vTOF/D");
  t1->Branch("sta4vTOF", &sta4vTOF, "sta4vTOF/D");

  t1->Branch("liqlf", &liqlf, "liqlf/D");
  t1->Branch("liqld", &liqld, "liqld/D");
  t1->Branch("liqrf", &liqrf, "liqrf/D");
  t1->Branch("liqrd", &liqrd, "liqrd/D");
  t1->Branch("liqlTOF", &liqlTOF, "liqlTOF/D");
  t1->Branch("liqrTOF", &liqrTOF, "liqrTOF/D");

  t1->Branch("brf", &brf, "brf/D");
   
   
  //=================== read file
  ifstream fp;
  fp.open(openFileName);

  if( !fp.is_open() ) {
    printf("******* cannot open dst file\n");
    return;
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
   //--------- gate
   vetogate = -1;
   //--------- GR
   grdE1 = TMath::QuietNaN();
   grdE2 = TMath::QuietNaN();
   grT1avg = TMath::QuietNaN();
   grT2avg = TMath::QuietNaN();
   grTOF1 = TMath::QuietNaN();
   grTOF2 = TMath::QuietNaN();

   grXC = TMath::QuietNaN();
   grthC = TMath::QuietNaN();
   grXAux = TMath::QuietNaN();

   //--------- BAND_TELE
   badEl = TMath::QuietNaN();
   badEr = TMath::QuietNaN();
   blo1 = TMath::QuietNaN();
   blo2 = TMath::QuietNaN();
   blo3 = TMath::QuietNaN();
   blo4 = TMath::QuietNaN();
   blo1Tavg = TMath::QuietNaN();
   blo2Tavg = TMath::QuietNaN();
   blo3Tavg = TMath::QuietNaN();
   blo4Tavg = TMath::QuietNaN();
   badElTOF = TMath::QuietNaN();
   badErTOF = TMath::QuietNaN();
   blo1TOF = TMath::QuietNaN();
   blo2TOF = TMath::QuietNaN();
   blo3TOF = TMath::QuietNaN();
   blo4TOF = TMath::QuietNaN();
   //--------- BAND_STACK
   sta1h = TMath::QuietNaN();
   sta2h = TMath::QuietNaN();
   sta1v = TMath::QuietNaN();
   sta2v = TMath::QuietNaN();
   sta3v = TMath::QuietNaN();
   sta4v = TMath::QuietNaN();
   sta1hTavg = TMath::QuietNaN();
   sta2hTavg = TMath::QuietNaN();
   sta1vTavg = TMath::QuietNaN();
   sta2vTavg = TMath::QuietNaN();
   sta3vTavg = TMath::QuietNaN();
   sta4vTavg = TMath::QuietNaN();
   
   sta1hTdif = TMath::QuietNaN();
   sta2hTdif = TMath::QuietNaN();
   sta1vTdif = TMath::QuietNaN();
   sta2vTdif = TMath::QuietNaN();
   sta3vTdif = TMath::QuietNaN();
   sta4vTdif = TMath::QuietNaN();
     
   sta1hTOF = TMath::QuietNaN();
   sta2hTOF = TMath::QuietNaN();
   sta1vTOF = TMath::QuietNaN();
   sta2vTOF = TMath::QuietNaN();
   sta3vTOF = TMath::QuietNaN();
   sta4vTOF = TMath::QuietNaN();
   sta_odd = TMath::QuietNaN();
   sta_even = TMath::QuietNaN();
   sta_ratio = TMath::QuietNaN();
   sta_sum = TMath::QuietNaN();
   //--------- BAND_LIQUID
   liqlf = TMath::QuietNaN();
   liqld = TMath::QuietNaN();
   liqrf = TMath::QuietNaN();
   liqrd = TMath::QuietNaN();
   //--------- RF
   grf = TMath::QuietNaN();
   brf = TMath::QuietNaN();

   //============================================= read from file
    //------------ CoinReg
    fp>>coinReg;
    //----------- GR
    fp>>grx;   if(grx == -1)  grx  = TMath::QuietNaN();
    fp>>grth;  if(grth == -1) grth = TMath::QuietNaN();
    fp>>gry;   if(gry == -1)  gry  = TMath::QuietNaN();
    fp>>grph;  if(grph == -1) grph = TMath::QuietNaN();    
    //----------- GR ADC, TDC
    fp>>gradc[0]; 
    fp>>gradc[1];
    fp>>gradc[2];
    fp>>gradc[3];
    fp>>grtdc[0]; if(grtdc[0] == -1) grtdc[0] = TMath::QuietNaN();
    fp>>grtdc[1]; if(grtdc[1] == -1) grtdc[1] = TMath::QuietNaN();  
    fp>>grtdc[2]; if(grtdc[2] == -1) grtdc[2] = TMath::QuietNaN();  
    fp>>grtdc[3]; if(grtdc[3] == -1) grtdc[3] = TMath::QuietNaN();  
    //------------ GR rf
    fp>>grrf;
    //------------ Stack ADC
    fp>>adc[0];  //if(adc[0]  < 5) adc[0] = TMath::QuietNaN();
    fp>>adc[1];  //if(adc[1]  < 5) adc[1] = TMath::QuietNaN();
    fp>>adc[2];  //if(adc[2]  < 5) adc[2] = TMath::QuietNaN();
    fp>>adc[3];  //if(adc[3]  < 5) adc[3] = TMath::QuietNaN();
    fp>>adc[4];  //if(adc[4]  < 5) adc[4] = TMath::QuietNaN();
    fp>>adc[5];  //if(adc[5]  < 5) adc[5] = TMath::QuietNaN();
    fp>>adc[6];  //if(adc[6]  < 5) adc[6] = TMath::QuietNaN();
    fp>>adc[7];  //if(adc[7]  < 5) adc[7] = TMath::QuietNaN();
    fp>>adc[8];  //if(adc[8]  < 5) adc[8] = TMath::QuietNaN();
    fp>>adc[9];  //if(adc[9]  < 5) adc[9] = TMath::QuietNaN();
    fp>>adc[10]; //if(adc[10] < 5) adc[10] = TMath::QuietNaN();
    fp>>adc[11]; //if(adc[11] < 5) adc[11] = TMath::QuietNaN();
    //------------ BAND telesope ADC, Liquid ADC
    fp>>adc2[0];
    fp>>adc2[1];
    fp>>adc2[2];
    fp>>adc2[3];
    fp>>adc2[4];
    fp>>adc2[5];
    fp>>adc2[6];
    fp>>adc2[7];
    fp>>adc2[8];
    fp>>adc2[9];
    fp>>adc2[10];
    fp>>adc2[11];
    fp>>adc2[12];
    fp>>adc2[13];
    //------------ STACK TDC
    fp>>tdc[0];  if(tdc[0]  == -1) tdc[0] = TMath::QuietNaN();
    fp>>tdc[1];  if(tdc[1]  == -1) tdc[1] = TMath::QuietNaN();
    fp>>tdc[2];  if(tdc[2]  == -1) tdc[2] = TMath::QuietNaN();
    fp>>tdc[3];  if(tdc[3]  == -1) tdc[3] = TMath::QuietNaN();
    fp>>tdc[4];  if(tdc[4]  == -1) tdc[4] = TMath::QuietNaN();
    fp>>tdc[5];  if(tdc[5]  == -1) tdc[5] = TMath::QuietNaN();
    fp>>tdc[6];  if(tdc[6]  == -1) tdc[6] = TMath::QuietNaN();
    fp>>tdc[7];  if(tdc[7]  == -1) tdc[7] = TMath::QuietNaN();
    fp>>tdc[8];  if(tdc[8]  == -1) tdc[8] = TMath::QuietNaN();
    fp>>tdc[9];  if(tdc[9]  == -1) tdc[9] = TMath::QuietNaN();
    fp>>tdc[10]; if(tdc[10] == -1) tdc[10] = TMath::QuietNaN();
    fp>>tdc[11]; if(tdc[11] == -1) tdc[11] = TMath::QuietNaN();
    //------------- BAND telesope, Liquid TDC
    fp>>tdc2[0]; if(tdc2[0]  == -1) tdc2[0] = TMath::QuietNaN(); 
    fp>>tdc2[1]; if(tdc2[1]  == -1) tdc2[1] = TMath::QuietNaN(); 
    fp>>tdc2[2]; if(tdc2[2]  == -1) tdc2[2] = TMath::QuietNaN(); 
    fp>>tdc2[3]; if(tdc2[3]  == -1) tdc2[3] = TMath::QuietNaN(); 
    fp>>tdc2[4]; if(tdc2[4]  == -1) tdc2[4] = TMath::QuietNaN(); 
    fp>>tdc2[5]; if(tdc2[5]  == -1) tdc2[5] = TMath::QuietNaN(); 
    fp>>tdc2[6]; if(tdc2[6]  == -1) tdc2[6] = TMath::QuietNaN(); 
    fp>>tdc2[7]; if(tdc2[7]  == -1) tdc2[7] = TMath::QuietNaN(); 
    fp>>tdc2[8]; if(tdc2[8]  == -1) tdc2[8] = TMath::QuietNaN(); 
    fp>>tdc2[9]; if(tdc2[9]  == -1) tdc2[9] = TMath::QuietNaN(); 
    fp>>tdc2[10];if(tdc2[10] == -1) tdc2[10] = TMath::QuietNaN();
    fp>>tdc2[11];if(tdc2[11] == -1) tdc2[11] = TMath::QuietNaN();
    //------------- Liquid 
    fp>>lrf[0];
    fp>>lrf[1];
    fp>>lrf[2];
    fp>>dummy;

    //===================================== 2ndary data processing
    //------------- Axuillary
    //grXC = grx - (20-3*grth*TMath::RadToDeg()) + 0.5e-2 * (gry; 
    grXC   = grx -380./2.3*grth*TMath::RadToDeg();
    grthC  = -0.4191*(grth-0.008770-0.52e-4*grx);
    grXAux = grXC - 140 * grthC* TMath::RadToDeg() + 50 ;  
    //_______________________________ Gate
   
    vetogate = 0;
    if( TMath::IsNaN(tdc2[8]) && TMath::IsNaN(tdc2[9])  ) vetogate += 1;
      
    //------------- RF
    grf = (grrf + rand->Uniform(0,1))*GR_CH2NS[0];
    brf = (lrf[1] + rand->Uniform(0,1))*LAS_CH2NS[0];
    //------------ GR 
    grdE1 = TMath::Sqrt(gradc[0]*gradc[1]);
    grdE2 = TMath::Sqrt(gradc[2]*gradc[3]);
  
    grT1avg = (grtdc[0]+36+grtdc[1])/2.*GR_CH2NS[0];
    grT2avg = (grtdc[2]+36+grtdc[3])/2.*GR_CH2NS[0];
  
    grTOF1 = grT1avg + (rand->Uniform(0,1)+rand->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
    grTOF2 = grT2avg + (rand->Uniform(0,1)+rand->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
      
    //----------- STACK 
    sta1h = TMath::Sqrt(adc[0]/112.*adc[1]/118.);
    sta2h = TMath::Sqrt(adc[2]/131.*adc[3]/109.);
    sta1v = TMath::Sqrt(adc[4]/116.*adc[5]/126.);
    sta2v = TMath::Sqrt(adc[6]/113.*adc[7]/123.);
    sta3v = TMath::Sqrt(adc[8]/122.*adc[9]/131.);
    sta4v = TMath::Sqrt(adc[10]/127.*adc[11]/124.);


    staMult = 0;
    
    sta_even = 6.4*(sta1h+sta2h);
    sta_odd = 6.4*(sta1v+sta2v+sta3v+sta4v);
    if (sta_even + sta_odd > 0.) {
      sta_sum = sta_even+sta_odd;
      sta_ratio = (sta_even-sta_odd)/(sta_even+sta_odd);
    }

    //printf("%f, %f \n", tdc[0],tdc[1]);
    sta1hTavg = (tdc[0] + tdc[1] )/2.*LAS_CH2NS[4];
    sta2hTavg = (tdc[2] + tdc[3] )/2.*LAS_CH2NS[5];
    sta1vTavg = (tdc[4] + tdc[5] )/2.*LAS_CH2NS[6];
    sta2vTavg = (tdc[6] + tdc[7] )/2.*LAS_CH2NS[7];
    sta3vTavg = (tdc[8] + tdc[9] )/2.*LAS_CH2NS[8];
    sta4vTavg = (tdc[10]+ tdc[11])/2.*LAS_CH2NS[9];

    sta1hTdif = (tdc[0] - tdc[1] + rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[4];   
    sta2hTdif = (tdc[2] - tdc[3] + rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[5];   
    sta1vTdif = (tdc[4] - tdc[5] + rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[6];   
    sta2vTdif = (tdc[6] - tdc[7] + rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[7];   
    sta3vTdif = (tdc[8] - tdc[9] + rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[8];   
    sta4vTdif = (tdc[10]- tdc[11]+ rand -> Uniform(0,1)- rand -> Uniform(0,1))*LAS_CH2NS[9];   
    
    sta1hTOF = sta1hTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[4] + STACK_TOF_OFFSET[0];
    sta2hTOF = sta2hTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[5] + STACK_TOF_OFFSET[1];
    sta1vTOF = sta1vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[6] + STACK_TOF_OFFSET[2];
    sta2vTOF = sta2vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[7] + STACK_TOF_OFFSET[3];
    sta3vTOF = sta3vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[8] + STACK_TOF_OFFSET[4];
    sta4vTOF = sta4vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[9] + STACK_TOF_OFFSET[5];
   
           
    //------------- Telescope      
    badEl = adc2[8];
    badEr = adc2[9];
    blo1 = TMath::Sqrt(adc2[0]*adc[1]);
    blo2 = TMath::Sqrt(adc2[2]*adc[3]);
    blo3 = TMath::Sqrt(adc2[4]*adc[5]);
    blo4 = TMath::Sqrt(adc2[6]*adc[7]);
      
    blo1Tavg = (tdc2[0]+tdc2[1])/2.*LAS_CH2NS[0];
    blo2Tavg = (tdc2[2]+tdc2[3])/2.*LAS_CH2NS[1];
    blo3Tavg = (tdc2[4]+tdc2[5])/2.*LAS_CH2NS[2];
    blo4Tavg = (tdc2[6]+tdc2[7])/2.*LAS_CH2NS[3];
      
    blo1TOF = TMath::QuietNaN();
    blo2TOF = TMath::QuietNaN();
    blo3TOF = TMath::QuietNaN();
    blo4TOF = TMath::QuietNaN();
   
    badElTOF = (tdc2[8] + rand->Uniform(0,1))*LAS_CH2NS[0] - brf ;
    badErTOF = (tdc2[9] + rand->Uniform(0,1))*LAS_CH2NS[0] - brf ;
    blo1TOF = blo1Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[0] + BLOCK_TOF_OFFSET[0];
    blo2TOF = blo2Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[1] + BLOCK_TOF_OFFSET[1];
    blo3TOF = blo3Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[2] + BLOCK_TOF_OFFSET[2];
    blo4TOF = blo4Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS[3] + BLOCK_TOF_OFFSET[3];
   
    //-------------- Liquid
    liqlf = adc2[10] + rand->Uniform(0,1);
    liqld = adc2[11] + rand->Uniform(0,1);
    liqrf = adc2[12] + rand->Uniform(0,1);
    liqrd = adc2[13] + rand->Uniform(0,1);
    liqlTOF = (tdc2[10] + rand->Uniform(0,1))*LAS_CH2NS[0] - brf ;
    liqrTOF = (tdc2[11] + rand->Uniform(0,1))*LAS_CH2NS[0] - brf ;

    //_______________________________ Physical
    /*
    //======= Stack X-pos
    Double_t staX1h = tdc[0]-tdc[1];
    Double_t staX2h = tdc[2]-tdc[3];
   
    //=======3He
    Double_t beta = GR_LENGTH/grTOF1;
    Double_t gamma = 1/TMath::Sqrt(1-beta*beta);
    Double_t redMomt = Brho * (1 + grx/X_D) * cVAC;
   
    //======= PID
    Double_t pidAOQ = Brho * (1 + grx/X_D) * cVAC / amu / beta/ gamma;
    Double_t pidZ   = beta* TMath::Sqrt(grdE1) ;// 1st approx
   
    //======= neutron
    Double_t beta_n = 0;
   
    //======= 3-body kinematics
    TVector3 vInc, vGR, vn;
    TLorentzVector pInc, pGR, pn;

    /**/
      
    //----------- Fill       
    f1->cd(); //set focus on this file
    t1->Fill(); 
      
    ///______________________________________________________________   
    clock.Stop("timer");
    Double_t time = clock.GetRealTime("timer");
    clock.Start("timer");

    if ( !shown ) {
      if (fmod(time, 10) < 1 ){
        printf( "%10d[%5.2f%%]|%6.1f|%3d min %5.2f sec\n", 
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

  }while(! fp.eof());

  f1->cd(); //set focus on this file
  f1->Write(); 
  f1->Close();
  
  printf("=======================================\n");
  printf("total number of event: %d \n", eventID);
  printf("............ done!\n");
}

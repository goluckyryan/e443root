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

  TString openFileName1;
  openFileName1.Form("./dstroot/run%04d_asci.dst", RunName);
  printf("input <====== %s \n", openFileName1.Data());
  
  TString openFileName2;
  openFileName2.Form("./dstroot/run%04d1_asci.dst", RunName); // 2nd file must be named runXXXX1_asci.dst
  printf("input <====== %s \n", openFileName2.Data());
  
  TString saveFileName;
  saveFileName.Form("run%04d.root", RunName);
  printf("output =====> %s , nEntries:%d\n", saveFileName.Data(), nEntries);
   
  TFile *f1 = new TFile (saveFileName,"recreate");
  TTree * t1 = new TTree("tree","tree");
  TRandom * rand = new TRandom();

  //================ Primary data
  Double_t gradc[4],grtdc[4],grrf,adc[12],adc2[14],tdc[12],tdc2[12],lrf[3],dummy1,dummy2;
  

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
  Double_t sta1hTem,sta2hTem,sta1vTem,sta2vTem,sta3vTem,sta4vTem;
  Double_t sta1hTavg,sta2hTavg,sta1vTavg,sta2vTavg,sta3vTavg,sta4vTavg;
  Double_t sta1hTdif,sta2hTdif,sta1vTdif,sta2vTdif,sta3vTdif,sta4vTdif;
  Int_t staMultH,staMultV;
  Int_t sta1hTri,sta2hTri,sta1vTri,sta2vTri,sta3vTri,sta4vTri;
  //--------- BAND liquid
  Double_t liqlf,liqld,liqrf,liqrd;
  //--------- GRRF & BANDRF
  Double_t grf,brf;
  Double_t lastgr;
  //-------------------------------------------------Physics
  Double_t grTOF1, grTOF2;
  Double_t badElTOF,badErTOF,blo1TOF,blo2TOF,blo3TOF,blo4TOF;
  Double_t sta1hTOF,sta2hTOF,sta1vTOF,sta2vTOF,sta3vTOF,sta4vTOF;
  Double_t sta1hTOFC,sta2hTOFC,sta1vTOFC,sta2vTOFC,sta3vTOFC,sta4vTOFC;
  Double_t liqlTOF,liqrTOF;
  Double_t sta_odd,sta_even,sta_ratio,sta_sum;

  Double_t grXC,grthC;
  //Double_t grXAux;
  Int_t event = -1;
  Double_t grx1;

  Double_t sta1hTOFS,sta2hTOFS,sta1vTOFS,sta2vTOFS,sta3vTOFS,sta4vTOFS;
  Double_t sta1hEk,sta2hEk,sta1vEk,sta2vEk,sta3vEk,sta4vEk;
  Double_t sta1hBeta,sta2hBeta,sta1vBeta,sta2vBeta,sta3vBeta,sta4vBeta;
  Double_t flightlength = 71;//cm



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
  //t1->Branch("grXAux", &grXAux, "grXAux/D");

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
  t1->Branch("adc2", adc2, "adc2[14]/D");
  t1->Branch("tdc", tdc, "tdc[12]/D");
  t1->Branch("tdc2", tdc2, "tdc2[12]/D");
  t1->Branch("staMultH", &staMultH, "staMultH/I");
  t1->Branch("staMultV", &staMultV, "staMultV/I");
  t1->Branch("sta1hTri", &sta1hTri, "sta1hTri/I");
  t1->Branch("sta2hTri", &sta2hTri, "sta2hTri/I");
  t1->Branch("sta1vTri", &sta1vTri, "sta1vTri/I");
  t1->Branch("sta2vTri", &sta2vTri, "sta2vTri/I");
  t1->Branch("sta3vTri", &sta3vTri, "sta3vTri/I");
  t1->Branch("sta4vTri", &sta4vTri, "sta4vTri/I");
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

  t1->Branch("sta1hTOFC", &sta1hTOFC, "sta1hTOFC/D");
  t1->Branch("sta2hTOFC", &sta2hTOFC, "sta2hTOFC/D");
  t1->Branch("sta1vTOFC", &sta1vTOFC, "sta1vTOFC/D");
  t1->Branch("sta2vTOFC", &sta2vTOFC, "sta2vTOFC/D");
  t1->Branch("sta3vTOFC", &sta3vTOFC, "sta3vTOFC/D");
  t1->Branch("sta4vTOFC", &sta4vTOFC, "sta4vTOFC/D");

  t1->Branch("liqlf", &liqlf, "liqlf/D");
  t1->Branch("liqld", &liqld, "liqld/D");
  t1->Branch("liqrf", &liqrf, "liqrf/D");
  t1->Branch("liqrd", &liqrd, "liqrd/D");
  t1->Branch("liqlTOF", &liqlTOF, "liqlTOF/D");
  t1->Branch("liqrTOF", &liqrTOF, "liqrTOF/D");

  t1->Branch("brf", &brf, "brf/D");
  //t1->Branch("event", &event, "event/I");
  //t1->Branch("grx1", &grx1, "grx1/D");
  t1->Branch("lastgr", &lastgr, "lastgr/D");

  t1->Branch("sta1hTOFS", &sta1hTOFS, "sta1hTOFS/D");
  t1->Branch("sta2hTOFS", &sta2hTOFS, "sta2hTOFS/D");
  t1->Branch("sta1vTOFS", &sta1vTOFS, "sta1vTOFS/D");
  t1->Branch("sta2vTOFS", &sta2vTOFS, "sta2vTOFS/D");
  t1->Branch("sta3vTOFS", &sta3vTOFS, "sta3vTOFS/D");
  t1->Branch("sta4vTOFS", &sta4vTOFS, "sta4vTOFS/D");
  t1->Branch("sta1hEk", &sta1hEk, "sta1hEk/D");
  t1->Branch("sta2hEk", &sta2hEk, "sta2hEk/D");
  t1->Branch("sta1vEk", &sta1vEk, "sta1vEk/D");
  t1->Branch("sta2vEk", &sta2vEk, "sta2vEk/D");
  t1->Branch("sta3vEk", &sta3vEk, "sta3vEk/D");
  t1->Branch("sta4vEk", &sta4vEk, "sta4vEk/D");
  t1->Branch("sta1hBeta", &sta1hBeta, "sta1hBeta/D");
  t1->Branch("sta2hBeta", &sta2hBeta, "sta2hBeta/D");
  t1->Branch("sta1vBeta", &sta1vBeta, "sta1vBeta/D");
  t1->Branch("sta2vBeta", &sta2vBeta, "sta2vBeta/D");
  t1->Branch("sta3vBeta", &sta3vBeta, "sta3vBeta/D");
  t1->Branch("sta4vBeta", &sta4vBeta, "sta4vBeta/D");
   
   
  //=================== read file
  ifstream fp1;
  fp1.open(openFileName1);

  if( !fp1.is_open() ) {
    printf("******* cannot open dst file : %s\n", openFileName1.Data());
    exit(-1);
  }

  ifstream fp2;
  fp2.open(openFileName2);

  if( !fp2.is_open() ) {
    printf("******* cannot open dst file : %s\n", openFileName2.Data());
    exit(-2);
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
   //grXAux = TMath::QuietNaN();

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
   sta1hTem = TMath::QuietNaN();
   sta2hTem = TMath::QuietNaN();
   sta1vTem = TMath::QuietNaN();
   sta2vTem = TMath::QuietNaN();
   sta3vTem = TMath::QuietNaN();
   sta4vTem = TMath::QuietNaN();
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

   sta1hTOFC = TMath::QuietNaN();
   sta2hTOFC = TMath::QuietNaN();
   sta1vTOFC = TMath::QuietNaN();
   sta2vTOFC = TMath::QuietNaN();
   sta3vTOFC = TMath::QuietNaN();
   sta4vTOFC = TMath::QuietNaN();

   sta_odd = TMath::QuietNaN();
   sta_even = TMath::QuietNaN();
   sta_ratio = TMath::QuietNaN();
   sta_sum = TMath::QuietNaN();
   sta1hTri = 0;
   sta2hTri = 0;
   sta1vTri = 0;
   sta2vTri = 0;
   sta3vTri = 0;
   sta4vTri = 0;
   staMultH = 0;
   staMultV = 0;
   //--------- BAND_LIQUID
   liqlf = TMath::QuietNaN();
   liqld = TMath::QuietNaN();
   liqrf = TMath::QuietNaN();
   liqrd = TMath::QuietNaN();
   //--------- RF
   grf = TMath::QuietNaN();
   brf = TMath::QuietNaN();
   //---------check
   event = -1;
   grx1 = TMath::QuietNaN();
   lastgr = TMath::QuietNaN();
   //--------physical
   sta1hTOFS = TMath::QuietNaN();
   sta2hTOFS = TMath::QuietNaN();
   sta1vTOFS = TMath::QuietNaN();
   sta2vTOFS = TMath::QuietNaN();
   sta3vTOFS = TMath::QuietNaN();
   sta4vTOFS = TMath::QuietNaN();
   sta1hBeta = TMath::QuietNaN();
   sta2hBeta = TMath::QuietNaN();
   sta1vBeta = TMath::QuietNaN();
   sta2vBeta = TMath::QuietNaN();
   sta3vBeta = TMath::QuietNaN();
   sta4vBeta = TMath::QuietNaN();
   sta1hEk = TMath::QuietNaN();
   sta2hEk = TMath::QuietNaN();
   sta1vEk = TMath::QuietNaN();
   sta2vEk = TMath::QuietNaN();
   sta3vEk = TMath::QuietNaN();
   sta4vEk = TMath::QuietNaN();

   //============================================= read from file
    //------------ CoinReg
    fp1>>coinReg;
    //----------- GR
    fp1>>grx;   if(grx == -1)  grx  = TMath::QuietNaN();
    fp1>>grth;  if(grth == -1) grth = TMath::QuietNaN();
    fp1>>gry;   if(gry == -1)  gry  = TMath::QuietNaN();
    fp1>>grph;  if(grph == -1) grph = TMath::QuietNaN();    
    //----------- GR ADC, TDC
    fp1>>gradc[0]; 
    fp1>>gradc[1];
    fp1>>gradc[2];
    fp1>>gradc[3];
    fp1>>grtdc[0]; if(grtdc[0] == -1) grtdc[0] = TMath::QuietNaN();
    fp1>>grtdc[1]; if(grtdc[1] == -1) grtdc[1] = TMath::QuietNaN();  
    fp1>>grtdc[2]; if(grtdc[2] == -1) grtdc[2] = TMath::QuietNaN();  
    fp1>>grtdc[3]; if(grtdc[3] == -1) grtdc[3] = TMath::QuietNaN();  
    //------------ GR rf
    fp1>>grrf;
    //------------ Stack ADC
    fp1>>adc[0];  //if(adc[0]  < 5) adc[0] = TMath::QuietNaN();
    fp1>>adc[1];  //if(adc[1]  < 5) adc[1] = TMath::QuietNaN();
    fp1>>adc[2];  //if(adc[2]  < 5) adc[2] = TMath::QuietNaN();
    fp1>>adc[3];  //if(adc[3]  < 5) adc[3] = TMath::QuietNaN();
    fp1>>adc[4];  //if(adc[4]  < 5) adc[4] = TMath::QuietNaN();
    fp1>>adc[5];  //if(adc[5]  < 5) adc[5] = TMath::QuietNaN();
    fp1>>adc[6];  //if(adc[6]  < 5) adc[6] = TMath::QuietNaN();
    fp1>>adc[7];  //if(adc[7]  < 5) adc[7] = TMath::QuietNaN();
    fp1>>adc[8];  //if(adc[8]  < 5) adc[8] = TMath::QuietNaN();
    fp1>>adc[9];  //if(adc[9]  < 5) adc[9] = TMath::QuietNaN();
    fp1>>adc[10]; //if(adc[10] < 5) adc[10] = TMath::QuietNaN();
    fp1>>adc[11]; //if(adc[11] < 5) adc[11] = TMath::QuietNaN();
    //------------ BAND telesope ADC, Liquid ADC
    fp1>>adc2[0];
    fp1>>adc2[1];
    fp1>>adc2[2];
    fp1>>adc2[3];
    fp1>>adc2[4];
    fp1>>adc2[5];
    fp1>>adc2[6];
    fp1>>adc2[7];
    fp1>>adc2[8];
    fp1>>adc2[9];
    fp1>>adc2[10];
    fp1>>adc2[11];
    fp1>>adc2[12];
    fp1>>adc2[13];
    //------------ STACK TDC
    fp1>>tdc[0];  if(tdc[0]  == -1) tdc[0] = TMath::QuietNaN();
    fp1>>tdc[1];  if(tdc[1]  == -1) tdc[1] = TMath::QuietNaN();
    fp1>>tdc[2];  if(tdc[2]  == -1) tdc[2] = TMath::QuietNaN();
    fp1>>tdc[3];  if(tdc[3]  == -1) tdc[3] = TMath::QuietNaN();
    fp1>>tdc[4];  if(tdc[4]  == -1) tdc[4] = TMath::QuietNaN();
    fp1>>tdc[5];  if(tdc[5]  == -1) tdc[5] = TMath::QuietNaN();
    fp1>>tdc[6];  if(tdc[6]  == -1) tdc[6] = TMath::QuietNaN();
    fp1>>tdc[7];  if(tdc[7]  == -1) tdc[7] = TMath::QuietNaN();
    fp1>>tdc[8];  if(tdc[8]  == -1) tdc[8] = TMath::QuietNaN();
    fp1>>tdc[9];  if(tdc[9]  == -1) tdc[9] = TMath::QuietNaN();
    fp1>>tdc[10]; if(tdc[10] == -1) tdc[10] = TMath::QuietNaN();
    fp1>>tdc[11]; if(tdc[11] == -1) tdc[11] = TMath::QuietNaN();
    //------------- BAND telesope, Liquid TDC
    fp1>>tdc2[0]; if(tdc2[0]  == -1) tdc2[0] = TMath::QuietNaN(); 
    fp1>>tdc2[1]; if(tdc2[1]  == -1) tdc2[1] = TMath::QuietNaN(); 
    fp1>>tdc2[2]; if(tdc2[2]  == -1) tdc2[2] = TMath::QuietNaN(); 
    fp1>>tdc2[3]; if(tdc2[3]  == -1) tdc2[3] = TMath::QuietNaN(); 
    fp1>>tdc2[4]; if(tdc2[4]  == -1) tdc2[4] = TMath::QuietNaN(); 
    fp1>>tdc2[5]; if(tdc2[5]  == -1) tdc2[5] = TMath::QuietNaN(); 
    fp1>>tdc2[6]; if(tdc2[6]  == -1) tdc2[6] = TMath::QuietNaN(); 
    fp1>>tdc2[7]; if(tdc2[7]  == -1) tdc2[7] = TMath::QuietNaN(); 
    fp1>>tdc2[8]; if(tdc2[8]  == -1) tdc2[8] = TMath::QuietNaN(); 
    fp1>>tdc2[9]; if(tdc2[9]  == -1) tdc2[9] = TMath::QuietNaN(); 
    fp1>>tdc2[10];if(tdc2[10] == -1) tdc2[10] = TMath::QuietNaN();
    fp1>>tdc2[11];if(tdc2[11] == -1) tdc2[11] = TMath::QuietNaN();
    //------------- LAS RF 
    fp1>>lrf[0];
    fp1>>lrf[1];
    fp1>>lrf[2];
    fp1>>dummy1;

    //------------- additional variables in 2nd dst file
    fp2>>event;
    fp2>>grx1;  if (grx1 == -1) grx1 = TMath::QuietNaN();
    fp2>>lastgr;
    fp2>>dummy2;

    //===================== Check matching

    if( TMath::Finite(grx) &&  grx1 != grx) {
      //printf(" grx not matching at eventID:%d, abort. \n", eventID);
      //exit(-3);
      printf(" grx not matching (%f:%f) at eventID:%d, fill with NAN\n", grx, grx1, eventID);
      lastgr = TMath::QuietNaN();
    }

    //===================================== 2ndary data processing
    //------------- Axuillary
    //grXC = grx -20 ; //same as graf_conv 
    //grthC  = -0.4191*(grth-8.770e-3-3.591e-5*grx); //same as graf_conv, should be equal to incident theta
    grthC  = -0.4191*(grth-8.770e-3-0.52e-4*grx);
    //grXC = grx -380./2.3*grthC*TMath::RadToDeg();//#1035 corrected for (d,3He)
    grXC = grx -10.*grthC*TMath::RadToDeg();//#1034 corrected for 12C(d,3He)
    //grXAux = grXC -380./2.3*grthC*TMath::RadToDeg();
    //_______________________________ Gate
   
    vetogate = 0;
    if( TMath::IsNaN(tdc2[8]) && TMath::IsNaN(tdc2[9])  ) vetogate += 1;
      
    //------------- RF
    grf = (grrf + rand->Uniform(0,1))*GR_CH2NS[0];
    //brf = (lrf[1] + rand->Uniform(0,1))*LAS_CH2NS;
    brf = lrf[1];
    //------------ GR 
    grdE1 = TMath::Sqrt(gradc[0]*gradc[1]);
    grdE2 = TMath::Sqrt(gradc[2]*gradc[3]);
  
    grT1avg = (grtdc[0]+36+grtdc[1])/2.*GR_CH2NS[0];
    grT2avg = (grtdc[2]+36+grtdc[3])/2.*GR_CH2NS[0];
  
    grTOF1 = grT1avg + (rand->Uniform(0,1)+rand->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
    grTOF2 = grT2avg + (rand->Uniform(0,1)+rand->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
      
    //----------- STACK 

    if (!TMath::IsNaN(tdc[0]) && !TMath::IsNaN(tdc[1])) sta1hTri = 1;
    if (!TMath::IsNaN(tdc[2]) && !TMath::IsNaN(tdc[3])) sta2hTri = 1;
    if (!TMath::IsNaN(tdc[4]) && !TMath::IsNaN(tdc[5])) sta1vTri = 1;
    if (!TMath::IsNaN(tdc[6]) && !TMath::IsNaN(tdc[7])) sta2vTri = 1;
    if (!TMath::IsNaN(tdc[8]) && !TMath::IsNaN(tdc[9])) sta3vTri = 1;
    if (!TMath::IsNaN(tdc[10]) && !TMath::IsNaN(tdc[11])) sta4vTri = 1;

    staMultH = sta1hTri +sta2hTri;
    staMultV = sta1vTri +sta2vTri + sta3vTri +sta4vTri;

    //rough calibration only for slewing corrrection
    sta1hTem = 6.4 * TMath::Sqrt(adc[0]/112.*adc[1]/118.);
    sta2hTem = 6.4 * TMath::Sqrt(adc[2]/131.*adc[3]/109.);
    sta1vTem = 6.4 * TMath::Sqrt(adc[4]/116.*adc[5]/126.);
    sta2vTem = 6.4 * TMath::Sqrt(adc[6]/113.*adc[7]/123.);
    sta3vTem = 6.4 * TMath::Sqrt(adc[8]/122.*adc[9]/131.);
    sta4vTem = 6.4 * TMath::Sqrt(adc[10]/127.*adc[11]/124.);

    if (adc[0]>0 && adc[1]>0){ 
    sta1h = 6.4 * TMath::Sqrt(adc[0]/112.*adc[1]/118.)/0.8130;
    }else if (adc[0]<0 && adc[1]<0){
    sta1h = -6.4 * TMath::Sqrt(adc[0]/112.*adc[1]/118.)/0.8130;
    }else{
      //sta1h = 6.4 * (adc[0]/112.+adc[1]/118.)/2./0.83338;
    //sta1h = 6.4 * TMath::Sqrt(-adc[0]/117.*adc[1]/120.)/0.83338 * (0.5 + rand->Uniform(0,1));
    sta1h=6.4 * TMath::Sqrt(-adc[0]/112.*adc[1]/118.)/0.8130;if ( (adc[0]/112.+adc[1]/118.) < 0) sta1h=-sta1h;
    }

    if (adc[2]>0 && adc[3]>0){ 
    sta2h = 6.4 * TMath::Sqrt(adc[2]/131.*adc[3]/109.)/0.8748 ;
    }else if (adc[2]<0 && adc[3]<0){
    sta2h = -6.4 * TMath::Sqrt(adc[2]/131.*adc[3]/109.)/0.8748 ;
    }else{
      //sta2h = 6.4 * (adc[2]/132.+adc[3]/114.)/2./0.9258 ;
    //sta2h = 6.4 * TMath::Sqrt(-adc[2]/132.*adc[3]/114.)/0.9258 * (0.5 + rand->Uniform(0,1));
    sta2h = 6.4 * TMath::Sqrt(-adc[2]/131.*adc[3]/109.)/0.8748;if ( (adc[2]/131.+adc[3]/109.) < 0) sta2h=-sta2h;
    }

    if (adc[4]>0 && adc[5]>0){ 
    sta1v = 6.4 * TMath::Sqrt(adc[4]/116.*adc[5]/126.)/0.8267;
    }else if (adc[4]<0 && adc[5]<0){
    sta1v = -6.4 * TMath::Sqrt(adc[4]/116.*adc[5]/126.)/0.8267;
    }else{
      //sta1v = 6.4 * (adc[4]/118.+adc[5]/120.)/2./0.837  ;
    //sta1v = 6.4 * TMath::Sqrt(-adc[4]/118.*adc[5]/120.)/0.837 * (0.5 + rand->Uniform(0,1));
    sta1v = 6.4 * TMath::Sqrt(-adc[4]/116.*adc[5]/126.)/0.8267;if ( (adc[4]/116.+adc[5]/126.) < 0) sta1v=-sta1v;
    }

    if (adc[6]>0 && adc[7]>0){ 
    sta2v = 6.4 * TMath::Sqrt(adc[6]/113.*adc[7]/123.)/0.7382;
    }else if (adc[6]<0 && adc[7]<0){
    sta2v = -6.4 * TMath::Sqrt(adc[6]/113.*adc[7]/123.)/0.7382;
    }else{
      //sta2v = 6.4 * (adc[6]/123.+adc[7]/130.)/2./0.72416;
    //sta2v = 6.4 * TMath::Sqrt(-adc[6]/123.*adc[7]/130.)/0.72416 * (0.5 + rand->Uniform(0,1));
    sta2v = 6.4 * TMath::Sqrt(-adc[6]/113.*adc[7]/123.)/0.7382;if ( (adc[6]/113.+adc[7]/123.) < 0) sta2v=-sta2v;
    }

    //if (0 < tdc[8] && tdc[8] < 400 && 0 < tdc[9] && tdc[9] < 400){ // offset #1223 proton case
    //sta3v = (6.4 * TMath::Sqrt(adc[8]/122.*adc[9]/131.)+3.0341)/0.7383;
    //}else 
    if (adc[8]>0 && adc[9]>0){ 
    sta3v = 6.4 * TMath::Sqrt(adc[8]/122.*adc[9]/131.)/0.7383;
    }else if (adc[8]<0 && adc[9]<0){
    sta3v = -6.4 * TMath::Sqrt(adc[8]/122.*adc[9]/131.)/0.7383;
    }else{
      //sta3v = 6.4 * (adc[8]/120.+adc[9]/129.)/2./0.8517 ;
    //sta3v = 6.4 * TMath::Sqrt(-adc[8]/120.*adc[9]/129.)/0.8517 * (0.5 + rand->Uniform(0,1));
    sta3v = 6.4 * TMath::Sqrt(-adc[8]/122.*adc[9]/131.)/0.7383;if ( (adc[8]/122.+adc[9]/131.) < 0) sta3v=-sta3v;
    }
    
    //if (0 < tdc[10] && tdc[10] < 400 && 0 < tdc[11] && tdc[11] < 400){ // offset #1223 proton case
    //sta4v = (6.4 * TMath::Sqrt(adc[10]/127.*adc[11]/124.)+2.5245)/0.7538;
    //}else 
    if (adc[10]>0 && adc[11]>0){ 
    sta4v = 6.4 * TMath::Sqrt(adc[10]/127.*adc[11]/124.)/0.7538;
    }else if (adc[10]<0 && adc[11]<0){
    sta4v = -6.4 * TMath::Sqrt(adc[10]/127.*adc[11]/124.)/0.7538;
    }else{
      //sta4v = 6.4 * (adc[10]/126.+adc[11]/132.)/2./0.83242;
    //sta4v = 6.4 * TMath::Sqrt(-adc[10]/126.*adc[11]/132.)/0.83242 * (0.5 + rand->Uniform(0,1));
    sta4v = 6.4 * TMath::Sqrt(-adc[10]/127.*adc[11]/124.)/0.7538;if ( (adc[10]/127.+adc[11]/124.) < 0) sta4v=-sta4v;
    }
    

    sta_even = sta1h+sta2h;
    sta_odd = (sta1v+sta2v+sta3v+sta4v)/1.08;     
    if (!( sta_odd + sta_even == 0.)) {
      sta_sum = sta_odd+sta_even;    
      sta_ratio = (sta_even-sta_odd)/(sta_even+sta_odd);
    }

    //printf("%f, %f \n", tdc[0],tdc[1]);

    //--------------RF:tdc slope fitting in #1035
    sta1hTavg = (tdc[0] * 0.9912 + tdc[1] * 0.9925 )/2. * LAS_CH2NS;     
    sta2hTavg = (tdc[2] * 1.091  + tdc[3] * 0.9946 )/2. * LAS_CH2NS; 
    sta1vTavg = (tdc[4] * 0.9889 + tdc[5] * 0.9865 )/2. * LAS_CH2NS;  
    sta2vTavg = (tdc[6] * 0.9886 + tdc[7] * 0.9916 )/2. * LAS_CH2NS;  
    sta3vTavg = (tdc[8] * 0.9758 + tdc[9] * 0.9936 )/2. * LAS_CH2NS;  
    sta4vTavg = (tdc[10]* 0.9832 + tdc[11]* 0.9924 )/2. * LAS_CH2NS;  

   //--------------RF:tdc slope fitting in #1223
    //sta1hTavg = (tdc[0] * 0.9831 + tdc[1] * 0.9881 )/2. * LAS_CH2NS;     
    //sta2hTavg = (tdc[2] * 1.091  + tdc[3] * 0.9865 )/2. * LAS_CH2NS; 
    //sta1vTavg = (tdc[4] * 0.9896 + tdc[5] * 0.9856 )/2. * LAS_CH2NS;  
    //sta2vTavg = (tdc[6] * 0.9864 + tdc[7] * 0.9916 )/2. * LAS_CH2NS;  
    //sta3vTavg = (tdc[8] * 0.9761 + tdc[9] * 0.9920 )/2. * LAS_CH2NS;  
    //sta4vTavg = (tdc[10]* 0.9807 + tdc[11]* 0.9912 )/2. * LAS_CH2NS;  
    
    sta1hTdif = (tdc[0] * 0.9912 - tdc[1] * 0.9925 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    sta2hTdif = (tdc[2] * 1.091  - tdc[3] * 0.9946 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    sta1vTdif = (tdc[4] * 0.9889 - tdc[5] * 0.9865 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    sta2vTdif = (tdc[6] * 0.9886 - tdc[7] * 0.9916 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    sta3vTdif = (tdc[8] * 0.9758 - tdc[9] * 0.9936 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    sta4vTdif = (tdc[10]* 0.9832 - tdc[11]* 0.9924 + rand -> Uniform(0,1) - rand -> Uniform(0,1))*LAS_CH2NS;   
    
    sta1hTOF = sta1hTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;
    sta2hTOF = sta2hTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;
    sta1vTOF = sta1vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;
    sta2vTOF = sta2vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;
    sta3vTOF = sta3vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;
    sta4vTOF = sta4vTavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS;

    //slewing correction only for #1035
    sta1hTOFC = sta1hTOF - 3.173068 / TMath::Sqrt(sta1hTem - 0.097028) + STACK_TOF_OFFSET[0];
    sta2hTOFC = sta2hTOF - 0.974083 / TMath::Sqrt(sta2hTem - 0.977325) + STACK_TOF_OFFSET[1];
    sta1vTOFC = sta1vTOF - 0.807234 / TMath::Sqrt(sta1vTem - 0.828953) + STACK_TOF_OFFSET[2];
    sta2vTOFC = sta2vTOF - 1.901365 / TMath::Sqrt(sta2vTem - 0.736045) + STACK_TOF_OFFSET[3];
    sta3vTOFC = sta3vTOF - 0.705086 / TMath::Sqrt(sta3vTem - 0.995722) + STACK_TOF_OFFSET[4];
    sta4vTOFC = sta4vTOF - 0.936808 / TMath::Sqrt(sta4vTem - 0.960274) + STACK_TOF_OFFSET[5];
           
    //------------- Telescope      
    badEl = adc2[8];
    badEr = adc2[9];
    blo1 = TMath::Sqrt(adc2[0]*adc2[1]);
    blo2 = TMath::Sqrt(adc2[2]*adc2[3]);
    blo3 = TMath::Sqrt(adc2[4]*adc2[5]);
    blo4 = TMath::Sqrt(adc2[6]*adc2[7]);
      
    blo1Tavg = (tdc2[0]+tdc2[1])/2.*LAS_CH2NS;
    blo2Tavg = (tdc2[2]+tdc2[3])/2.*LAS_CH2NS;
    blo3Tavg = (tdc2[4]+tdc2[5])/2.*LAS_CH2NS;
    blo4Tavg = (tdc2[6]+tdc2[7])/2.*LAS_CH2NS;
      
    blo1TOF = TMath::QuietNaN();
    blo2TOF = TMath::QuietNaN();
    blo3TOF = TMath::QuietNaN();
    blo4TOF = TMath::QuietNaN();
   
    badElTOF = (tdc2[8] + rand->Uniform(0,1))*LAS_CH2NS - brf ;
    badErTOF = (tdc2[9] + rand->Uniform(0,1))*LAS_CH2NS - brf ;
    blo1TOF = blo1Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS + BLOCK_TOF_OFFSET[0];
    blo2TOF = blo2Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS + BLOCK_TOF_OFFSET[1];
    blo3TOF = blo3Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS + BLOCK_TOF_OFFSET[2];
    blo4TOF = blo4Tavg + ((rand->Uniform(0,1)+rand->Uniform(0,1))/2. - (lrf[1] + rand->Uniform(0,1)))*LAS_CH2NS + BLOCK_TOF_OFFSET[3];
   
    //-------------- Liquid
    liqlf = adc2[10] + rand->Uniform(0,1);
    liqld = adc2[11] + rand->Uniform(0,1);
    liqrf = adc2[12] + rand->Uniform(0,1);
    liqrd = adc2[13] + rand->Uniform(0,1);
    liqlTOF = (tdc2[10] + rand->Uniform(0,1))*LAS_CH2NS - brf ;
    liqrTOF = (tdc2[11] + rand->Uniform(0,1))*LAS_CH2NS - brf ;

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

    */

  /*  //======proton beam:  kinetic energy by TOF 
      
    if(sta1hTOFC < 100){
      sta1hTOFS = sta1hTOFC - 59.364;// p beam
    }else{
      sta1hTOFS = sta1hTOFC - 59.364*2;
    }

    if(sta2hTOFC < 100){
      sta2hTOFS = sta2hTOFC - 59.364;// p beam
    }else{
      sta2hTOFS = sta2hTOFC - 59.364*2;
    }

    if(sta1vTOFC < 100){
      sta1vTOFS = sta1vTOFC - 59.364;// p beam
    }else{
      sta1vTOFS = sta1vTOFC - 59.364*2;
    }

    if(sta2vTOFC < 100){
      sta2vTOFS = sta2vTOFC - 59.364;// p beam
    }else{
      sta2vTOFS = sta2vTOFC - 59.364*2;
    }

    if(sta3vTOFC < 100){
      sta3vTOFS = sta3vTOFC - 59.364;// p beam
    }else{
      sta3vTOFS = sta3vTOFC - 59.364*2;
    }

    if(sta4vTOFC < 100){
      sta4vTOFS = sta4vTOFC - 59.364;// p beam
    }else{
      sta4vTOFS = sta4vTOFC - 59.364*2;
    }*/

      //======deuteron beam:  kinetic energy by TOF 
    if(sta1hTOFC < 100){
      sta1hTOFS = sta1hTOFC;// d beam
    }else{
      sta1hTOFS = sta1hTOFC - 98.850;
    }

    if(sta2hTOFC < 100){
      sta2hTOFS = sta2hTOFC;// d beam
    }else{
      sta2hTOFS = sta2hTOFC - 98.850;
    }

    if(sta1vTOFC < 100){
      sta1vTOFS = sta1vTOFC;// d beam
    }else{
      sta1vTOFS = sta1vTOFC - 98.850;
    }

    if(sta2vTOFC < 100){
      sta2vTOFS = sta2vTOFC;// d beam
    }else{
      sta2vTOFS = sta2vTOFC - 98.850;
    }

    if(sta3vTOFC < 100){
      sta3vTOFS = sta3vTOFC;// d beam
    }else{
      sta3vTOFS = sta3vTOFC - 98.850;
    }

    if(sta4vTOFC < 100){
      sta4vTOFS = sta4vTOFC;// d beam
    }else{
      sta4vTOFS = sta4vTOFC - 98.850;
      }

    sta1hBeta = TMath::Sqrt(flightlength*flightlength+4*4)/sta1hTOFS/29.979;
    sta2hBeta = TMath::Sqrt(flightlength*flightlength+4*4)/sta2hTOFS/29.979;
    sta1vBeta = TMath::Sqrt(flightlength*flightlength+12*12)/sta1vTOFS/29.979;
    sta2vBeta = TMath::Sqrt(flightlength*flightlength+4*4)/sta2vTOFS/29.979;
    sta3vBeta = TMath::Sqrt(flightlength*flightlength+4*4)/sta3vTOFS/29.979;
    sta4vBeta = TMath::Sqrt(flightlength*flightlength+12*12)/sta4vTOFS/29.979;

    sta1hEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta1hBeta,2)) - 1);
    sta2hEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta2hBeta,2)) - 1);
    sta1vEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta1vBeta,2)) - 1);
    sta2vEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta2vBeta,2)) - 1);
    sta3vEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta3vBeta,2)) - 1);
    sta4vEk = mn * (1/TMath::Sqrt(1-TMath::Power(sta4vBeta,2)) - 1);




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

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
#include "TBenchmark.h"
#include "constant.h"

using namespace std;

void dst2root(TString openFileName){ //the file name should be "../dstroot/runXXXX_asci.dst"
   gROOT->ProcessLine(".!date");
   gStyle->SetOptStat(0);
   
   TString saveFileName = openFileName;
   printf("%s \n", saveFileName.Data());
   saveFileName.Remove(0,11); //remove head
   Int_t len = saveFileName.Length();
   saveFileName.Remove(len  - 9); // remove tail
   saveFileName = saveFileName + ".root";
   printf("%s \n", saveFileName.Data());
   
   TFile *f1 = new TFile (saveFileName,"recreate");
   TTree * t1 = new TTree("tree","tree");
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

   Double_t grXC;


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

   t1->Branch("sta_odd", &sta_odd, "sta_odd/D");
   t1->Branch("sta_even", &sta_even, "sta_even/D");
   t1->Branch("sta_ratio", &sta_ratio, "sta_ratio/D");
   t1->Branch("sta1hTavg", &sta1hTavg, "sta1hTavg/D");
   t1->Branch("sta2hTavg", &sta2hTavg, "sta2hTavg/D");
   t1->Branch("sta1vTavg", &sta1vTavg, "sta1vTavg/D");
   t1->Branch("sta2vTavg", &sta2vTavg, "sta2vTavg/D");
   t1->Branch("sta3vTavg", &sta3vTavg, "sta3vTavg/D");
   t1->Branch("sta4vTavg", &sta4vTavg, "sta4vTavg/D");
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

   string line;
   Int_t gradc[4],grtdc[4],grrf,adc[12],adc2[14],tdc[12],tdc2[12],lrf[3],dummy;
   
   TBenchmark clock;
   Bool_t shown;
   
   clock.Reset();
   clock.Start("timer");
   shown = 0;
   
   do{
   
      
      //if( eventID > 10) break;
      
      eventID ++;
      
      fp>>coinReg;
      
      //----------- GR
      fp>>grx;
      fp>>gry;
      fp>>grth;
      fp>>grph;       
      //----------- GR ADC, TDC
      fp>>gradc[0]; 
      fp>>gradc[1];
      fp>>gradc[2];
      fp>>gradc[3];
      fp>>grtdc[0];
      fp>>grtdc[1];
      fp>>grtdc[2];
      fp>>grtdc[3];
      //------------ GR rf
      fp>>grrf;
      //------------ Stack ADC
      fp>>adc[0];
      fp>>adc[1];
      fp>>adc[2];
      fp>>adc[3];
      fp>>adc[4];
      fp>>adc[5];
      fp>>adc[6];
      fp>>adc[7];
      fp>>adc[8];
      fp>>adc[9];
      fp>>adc[10];
      fp>>adc[11];
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
      fp>>tdc[0];
      fp>>tdc[1];
      fp>>tdc[2];
      fp>>tdc[3];
      fp>>tdc[4];
      fp>>tdc[5];
      fp>>tdc[6];
      fp>>tdc[7];
      fp>>tdc[8];
      fp>>tdc[9];
      fp>>tdc[10];
      fp>>tdc[11];
      //------------- BAND telesope, Liquid TDC
      fp>>tdc2[0];
      fp>>tdc2[1];
      fp>>tdc2[2];
      fp>>tdc2[3];
      fp>>tdc2[4];
      fp>>tdc2[5];
      fp>>tdc2[6];
      fp>>tdc2[7];
      fp>>tdc2[8];
      fp>>tdc2[9];
      fp>>tdc2[10];
      fp>>tdc2[11];
      //------------- Liquid 
      fp>>lrf[0];
      fp>>lrf[1];
      fp>>lrf[2];
      //fp>>dummy;
      
      //================================ 2ndary data processing
      //------------- Axuillary 
      grXC = grx - 380./2.3*grth*TMath::RadToDeg();
      
      //_______________________________ Gate
   
      vetogate = 0;
      if((tdc2[8]<0)&&tdc2[9]<0) vetogate += 1;
      
      //------------- RF
      grf = (grrf + gRandom->Uniform(0,1))*GR_CH2NS[0];
      brf = (lrf[1] + gRandom->Uniform(0,1))*LAS_CH2NS[0];
      //------------ GR 
      grdE1 = TMath::Sqrt(gradc[0]*gradc[1]);
      grdE2 = TMath::Sqrt(gradc[2]*gradc[3]);
  
      grT1avg = (grtdc[0]+36+grtdc[1])/2.*GR_CH2NS[0];
      grT2avg = (grtdc[2]+36+grtdc[3])/2.*GR_CH2NS[0];
  
      grTOF1 = grT1avg + (gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
      grTOF2 = grT2avg + (gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
      
      //----------- STACK 
      sta1h = sqrt(adc[0]/112.*adc[1]/118.);
      sta2h = sqrt(adc[2]/131.*adc[3]/109.);
      sta1v = sqrt(adc[4]/116.*adc[5]/126.);
      sta2v = sqrt(adc[6]/113.*adc[7]/123.);
      sta3v = sqrt(adc[8]/122.*adc[9]/131.);
      sta4v = sqrt(adc[10]/127.*adc[11]/124.);
      
      sta_even = 6.4*(sta1h+sta2h);
      sta_odd = 6.4*(sta1v+sta2v+sta3v+sta4v);
      if (sta_even + sta_odd > 0.) {
          sta_sum = sta_even+sta_odd;
          sta_ratio = (sta_even-sta_odd)/(sta_even+sta_odd);
      }
      
      sta1hTavg = (tdc[0] + tdc[1] )/2.*LAS_CH2NS[4];
      sta2hTavg = (tdc[2] + tdc[3] )/2.*LAS_CH2NS[5];
      sta1vTavg = (tdc[4] + tdc[5] )/2.*LAS_CH2NS[6];
      sta2vTavg = (tdc[6] + tdc[7] )/2.*LAS_CH2NS[7];
      sta3vTavg = (tdc[8] + tdc[9] )/2.*LAS_CH2NS[8];
      sta4vTavg = (tdc[10]+ tdc[11])/2.*LAS_CH2NS[9];
      
      if(sta1hTavg != -LAS_CH2NS[4] ) sta1hTOF = sta1hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[4] + STACK_TOF_OFFSET[0];
      if(sta2hTavg < 700*LAS_CH2NS[5] ) sta2hTOF = sta2hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[5] + STACK_TOF_OFFSET[1];
      if(sta1vTavg != -LAS_CH2NS[6] ) sta1vTOF = sta1vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[6] + STACK_TOF_OFFSET[2];
      if(sta2vTavg != -LAS_CH2NS[7] ) sta2vTOF = sta2vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[7] + STACK_TOF_OFFSET[3];
      if(sta3vTavg != -LAS_CH2NS[8] ) sta3vTOF = sta3vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[8] + STACK_TOF_OFFSET[4];
      if(sta4vTavg != -LAS_CH2NS[9] ) sta4vTOF = sta4vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[9] + STACK_TOF_OFFSET[5];
   
           
      //------------- Telescope      
      badEl = adc2[8];
      badEr = adc2[9];
      blo1 = TMath::Sqrt(adc2[0]*adc[1]);
      blo2 = TMAth::Sqrt(adc2[2]*adc[3]);
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
   
      if( tdc2[8] != -1 ) badElTOF = (tdc2[8] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
      if( tdc2[9] != -1 ) badErTOF = (tdc2[9] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
      if(blo1Tavg != -LAS_CH2NS[0] )blo1TOF = blo1Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[0] + BLOCK_TOF_OFFSET[0];
      if(blo2Tavg != -LAS_CH2NS[1] )blo2TOF = blo2Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[1] + BLOCK_TOF_OFFSET[1];
      if(blo3Tavg != -LAS_CH2NS[2] )blo3TOF = blo3Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[2] + BLOCK_TOF_OFFSET[2];
      if(blo4Tavg != -LAS_CH2NS[3] )blo4TOF = blo4Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[3] + BLOCK_TOF_OFFSET[3];
   
   
      //-------------- Liquid
      liqlf = adc2[10];
      liqld = adc2[11];
      liqrf = adc2[12];
      liqrd = adc2[13];
      liqlTOF = (tdc2[10] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
      liqrTOF = (tdc2[11] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
      
      //----------- Fill       
      f1->cd(); //set focus on this file
      t1->Fill(); 
      
      ///______________________________________________________________   
      clock.Stop("timer");
      Double_t time = clock.GetRealTime("timer");
      clock.Start("timer");

      if ( !shown ) {
         if (fmod(time, 10) < 1 ){
            printf( "%10d[%2d%%]|%3d min %5.2f sec | expect:%5.1fmin %10d\n", 
               entry, 
               TMath::Nint((entry+1)*100./totnumEntry),
               TMath::Floor(time/60), time - TMath::Floor(time/60)*60,
               totnumEntry*time/(entry+1)/60.,
               count);
               shown = 1;
         }
      }else{
         if (fmod(time, 10) > 9 ){
            shown = 0;
         }
      }

    
//nextLine:
//      getline(fp, line, '\r'); // read fp and store in line, and next line


   }while(! fp.eof());
   
  
  printf("=======================================\n");
   
   
   /*//======================================================================
   TFile *f1 = new TFile ("test.root","read");
   TTree *t1 = (TTree*)f1->Get("t1");
   
   */
   
   printf("............ done!\n");
}

#define Selector_Sort_cxx
// The class definition in Selector_Sort.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("Selector_Sort.C")
// Root > T->Process("Selector_Sort.C","some options")
// Root > T->Process("Selector_Sort.C+")
//

#include "Selector_Sort.h"
#include <TH2.h>
#include <TStyle.h>


void Selector_Sort::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void Selector_Sort::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Selector_Sort::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Selector_Sort::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   
   
   //_______________________________ Setting
   //Double_t Brho = BRHO_3HE;
   
   
   
  //_______________________________ Get Branch Entry
   b_event->GetEntry(entry,0);
   b_grx->GetEntry(entry,0);
   b_gry->GetEntry(entry,0);
   b_grth->GetEntry(entry,0);
   b_grph->GetEntry(entry,0);
   b_gradc->GetEntry(entry,0);
   b_grtdc->GetEntry(entry,0);
   b_grrf->GetEntry(entry,0);
   b_adc->GetEntry(entry,0);
   b_adc2->GetEntry(entry,0);
   b_tdc->GetEntry(entry,0);
   b_tdc2->GetEntry(entry,0);
   b_lrf->GetEntry(entry,0);
   
   //_______________________________ Filter primary data
   //if( coinReg <6 ) return kTRUE;
   if( grth == -1) return kTRUE;
   //count ++;

   //_______________________________ eventID
   eventID = entry;

   //_______________________________  coinReg
   //coinReg = TMath::QuietNaN();
   coinReg = event;


   //_______________________________ RF
   grf = (grrf + gRandom->Uniform(0,1))*GR_CH2NS[0];
   brf = (lrf[1] + gRandom->Uniform(0,1))*LAS_CH2NS[0];
   
   
   //_______________________________ GR
   grdE1 = sqrt(gradc[0]*gradc[1]);
   grdE2 = sqrt(gradc[2]*gradc[3]);
   
   grT1avg = (grtdc[0]+36+grtdc[1])/2.*GR_CH2NS[0];
   grT2avg = (grtdc[2]+36+grtdc[3])/2.*GR_CH2NS[0];
   grTOF1 = grT1avg + (gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
   grTOF2 = grT2avg + (gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.*GR_CH2NS[0] - grf + GR_TOF1_OFFSET;
  
   grXC = grx - 380./2.3*grth*TMath::RadToDeg();

   //_______________________________ Gate
   
   vetogate = 0;
   if((tdc2[8]<0)&&tdc2[9]<0) vetogate += 1;
   
//   blockgate = 0;
//   stackgate = 0;
//   if((tdc2[0]>0)&&(tdc2[0]<250)&&(tdc2[1]>0)&&(tdc2[1]<250)) blockgate += 1;
//   if((tdc2[2]>0)&&(tdc2[2]<250)&&(tdc2[3]>0)&&(tdc2[3]<250)) blockgate += 2;
//   if((tdc2[4]>0)&&(tdc2[4]<250)&&(tdc2[5]>0)&&(tdc2[5]<250)) blockgate += 4;
//   if((tdc2[6]>0)&&(tdc2[6]<250)&&(tdc2[7]>0)&&(tdc2[7]<250)) blockgate += 8;

//   if((tdc[0]>0)&&(tdc[0]<400)&&(tdc[1]>0)&&(tdc[1]<400)) stackgate += 1;
//   if((tdc[2]>0)&&(tdc[2]<300)&&(tdc[3]>0)&&(tdc[3]<350)) stackgate += 2;
//   if((tdc[4]>0)&&(tdc[4]<400)&&(tdc[5]>0)&&(tdc[5]<400)) stackgate += 4;
//   if((tdc[6]>0)&&(tdc[6]<400)&&(tdc[7]>0)&&(tdc[7]<400)) stackgate += 8;
//   if((tdc[8]>0)&&(tdc[8]<400)&&(tdc[9]>0)&&(tdc[9]<400)) stackgate += 16;
//   if((tdc[10]>0)&&(tdc[10]<400)&&(tdc[11]>0)&&(tdc[11]<400)) stackgate += 32;
   
   
   //_______________________________ BAND telesope
   badEl = adc2[8];
   badEr = adc2[9];
   blo1 = sqrt(adc2[0]*adc[1]);
   blo2 = sqrt(adc2[2]*adc[3]);
   blo3 = sqrt(adc2[4]*adc[5]);
   blo4 = sqrt(adc2[6]*adc[7]);
   blo1Tavg = (tdc2[0]+tdc2[1])/2.*LAS_CH2NS[0];
   blo2Tavg = (tdc2[2]+tdc2[3])/2.*LAS_CH2NS[1];
   blo3Tavg = (tdc2[4]+tdc2[5])/2.*LAS_CH2NS[2];
   blo4Tavg = (tdc2[6]+tdc2[7])/2.*LAS_CH2NS[3];
   
//   if (tdc2[8]>0) badElTOF = tdc2[8]*LAS_CH2NS-brf+BLOCK_TOF_OFFSET;
//   if (tdc2[9]>0) badErTOF = tdc2[9]*LAS_CH2NS-brf+BLOCK_TOF_OFFSET;
//   if (blockgate&1) blo1TOF = blo1Tavg-brf+BLOCK_TOF_OFFSET;
//   if (blockgate&2) blo2TOF = blo2Tavg-brf+BLOCK_TOF_OFFSET;
//   if (blockgate&4) blo3TOF = blo3Tavg-brf+BLOCK_TOF_OFFSET;
//   if (blockgate&8) blo4TOF = blo4Tavg-brf+BLOCK_TOF_OFFSET;
   
   blo1TOF = TMath::QuietNaN();
   blo2TOF = TMath::QuietNaN();
   blo3TOF = TMath::QuietNaN();
   blo4TOF = TMath::QuietNaN();
   
   if( tdc2[8] != -1 ) badElTOF = (tdc2[8] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
   if( tdc2[9] != -1 ) badErTOF = (tdc2[9] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
   
//   if(blo1Tavg != -LAS_CH2NS[0] )blo1TOF = blo1Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[0] + BLOCK_TOF_OFFSET[0];
//   if(blo2Tavg != -LAS_CH2NS[1] )blo2TOF = blo2Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[1] + BLOCK_TOF_OFFSET[1];
//   if(blo3Tavg != -LAS_CH2NS[2] )blo3TOF = blo3Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[2] + BLOCK_TOF_OFFSET[2];
//   if(blo4Tavg != -LAS_CH2NS[3] )blo4TOF = blo4Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[3] + BLOCK_TOF_OFFSET[3];
   
   
   if(blo1Tavg != -LAS_CH2NS[0] )blo1TOF = blo1Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[0] ;//+ BLOCK_TOF_OFFSET[0];
   if(blo2Tavg != -LAS_CH2NS[1] )blo2TOF = blo2Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[1] ;//+ BLOCK_TOF_OFFSET[1];
   if(blo3Tavg != -LAS_CH2NS[2] )blo3TOF = blo3Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[2] ;//+ BLOCK_TOF_OFFSET[2];
   if(blo4Tavg != -LAS_CH2NS[3] )blo4TOF = blo4Tavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[3] ;//+ BLOCK_TOF_OFFSET[3];
   
   
   //_______________________________ BAND stack
//   if (adc[0]>-5 && adc[1]>-5) sta1h = sqrt(adc[0]/112.*adc[1]/118.);
//   if (adc[2]>-5 && adc[3]>-5) sta2h = sqrt(adc[2]/131.*adc[3]/109.);
//   if (adc[4]>-5 && adc[5]>-5) sta1v = sqrt(adc[4]/116.*adc[5]/126.);
//   if (adc[6]>-5 && adc[7]>-5) sta2v = sqrt(adc[6]/113.*adc[7]/123.);
//   if (adc[8]>-5 && adc[9]>-5) sta3v = sqrt(adc[8]/122.*adc[9]/131.);
//   if (adc[10]>-5 && adc[11]>-5) sta4v = sqrt(adc[10]/127.*adc[11]/124.);
   
//   sta1h = sqrt(TMath::Max(adc[0] ,0)/112.*TMath::Max(adc[1] ,0)/118.);
//   sta2h = sqrt(TMath::Max(adc[2] ,0)/131.*TMath::Max(adc[3] ,0)/109.);
//   sta1v = sqrt(TMath::Max(adc[4] ,0)/116.*TMath::Max(adc[5] ,0)/126.);
//   sta2v = sqrt(TMath::Max(adc[6] ,0)/113.*TMath::Max(adc[7] ,0)/123.);
//   sta3v = sqrt(TMath::Max(adc[8] ,0)/122.*TMath::Max(adc[9] ,0)/131.);
//   sta4v = sqrt(TMath::Max(adc[10],0)/127.*TMath::Max(adc[11],0)/124.);


   sta1h = sqrt(adc[0]/112.*adc[1]/118.);
   sta2h = sqrt(adc[2]/131.*adc[3]/109.);
   sta1v = sqrt(adc[4]/116.*adc[5]/126.);
   sta2v = sqrt(adc[6]/113.*adc[7]/123.);
   sta3v = sqrt(adc[8]/122.*adc[9]/131.);
   sta4v = sqrt(adc[10]/127.*adc[11]/124.);
   
   sta_even = 6.4*(sta1h+sta2h);
   sta_odd = 6.4*(sta1v+sta2v+sta3v+sta4v);
//   if (sta_even+sta_odd>0.&& stackgate) sta_ratio = (sta_even-sta_odd)/(sta_even+sta_odd);
   if (sta_even+sta_odd>0.) {
        sta_sum = sta_even+sta_odd;
        sta_ratio = (sta_even-sta_odd)/(sta_even+sta_odd);
   }
   sta1hTavg = (tdc[0]+tdc[1])/2.*LAS_CH2NS[4];
   sta2hTavg = (tdc[2]+tdc[3])/2.*LAS_CH2NS[5];
   sta1vTavg = (tdc[4]+tdc[5])/2.*LAS_CH2NS[6];
   sta2vTavg = (tdc[6]+tdc[7])/2.*LAS_CH2NS[7];
   sta3vTavg = (tdc[8]+tdc[9])/2.*LAS_CH2NS[8];
   sta4vTavg = (tdc[10]+tdc[11])/2.*LAS_CH2NS[9];
   //if((tdc[0]>0)&&(tdc[0]<60)&&(tdc[1]>0)&&(tdc[1]<60)) sta1hTOF = sta1hTavg-brf+STACK_TOF_OFFSET;
   //if((tdc[2]>0)&&(tdc[2]<60)&&(tdc[3]>0)&&(tdc[3]<60)) sta2hTOF = sta2hTavg-brf+STACK_TOF_OFFSET;
   //if((tdc[4]>0)&&(tdc[4]<100)&&(tdc[5]>0)&&(tdc[5]<100)) sta1vTOF = sta1vTavg-brf+STACK_TOF_OFFSET;
   //if((tdc[6]>0)&&(tdc[6]<120)&&(tdc[7]>0)&&(tdc[7]<120)) sta2vTOF = sta2vTavg-brf+STACK_TOF_OFFSET;
   //if((tdc[8]>0)&&(tdc[8]<110)&&(tdc[9]>0)&&(tdc[9]<110)) sta3vTOF = sta3vTavg-brf+STACK_TOF_OFFSET;
   //if((tdc[10]>0)&&(tdc[10]<120)&&(tdc[11]>0)&&(tdc[11]<120)) sta4vTOF = sta4vTavg-brf+STACK_TOF_OFFSET;
   
//   if(sta1hTavg != -LAS_CH2NS[4] ) sta1hTOF = sta1hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[4] + STACK_TOF_OFFSET[0];
//   if(sta2hTavg < 700*LAS_CH2NS[5] ) sta2hTOF = sta2hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[5] + STACK_TOF_OFFSET[1];
//   if(sta1vTavg != -LAS_CH2NS[6] ) sta1vTOF = sta1vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[6] + STACK_TOF_OFFSET[2];
//   if(sta2vTavg != -LAS_CH2NS[7] ) sta2vTOF = sta2vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[7] + STACK_TOF_OFFSET[3];
//   if(sta3vTavg != -LAS_CH2NS[8] ) sta3vTOF = sta3vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[8] + STACK_TOF_OFFSET[4];
//   if(sta4vTavg != -LAS_CH2NS[9] ) sta4vTOF = sta4vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2. - (lrf[1] + gRandom->Uniform(0,1)))*LAS_CH2NS[9] + STACK_TOF_OFFSET[5];
   
   if(sta1hTavg != -LAS_CH2NS[4]   ) sta1hTOF = sta1hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[4] ;//+ STACK_TOF_OFFSET[0];
   if(sta2hTavg < 700*LAS_CH2NS[5] ) sta2hTOF = sta2hTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[5] ;//+ STACK_TOF_OFFSET[1];
   if(sta1vTavg != -LAS_CH2NS[6]   ) sta1vTOF = sta1vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[6] ;//+ STACK_TOF_OFFSET[2];
   if(sta2vTavg != -LAS_CH2NS[7]   ) sta2vTOF = sta2vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[7] ;//+ STACK_TOF_OFFSET[3];
   if(sta3vTavg != -LAS_CH2NS[8]   ) sta3vTOF = sta3vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[8] ;//+ STACK_TOF_OFFSET[4];
   if(sta4vTavg != -LAS_CH2NS[9]   ) sta4vTOF = sta4vTavg + ((gRandom->Uniform(0,1)+gRandom->Uniform(0,1))/2.)*LAS_CH2NS[9] ;//+ STACK_TOF_OFFSET[5];
   
   //_______________________________ BAND liquid
   liqlf = adc2[10];
   liqld = adc2[11];
   liqrf = adc2[12];
   liqrd = adc2[13];
   liqlTOF = (tdc2[10] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;
   liqrTOF = (tdc2[11] + gRandom->Uniform(0,1))*LAS_CH2NS[0] - brf ;

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
   //_______________________________ Filter Phsyical data
    
   
   /**/
   count ++;
   
   ///______________________________________________________________   
   saveFile->cd(); //set focus on this file
   newTree->Fill();  

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

  
   return kTRUE;
}

void Selector_Sort::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Selector_Sort::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   saveFile->cd(); //set focus on this file
   newTree->Write(); 
   saveFile->Close();

   printf("-------------- done. %s, %d\n", saveFileName.Data(), count);

}

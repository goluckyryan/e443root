#define Selector_disc_cxx
// The class definition in Selector_disc.h has been generated automatically
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
// Root > T->Process("Selector_disc.C")
// Root > T->Process("Selector_disc.C","some options")
// Root > T->Process("Selector_disc.C+")
//

#include "Selector_disc.h"
#include <TH2.h>
#include <TStyle.h>


void Selector_disc::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void Selector_disc::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Selector_disc::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Selector_disc::GetEntry() or TBranch::GetEntry()
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


  //_____________________________
  b_coinReg->GetEntry(entry);
  b_vetogate->GetEntry(entry); 
  b_grx->GetEntry(entry);      
  b_gry->GetEntry(entry);
  b_grth->GetEntry(entry);
  b_grph->GetEntry(entry);
  b_grdE1->GetEntry(entry);
  b_grTOF1->GetEntry(entry);
  b_grXC->GetEntry(entry);
  b_grthC->GetEntry(entry);
  b_grXAux->GetEntry(entry);
  b_liqlf->GetEntry(entry);
  b_liqld->GetEntry(entry);
  b_liqrf->GetEntry(entry);
  b_liqrd->GetEntry(entry);
  b_liqlTOF->GetEntry(entry);
  b_liqrTOF->GetEntry(entry);


  //______________________________

  ratio1 = Findratio(1, liqld, liqlf);
  ratio2 = Findratio(2, liqrd, liqrf);

      
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

Double_t Selector_disc::Findratio(Int_t id, Double_t Xpos, Double_t Ypos)
{
  Double_t ratio = TMath::QuietNaN();
  
  if( id == 1){ // Left

    if( 6 < Xpos && Xpos <= 40){
      Double_t m1 = polL1->Derivative(Xpos); 
      Double_t m2 = polL2->Derivative(Xpos);

      Double_t c1 = polL1->Eval(Xpos)-m1*Xpos; 
      Double_t c2 = polL2->Eval(Xpos)-m2*Xpos;

      Double_t xp = -(c1-c2)/(m1-m2);
      Double_t yp = m1*xp + c1;

      Double_t slop1 = (Ypos-yp)/(Xpos-xp);
      Double_t slop0 = (polL2->Eval(Xpos)-yp)/(Xpos-xp);

      ratio = (slop1/slop0 -1)/polLN(Xpos) +1;
      
    }else if( Xpos > 40 ){

      ratio = (Ypos+213.16)/(Xpos+228.1)/1.11;

    }

  }else{ // Right

    if( 6 < Xpos && Xpos <= 40){

      Double_t m1 = polR1->Derivative(Xpos); 
      Double_t m2 = polR2->Derivative(Xpos);

      Double_t c1 = polR1->Eval(Xpos)-m1*Xpos; 
      Double_t c2 = polR2->Eval(Xpos)-m2*Xpos;

      Double_t xp = -(c1-c2)/(m1-m2);
      Double_t yp = m1*xp + c1;

      Double_t slop1 = (liqlf-yp)/(Xpos-xp);
      Double_t slop0 = (polR2->Eval(Xpos)-yp)/(Xpos-xp);

      ratio = (slop1/slop0 -1)/polRN(Xpos) +1;
      
    }else if( Xpos > 40 ){

      ratio = (Ypos+197.77)/(Xpos+216.14)*138./155.;
    }

  }
  return ratio;
}

void Selector_disc::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Selector_disc::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

   saveFile->cd(); //set focus on this file
   newTree->Write(); 
   saveFile->Close();

   printf("-------------- done. %s, %d\n", saveFileName.Data(), count);
}

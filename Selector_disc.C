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

  //printf("%d, Start Finding distance\n", entry);

  FindDistance(1, liqld, liqlf);
  
  dist1X = result[0];
  dist1Y = result[1];
  dist1  = result[2];

  //printf("%d, (%f,%f), dist1X : %f, dist : %f \n", entry, liqld, liqlf, dist1X, dist1);
  
  FindDistance(2, liqrd, liqrf); 
 
  dist2X = result[0];
  dist2Y = result[1];
  dist2  = result[2];
      
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
               (Int_t) TMath::Nint((entry+1)*100./totnumEntry),
               (Int_t) TMath::Floor(time/60), time - TMath::Floor(time/60)*60,
               totnumEntry*time/(entry+1)/60.,
               (int)count);
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

  //TF1 * pol1;  for simpler code
  //TF1 * pol1;
  //TF1 * polN;
  // 
  ////printf("FR=========================== \n");
  //if( id == 1){
  //  pol1 = polL1;
  //  pol2 = polL2;
  //  polN = polLN;
  //}else{
  //  pol1 = polR1;
  //  pol2 = polR2;
  //  polN = polRN;
  //}
  

  if( id == 1){ // Left

    if( 6 < Xpos && Xpos <= 40){

      if( Ypos > 3.75 * Xpos + 25) return ratio;

      Double_t m1 = polL1->Derivative(Xpos); 
      Double_t m2 = polL2->Derivative(Xpos);

      Double_t c1 = polL1->Eval(Xpos)-m1*Xpos; 
      Double_t c2 = polL2->Eval(Xpos)-m2*Xpos;

      Double_t xp = -(c1-c2)/(m1-m2);
      Double_t yp = m1*xp + c1;

      Double_t slop1 = (Ypos-yp)/(Xpos-xp);
      Double_t slop0 = (polL2->Eval(Xpos)-yp)/(Xpos-xp);

      ratio = (slop1/slop0 -1)/polLN->Eval(Xpos) +1;
      
    }else if( Xpos > 40 ){

      if( Ypos > 5./3. * Xpos + 200./3.) return ratio;

      ratio = (Ypos+213.16)/(Xpos+228.1)/1.11;

    }

  }else{ // Right

    if( 6 < Xpos && Xpos <= 40){

      if( Ypos > 4. * Xpos + 20) return ratio;

      Double_t m1 = polR1->Derivative(Xpos); 
      Double_t m2 = polR2->Derivative(Xpos);

      Double_t c1 = polR1->Eval(Xpos)-m1*Xpos; 
      Double_t c2 = polR2->Eval(Xpos)-m2*Xpos;

      Double_t xp = -(c1-c2)/(m1-m2);
      Double_t yp = m1*xp + c1;

      Double_t slop1 = (liqlf-yp)/(Xpos-xp);
      Double_t slop0 = (polR2->Eval(Xpos)-yp)/(Xpos-xp);

      ratio = (slop1/slop0 -1)/polRN->Eval(Xpos) +1;
      
    }else if( Xpos > 40 ){

      if( Ypos > 5./3. * Xpos + 200./3.) return ratio;

      ratio = (Ypos+197.77)/(Xpos+216.14)*138./155.;
    }

  }
  return ratio;
}

void Selector_disc::FindDistance(Int_t id, Double_t Xpos, Double_t Ypos){
  Double_t dist = TMath::QuietNaN();
  Int_t searchDir = 0;
  //TProfile * htemp;
  TF1 * gtemp;
  TF1 * ftemp;

  //printf("FD=========================== \n");
  if( id == 1){
    gtemp = g1;
    ftemp = f1;
  }else{
    gtemp = g2;
    ftemp = f2;
  }

  Double_t tempY  ;
  Double_t tempDX ;

  if( Xpos <= xCut){

    //tempY  = htemp->Interpolate(Xpos); 
    //tempDX = htemp->Interpolate(Xpos+1)- htempY;

    tempY  = gtemp->Eval(Xpos); 
    tempDX = gtemp->Derivative(Xpos);


  }else{

    tempY  = ftemp->Eval(Xpos); 
    tempDX = ftemp->Derivative(Xpos);

  }
  
  if( tempY == Ypos){
    dist = 0.0;
  }else{
    if(tempDX == 0 ){

      dist = TMath::Abs(Ypos - tempY);

    }else if(tempDX > 0 ){

      if( Ypos > tempY){
        searchDir = 1; // forward
      }else{
        searchDir = -1;
      }

    }else{

      if( Ypos > tempY){
        searchDir = -1; // forward
      }else{
        searchDir = 1;
      }

    }
  }

  Double_t xScan = Xpos;
  Double_t yScan = TMath::QuietNaN();

  if( searchDir != 0){
    Double_t step =0.01;
    Double_t tempDist = 100000000.;
    Int_t counti = 0;
    do{
      dist = tempDist;
      xScan += searchDir*step;
      if( xScan <= xCut){
        //yScan = htemp->Interpolate(xScan);
        yScan = gtemp->Eval(xScan);
      }else{
        yScan = ftemp->Eval(xScan);
      }

      tempDist = TMath::Sqrt(TMath::Power(Xpos-xScan,2)+TMath::Power(Ypos-yScan,2));
      counti ++;
      //printf("%d, %4f, dTemp:%f, dist:%f \n",counti,  xScan, tempDist, dist);
    }while (dist > tempDist);
    xScan -= searchDir*step;
  }
  //printf("============= xScan:%f, dist:%f \n", xScan, dist);

  result[0] = xScan;
  //result[1] = yScan;
  if( xScan <= xCut){
    //result[1] = htemp->Interpolate(xScan);
    result[1] = gtemp->Eval(xScan);
  }else{
    result[1] = ftemp->Eval(xScan);
  }
  result[2] = TMath::Sign(dist, Ypos-tempY);

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

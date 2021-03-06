//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec  8 20:43:16 2015 by ROOT version 5.34/32
// from TTree t1/e443 events
// found on file: dstroot/run1035.root
//////////////////////////////////////////////////////////

#ifndef Selector_Sort_h
#define Selector_Sort_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "constant.h"

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Selector_Sort : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   TBenchmark clock;
   Bool_t shown;

   Int_t count;

   //=============================== Store in new ROOT file
   TFile * saveFile;
   TTree * newTree;
   
   TString saveFileName;
   
   Int_t totnumEntry;

   //---------eventID
   Int_t eventID;
   //---------coinReg
   Int_t coinReg;
   //--------- gate
   Int_t gate;
   Int_t vetogate, blockgate,stackgate;
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
   
   //------------------------------------------- Gate Detail
   //TCutG * pid3He;
   Int_t vetogate, blockgate, stackgate;
   
   //=============================== Declaration of leaf types
   Int_t           event;
   Double_t        grx;
   Double_t        grth;
   Double_t        gry;
   Double_t        grph;
   Int_t           gradc[4];
   Int_t           grtdc[4];
   Int_t           grrf;
   Int_t           adc[12];
   Int_t           adc2[14];
   Int_t           tdc[12];
   Int_t           tdc2[12];
   Int_t           lrf[3];

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_grx;   //!
   TBranch        *b_grth;   //!
   TBranch        *b_gry;   //!
   TBranch        *b_grph;   //!
   TBranch        *b_gradc;   //!
   TBranch        *b_grtdc;   //!
   TBranch        *b_grrf;   //!
   TBranch        *b_adc;   //!
   TBranch        *b_adc2;   //!
   TBranch        *b_tdc;   //!
   TBranch        *b_tdc2;   //!
   TBranch        *b_lrf;   //!

   Selector_Sort(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Selector_Sort() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(Selector_Sort,0);
};

#endif

#ifdef Selector_Sort_cxx
void Selector_Sort::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   totnumEntry = tree->GetEntries();
  
   //======================= Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("grx", &grx, &b_grx);
   fChain->SetBranchAddress("grth", &grth, &b_grth);
   fChain->SetBranchAddress("gry", &gry, &b_gry);
   fChain->SetBranchAddress("grph", &grph, &b_grph);
   fChain->SetBranchAddress("gradc", gradc, &b_gradc);
   fChain->SetBranchAddress("grtdc", grtdc, &b_grtdc);
   fChain->SetBranchAddress("grrf", &grrf, &b_grrf);
   fChain->SetBranchAddress("adc", adc, &b_adc);
   fChain->SetBranchAddress("adc2", adc2, &b_adc2);
   fChain->SetBranchAddress("tdc", tdc, &b_tdc);
   fChain->SetBranchAddress("tdc2", tdc2, &b_tdc2);
   fChain->SetBranchAddress("lrf", lrf, &b_lrf);

   //================================Store in New ROOT file
   saveFileName = fChain->GetDirectory()->GetName();
   saveFileName = "S_"+saveFileName;
   
   totnumEntry = tree->GetEntries();
   
   printf("Converting %s ------> %s , total Entry : %d \n", fChain->GetDirectory()->GetName(), saveFileName.Data(), totnumEntry);
   
   saveFile = new TFile( saveFileName,"recreate");
   newTree =  new TTree("tree","tree");
   
   clock.Reset();
   clock.Start("timer");
   shown = 0;

  
   //================================= Initialize Tree variable
   eventID = 0;
   //---------coinReg
   coinReg = -1;
   //--------- gate
   gate = -1;
   vetogate = -1;
   blockgate = -1;
   stackgate = -1;
   //--------- GR
   //GRX = TMath::QuietNaN();
   //GRY = TMath::QuietNaN();
   //GRTH = TMath::QuietNaN();
   //GRPH = TMath::QuietNaN();

   grdE1 = TMath::QuietNaN();
   grdE2 = TMath::QuietNaN();
   grT1avg = TMath::QuietNaN();
   grT2avg = TMath::QuietNaN();
   grTOF1 = TMath::QuietNaN();
   grTOF2 = TMath::QuietNaN();

   grXC = TMath::QuietNaN();

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
 
 
   //=================================== Declare Leaf of Tree
   newTree->Branch("eventID", &eventID, "eventID/I");
   newTree->Branch("coinReg", &coinReg, "coinReg/I");
   newTree->Branch("gate", &gate, "gate/I");
   newTree->Branch("grx", &grx, "grx/D");
   newTree->Branch("gry", &gry, "gry/D");
   newTree->Branch("grth", &grth, "grth/D");
   newTree->Branch("grph", &grph, "grph/D");

   newTree->Branch("grdE1", &grdE1, "grdE1/D");
   newTree->Branch("grdE2", &grdE2, "grdE2/D");
   newTree->Branch("grTOF1", &grTOF1, "grTOF1/D");
   newTree->Branch("grTOF2", &grTOF2, "grTOF2/D");
   
   newTree->Branch("grf", &grf, "grf/D");
   newTree->Branch("grXC", &grXC, "grXC/D");

   newTree->Branch("badEl", &badEl, "badEl/D");
   newTree->Branch("badEr", &badEr, "badEr/D");
   newTree->Branch("vetogate", &vetogate, "vetogate/I");
   newTree->Branch("blo1", &blo1, "blo1/D");
   newTree->Branch("blo2", &blo2, "blo2/D");
   newTree->Branch("blo3", &blo3, "blo3/D");
   newTree->Branch("blo4", &blo4, "blo4/D");
   newTree->Branch("badElTOF", &badElTOF, "badElTOF/D");
   newTree->Branch("badErTOF", &badErTOF, "badErTOF/D");
   newTree->Branch("blo1Tavg", &blo1Tavg, "blo1Tavg/D");
   newTree->Branch("blo2Tavg", &blo2Tavg, "blo2Tavg/D");
   newTree->Branch("blo3Tavg", &blo3Tavg, "blo3Tavg/D");
   newTree->Branch("blo4Tavg", &blo4Tavg, "blo4Tavg/D");
   newTree->Branch("blo1TOF", &blo1TOF, "blo1TOF/D");
   newTree->Branch("blo2TOF", &blo2TOF, "blo2TOF/D");
   newTree->Branch("blo3TOF", &blo3TOF, "blo3TOF/D");
   newTree->Branch("blo4TOF", &blo4TOF, "blo4TOF/D");

   newTree->Branch("sta_odd", &sta_odd, "sta_odd/D");
   newTree->Branch("sta_even", &sta_even, "sta_even/D");
   newTree->Branch("sta_ratio", &sta_ratio, "sta_ratio/D");
   newTree->Branch("sta1hTavg", &sta1hTavg, "sta1hTavg/D");
   newTree->Branch("sta2hTavg", &sta2hTavg, "sta2hTavg/D");
   newTree->Branch("sta1vTavg", &sta1vTavg, "sta1vTavg/D");
   newTree->Branch("sta2vTavg", &sta2vTavg, "sta2vTavg/D");
   newTree->Branch("sta3vTavg", &sta3vTavg, "sta3vTavg/D");
   newTree->Branch("sta4vTavg", &sta4vTavg, "sta4vTavg/D");
   newTree->Branch("sta1hTOF", &sta1hTOF, "sta1hTOF/D");
   newTree->Branch("sta2hTOF", &sta2hTOF, "sta2hTOF/D");
   newTree->Branch("sta1vTOF", &sta1vTOF, "sta1vTOF/D");
   newTree->Branch("sta2vTOF", &sta2vTOF, "sta2vTOF/D");
   newTree->Branch("sta3vTOF", &sta3vTOF, "sta3vTOF/D");
   newTree->Branch("sta4vTOF", &sta4vTOF, "sta4vTOF/D");

   newTree->Branch("liqlf", &liqlf, "liqlf/D");
   newTree->Branch("liqld", &liqld, "liqld/D");
   newTree->Branch("liqrf", &liqrf, "liqrf/D");
   newTree->Branch("liqrd", &liqrd, "liqrd/D");
   newTree->Branch("liqlTOF", &liqlTOF, "liqlTOF/D");
   newTree->Branch("liqrTOF", &liqrTOF, "liqrTOF/D");

   newTree->Branch("brf", &brf, "brf/D");
   
}

Bool_t Selector_Sort::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Selector_Sort_cxx

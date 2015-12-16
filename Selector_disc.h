//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 16 18:13:35 2015 by ROOT version 5.34/32
// from TTree tree/tree
// found on file: run1035.root
//////////////////////////////////////////////////////////

#ifndef Selector_disc_h
#define Selector_disc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class Selector_disc : public TSelector {
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

   TF1 * polL1;
   TF1 * polL2;
   TF1 * polLN;
   
   TF1 * polR1;
   TF1 * polR2;
   TF1 * polRN;

   //======================== newTree variable
   //   Double_t liqlf,liqld,liqrf,liqrd;
   //   Double_t liqlTOF, liqrTOF;

   //Double_t xp,yp;
   //Double_t slop1, slop0;
   Double_t ratio1, ratio2;

   // Declaration of leaf types
  // Int_t           eventID;
   Int_t           coinReg;
   Int_t           vetogate;
   Double_t        grx;
   Double_t        gry;
   Double_t        grth;
   Double_t        grph;
   Double_t        grdE1;
   // Double_t        grdE2;
   Double_t        grTOF1;
   // Double_t        grTOF2;
  // Double_t        grf;
   Double_t        grXC;
   Double_t        grthC;
   Double_t        grXAux;
  // Double_t        badEl;
  // Double_t        badEr;
  // Double_t        blo1;
  // Double_t        blo2;
  // Double_t        blo3;
  // Double_t        blo4;
  // Double_t        badElTOF;
  // Double_t        badErTOF;
  // Double_t        blo1Tavg;
  // Double_t        blo2Tavg;
  // Double_t        blo3Tavg;
  // Double_t        blo4Tavg;
  // Double_t        blo1TOF;
  // Double_t        blo2TOF;
  // Double_t        blo3TOF;
  // Double_t        blo4TOF;
  // Double_t        adc[12];
  // Double_t        tdc[12];
  // Int_t           staMult;
  // Double_t        sta1h;
  // Double_t        sta2h;
  // Double_t        sta1v;
  // Double_t        sta2v;
  // Double_t        sta3v;
  // Double_t        sta4v;
  // Double_t        sta_odd;
  // Double_t        sta_even;
  // Double_t        sta_sum;
  // Double_t        sta_ratio;
  // Double_t        sta1hTavg;
  // Double_t        sta2hTavg;
  // Double_t        sta1vTavg;
  // Double_t        sta2vTavg;
  // Double_t        sta3vTavg;
  // Double_t        sta4vTavg;
  // Double_t        sta1hTdif;
  // Double_t        sta2hTdif;
  // Double_t        sta1vTdif;
  // Double_t        sta2vTdif;
  // Double_t        sta3vTdif;
  // Double_t        sta4vTdif;
  // Double_t        sta1hTOF;
  // Double_t        sta2hTOF;
  // Double_t        sta1vTOF;
  // Double_t        sta2vTOF;
  // Double_t        sta3vTOF;
  // Double_t        sta4vTOF;
   Double_t        liqlf;
   Double_t        liqld;
   Double_t        liqrf;
   Double_t        liqrd;
   Double_t        liqlTOF;
   Double_t        liqrTOF;
  // Double_t        brf;

   // List of branches
   // TBranch        *b_eventID;   //!
    TBranch        *b_coinReg;   //!
    TBranch        *b_vetogate;   //!
    TBranch        *b_grx;   //!
    TBranch        *b_gry;   //!
    TBranch        *b_grth;   //!
    TBranch        *b_grph;   //!
    TBranch        *b_grdE1;   //!
   // TBranch        *b_grdE2;   //!
    TBranch        *b_grTOF1;   //!
   // TBranch        *b_grTOF2;   //!
   // TBranch        *b_grf;   //!
    TBranch        *b_grXC;   //!
    TBranch        *b_grthC;   //!
    TBranch        *b_grXAux;   //!
   // TBranch        *b_badEl;   //!
   // TBranch        *b_badEr;   //!
   // TBranch        *b_blo1;   //!
   // TBranch        *b_blo2;   //!
   // TBranch        *b_blo3;   //!
   // TBranch        *b_blo4;   //!
   // TBranch        *b_badElTOF;   //!
   // TBranch        *b_badErTOF;   //!
   // TBranch        *b_blo1Tavg;   //!
   // TBranch        *b_blo2Tavg;   //!
   // TBranch        *b_blo3Tavg;   //!
   // TBranch        *b_blo4Tavg;   //!
   // TBranch        *b_blo1TOF;   //!
   // TBranch        *b_blo2TOF;   //!
   // TBranch        *b_blo3TOF;   //!
   // TBranch        *b_blo4TOF;   //!
   // TBranch        *b_adc;   //!
   // TBranch        *b_tdc;   //!
   // TBranch        *b_staMult;   //!
   // TBranch        *b_sta1h;   //!
   // TBranch        *b_sta2h;   //!
   // TBranch        *b_sta1v;   //!
   // TBranch        *b_sta2v;   //!
   // TBranch        *b_sta3v;   //!
   // TBranch        *b_sta4v;   //!
   // TBranch        *b_sta_odd;   //!
   // TBranch        *b_sta_even;   //!
   // TBranch        *b_sta_sum;   //!
   // TBranch        *b_sta_ratio;   //!
   // TBranch        *b_sta1hTavg;   //!
   // TBranch        *b_sta2hTavg;   //!
   // TBranch        *b_sta1vTavg;   //!
   // TBranch        *b_sta2vTavg;   //!
   // TBranch        *b_sta3vTavg;   //!
   // TBranch        *b_sta4vTavg;   //!
   // TBranch        *b_sta1hTdif;   //!
   // TBranch        *b_sta2hTdif;   //!
   // TBranch        *b_sta1vTdif;   //!
   // TBranch        *b_sta2vTdif;   //!
   // TBranch        *b_sta3vTdif;   //!
   // TBranch        *b_sta4vTdif;   //!
   // TBranch        *b_sta1hTOF;   //!
   // TBranch        *b_sta2hTOF;   //!
   // TBranch        *b_sta1vTOF;   //!
   // TBranch        *b_sta2vTOF;   //!
   // TBranch        *b_sta3vTOF;   //!
   // TBranch        *b_sta4vTOF;   //!
   TBranch        *b_liqlf;   //!
   TBranch        *b_liqld;   //!
   TBranch        *b_liqrf;   //!
   TBranch        *b_liqrd;   //!
   TBranch        *b_liqlTOF;   //!
   TBranch        *b_liqrTOF;   //!
   //TBranch        *b_brf;   //!

   Selector_disc(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Selector_disc() { }
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
   virtual Double_t  Findratio(Int_t id, Double_t Xpos, Double_t Ypos);

   ClassDef(Selector_disc,0);
};

#endif

#ifdef Selector_disc_cxx
void Selector_disc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

  // fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("coinReg", &coinReg, &b_coinReg);
   fChain->SetBranchAddress("vetogate", &vetogate, &b_vetogate);
   fChain->SetBranchAddress("grx", &grx, &b_grx);
   fChain->SetBranchAddress("gry", &gry, &b_gry);
   fChain->SetBranchAddress("grth", &grth, &b_grth);
   fChain->SetBranchAddress("grph", &grph, &b_grph);
   fChain->SetBranchAddress("grdE1", &grdE1, &b_grdE1);
   // fChain->SetBranchAddress("grdE2", &grdE2, &b_grdE2);
   fChain->SetBranchAddress("grTOF1", &grTOF1, &b_grTOF1);
   //fChain->SetBranchAddress("grTOF2", &grTOF2, &b_grTOF2);
   // fChain->SetBranchAddress("grf", &grf, &b_grf);
   fChain->SetBranchAddress("grXC", &grXC, &b_grXC);
   fChain->SetBranchAddress("grthC", &grthC, &b_grthC);
   fChain->SetBranchAddress("grXAux", &grXAux, &b_grXAux);
  // fChain->SetBranchAddress("badEl", &badEl, &b_badEl);
  // fChain->SetBranchAddress("badEr", &badEr, &b_badEr);
  // fChain->SetBranchAddress("blo1", &blo1, &b_blo1);
  // fChain->SetBranchAddress("blo2", &blo2, &b_blo2);
  // fChain->SetBranchAddress("blo3", &blo3, &b_blo3);
  // fChain->SetBranchAddress("blo4", &blo4, &b_blo4);
  // fChain->SetBranchAddress("badElTOF", &badElTOF, &b_badElTOF);
  // fChain->SetBranchAddress("badErTOF", &badErTOF, &b_badErTOF);
  // fChain->SetBranchAddress("blo1Tavg", &blo1Tavg, &b_blo1Tavg);
  // fChain->SetBranchAddress("blo2Tavg", &blo2Tavg, &b_blo2Tavg);
  // fChain->SetBranchAddress("blo3Tavg", &blo3Tavg, &b_blo3Tavg);
  // fChain->SetBranchAddress("blo4Tavg", &blo4Tavg, &b_blo4Tavg);
  // fChain->SetBranchAddress("blo1TOF", &blo1TOF, &b_blo1TOF);
  // fChain->SetBranchAddress("blo2TOF", &blo2TOF, &b_blo2TOF);
  // fChain->SetBranchAddress("blo3TOF", &blo3TOF, &b_blo3TOF);
  // fChain->SetBranchAddress("blo4TOF", &blo4TOF, &b_blo4TOF);
  // fChain->SetBranchAddress("adc", adc, &b_adc);
  // fChain->SetBranchAddress("tdc", tdc, &b_tdc);
  // fChain->SetBranchAddress("staMult", &staMult, &b_staMult);
  // fChain->SetBranchAddress("sta1h", &sta1h, &b_sta1h);
  // fChain->SetBranchAddress("sta2h", &sta2h, &b_sta2h);
  // fChain->SetBranchAddress("sta1v", &sta1v, &b_sta1v);
  // fChain->SetBranchAddress("sta2v", &sta2v, &b_sta2v);
  // fChain->SetBranchAddress("sta3v", &sta3v, &b_sta3v);
  // fChain->SetBranchAddress("sta4v", &sta4v, &b_sta4v);
  // fChain->SetBranchAddress("sta_odd", &sta_odd, &b_sta_odd);
  // fChain->SetBranchAddress("sta_even", &sta_even, &b_sta_even);
  // fChain->SetBranchAddress("sta_sum", &sta_sum, &b_sta_sum);
  // fChain->SetBranchAddress("sta_ratio", &sta_ratio, &b_sta_ratio);
  // fChain->SetBranchAddress("sta1hTavg", &sta1hTavg, &b_sta1hTavg);
  // fChain->SetBranchAddress("sta2hTavg", &sta2hTavg, &b_sta2hTavg);
  // fChain->SetBranchAddress("sta1vTavg", &sta1vTavg, &b_sta1vTavg);
  // fChain->SetBranchAddress("sta2vTavg", &sta2vTavg, &b_sta2vTavg);
  // fChain->SetBranchAddress("sta3vTavg", &sta3vTavg, &b_sta3vTavg);
  // fChain->SetBranchAddress("sta4vTavg", &sta4vTavg, &b_sta4vTavg);
  // fChain->SetBranchAddress("sta1hTdif", &sta1hTdif, &b_sta1hTdif);
  // fChain->SetBranchAddress("sta2hTdif", &sta2hTdif, &b_sta2hTdif);
  // fChain->SetBranchAddress("sta1vTdif", &sta1vTdif, &b_sta1vTdif);
  // fChain->SetBranchAddress("sta2vTdif", &sta2vTdif, &b_sta2vTdif);
  // fChain->SetBranchAddress("sta3vTdif", &sta3vTdif, &b_sta3vTdif);
  // fChain->SetBranchAddress("sta4vTdif", &sta4vTdif, &b_sta4vTdif);
  // fChain->SetBranchAddress("sta1hTOF", &sta1hTOF, &b_sta1hTOF);
  // fChain->SetBranchAddress("sta2hTOF", &sta2hTOF, &b_sta2hTOF);
  // fChain->SetBranchAddress("sta1vTOF", &sta1vTOF, &b_sta1vTOF);
  // fChain->SetBranchAddress("sta2vTOF", &sta2vTOF, &b_sta2vTOF);
  // fChain->SetBranchAddress("sta3vTOF", &sta3vTOF, &b_sta3vTOF);
  // fChain->SetBranchAddress("sta4vTOF", &sta4vTOF, &b_sta4vTOF);
   fChain->SetBranchAddress("liqlf", &liqlf, &b_liqlf);
   fChain->SetBranchAddress("liqld", &liqld, &b_liqld);
   fChain->SetBranchAddress("liqrf", &liqrf, &b_liqrf);
   fChain->SetBranchAddress("liqrd", &liqrd, &b_liqrd);
   fChain->SetBranchAddress("liqlTOF", &liqlTOF, &b_liqlTOF);
   fChain->SetBranchAddress("liqrTOF", &liqrTOF, &b_liqrTOF);
   // fChain->SetBranchAddress("brf", &brf, &b_brf);

   //================================Store in New ROOT file
   totnumEntry = tree->GetEntries();

   count = 0;

   saveFileName = fChain->GetDirectory()->GetName();
   saveFileName = "X_"+saveFileName;
   
   printf("Converting %s ------> %s , total Entry : %d \n", fChain->GetDirectory()->GetName(), saveFileName.Data(), totnumEntry);
   
   saveFile = new TFile( saveFileName,"recreate");
   newTree =  new TTree("tree","tree");
   
   clock.Reset();
   clock.Start("timer");
   shown = 0;

   //================================= Initialize Tree variable
  // xp   =TMath::QuietNaN();
  // yp   =TMath::QuietNaN();
  // slop1=TMath::QuietNaN();
  // slop0=TMath::QuietNaN();
   ratio1=TMath::QuietNaN();
   ratio2=TMath::QuietNaN();

   //=================================== Declare Leaf of Tree
   newTree->Branch("coinReg", &coinReg, "coinReg/I");
   newTree->Branch("vetogate", &vetogate, "vetogate/I");
  
   newTree->Branch("grx", &grx, "grx/D");
   newTree->Branch("gry", &gry, "gry/D");
   newTree->Branch("grth", &grth, "grth/D");
   newTree->Branch("grph", &grph, "grph/D");
  
   newTree->Branch("grdE1", &grdE1, "grdE1/D");
   newTree->Branch("grTOF1", &grTOF1, "grTOF1/D");

   newTree->Branch("grXC", &grXC, "grXC/D");
   newTree->Branch("grthC", &grthC, "grthC/D");
   newTree->Branch("grXAux", &grXAux, "grXAux/D");

   
   newTree->Branch("liqlf", &liqlf, "liqlf/D");
   newTree->Branch("liqld", &liqld, "liqld/D");
   newTree->Branch("liqrf", &liqrf, "liqrf/D");
   newTree->Branch("liqrd", &liqrd, "liqrd/D");
   newTree->Branch("liqlTOF", &liqlTOF, "liqlTOF/D");
   newTree->Branch("liqrTOF", &liqrTOF, "liqrTOF/D");
   ;
   newTree->Branch("ratio1", &ratio1, "ratio1/D");
   newTree->Branch("ratio2", &ratio2, "ratio2/D");

   //============================= Discrimination function

   Double_t paraL1[7] = {-23.7487397849174,12.4173924736471,-0.828505201407496,3.21874524544915e-2,-6.68727859698445e-4,7.02202234122151e-6,-2.927396710109e-8};   
   polL1 = new TF1 ("polL1", "pol6(0)", 0, 40);
   polL1->SetParameters(paraL1);

   Double_t paraL2[7] = {0.172266977459559,5.9637634482046,-0.304035887742127,1.07119395996909e-2,-2.06969268428427e-4,2.0368291152803e-6,-7.95796852547104e-9};
   polL2 = new TF1 ("polL2", "pol6(0)", 0, 40);
   polL2->SetParameters(paraL2);

   Double_t paraLN[7] = {17.1403933439602,2.82668091227555,-0.783608859798423,5.7686415838032e-2,-1.90037901151952e-3,2.95287817182544e-5,-1.77117413818978e-7};
   polLN = new TF1 ("polLN", "pol6(0)", 0, 40);
   polLN->SetParameters(paraLN);

   
   Double_t paraR1[7] = {-60.9368197474984,17.2229169053041,-0.963291218034542,2.83802883192005e-2,-3.65675700383337e-4,7.63218683391687e-7,1.44952843894597e-8};   
   polR1 = new TF1 ("polR1", "pol6(0)", 0, 40);
   polR1->SetParameters(paraR1);

   Double_t paraR2[7] = {-32.5487396988204,11.9170780214638,-0.773943598797264,0.0332804686513151,-8.4061488205942e-4,1.1347308205463e-5,-6.26354798929918e-8};
   polR2 = new TF1 ("polR2", "pol6(0)", 0, 40);
   polR2->SetParameters(paraR2);

   Double_t paraRN[7] = {2.2288133607634,2.2787473254491,-0.16022350090044,-8.96922502345884e-3,1.02117809871895e-3,-2.83696246833626e-5,2.51396547312939e-7};
   polRN = new TF1 ("polRN", "pol6(0)", 0, 40);
   polRN->SetParameters(paraRN);
   
}

Bool_t Selector_disc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Selector_disc_cxx

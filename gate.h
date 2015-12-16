void gate(TTree* tree = (TTree*)_file0->Get("tree")){
  //gROOT->Reset();
        gROOT->ProcessLine(".!date");
        gROOT->ProcessLine(".L ~/ana/yulei/rootmacro/constant.h");
        
        
        //======================================================== InPut setting
        //char * rootfile = "run1035.root";
        //Int_t Div[2] = {2,2};  //x,y
        //Int_t size[2] = {400,400}; //x,y
        
        //Double_t BGscale = 1.0;

	//====================================== Load root file
        //TFile *f0 = new TFile (rootfile, "read"); 
        //TTree *tree = (TTree*)f0->Get("tree");
        //printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
        gStyle->SetOptStat(1112211);
        gStyle->SetOptFit(1);

        tree->Print();
        
        //======================================================== Browser or Canvas
        //     TCanvas * cAna = new TCanvas("cAna", "cAna", 200, 0 , size[0]*Div[0], size[1]*Div[1]);
        //cAna->Divide(Div[0],Div[1]);
        //cAna->cd(1);    

      
        //======================================================== Cut/Gate
        TString gateStr; 
        
        gateStr.Form("0<blo1Tavg && blo1Tavg<250*%f",LAS_CH2NS[0]); TCut gateBlo1 = gateStr;
        gateStr.Form("0<blo2Tavg && blo2Tavg<250*%f",LAS_CH2NS[1]); TCut gateBlo2 = gateStr;
        gateStr.Form("0<blo3Tavg && blo3Tavg<250*%f",LAS_CH2NS[2]); TCut gateBlo3 = gateStr;
        gateStr.Form("0<blo4Tavg && blo4Tavg<250*%f",LAS_CH2NS[3]); TCut gateBlo4 = gateStr;
        
        gateStr.Form("0<sta1hTavg && sta1hTavg<400*%f",LAS_CH2NS[4]); TCut gateSta1h = gateStr;
        gateStr.Form("0<sta2hTavg && sta2hTavg<400*%f",LAS_CH2NS[5]); TCut gateSta2h = gateStr;
        gateStr.Form("0<sta1vTavg && sta1vTavg<400*%f",LAS_CH2NS[6]); TCut gateSta1v = gateStr;
        gateStr.Form("0<sta2vTavg && sta2vTavg<400*%f",LAS_CH2NS[7]); TCut gateSta2v = gateStr;
        gateStr.Form("0<sta3vTavg && sta3vTavg<400*%f",LAS_CH2NS[8]); TCut gateSta3v = gateStr;
        gateStr.Form("0<sta4vTavg && sta4vTavg<400*%f",LAS_CH2NS[9]); TCut gateSta4v = gateStr;
        
        TCut gateN = "";

        TCut gateXC_cent = "-171.79 < grXC && grXC < 107.67";
        TCut gateXC_side1 = "-331.53 < grXC && grXC < -171.79";
        TCut gateXC_side2 = "107.67 < grXC && grXC < 267.41";        
   
        //------- Graphical Cut
        TCutG * gate3He_a = new TCutG("cut3He_a", 5);
        gate3He_a->SetVarX("grTOF1");
        gate3He_a->SetVarY("grdE1");
        gate3He_a->SetPoint(0, 183.7, 305.8);
        gate3He_a->SetPoint(1, 188.4, 211.6);
        gate3He_a->SetPoint(2, 214.2, 195.4);
        gate3He_a->SetPoint(3, 225.9, 304.6);
        gate3He_a->SetPoint(4, 208.6, 343.0);
        gate3He_a->SetPoint(5, 183.7, 305.8);
	
        TCutG * gate3He_b = new TCutG("cut3He_b", 5);
        gate3He_b->SetVarX("grTOF1");
        gate3He_b->SetVarY("grdE1");
        gate3He_b->SetPoint(0, 183.7+99, 305.8);
        gate3He_b->SetPoint(1, 188.4+99, 211.6);
        gate3He_b->SetPoint(2, 214.2+99, 195.4);
        gate3He_b->SetPoint(3, 225.9+99, 304.6);
        gate3He_b->SetPoint(4, 208.6+99, 343.0);
        gate3He_b->SetPoint(5, 183.7+99, 305.8);

         //------- complex gate
        TCut gateBloAll  = "coinReg == 16" && (gateBlo1  || gateBlo2  || gateBlo3  || gateBlo4);
        TCut gateStaAll  = "vetogate == 1 && coinReg == 16 && (cut3He_a || cut3He_b) && TMath::Finite(grXC)" && (gateSta1h || gateSta2h || gateSta1v || gateSta1v || gateSta2v || gateSta3v || gateSta4v);
        TCut gateStaAll  = "vetogate == 1 && coinReg == 16 && (cut3He_a || cut3He_b) " +  (gateSta1h || gateSta2h || gateSta1v || gateSta1v || gateSta2v || gateSta3v || gateSta4v);

     
}

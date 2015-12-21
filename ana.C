{
        gROOT->Reset();
        gROOT->ProcessLine(".!date");
        gROOT->ProcessLine(".L constant.h");
        gROOT->ProcessLine(".L nuclei_mass.h");
        gROOT->ProcessLine(".L Fit_2Gauss.c");
        gROOT->ProcessLine(".L Fit_2Gauss_sub.c");
        //======================================================== InPut setting
        char * rootfile = "Y_run1035.root";
        
        Double_t BGscale = 2.0;

        //-------- Canvas control
        Bool_t analysis = 0;
        Int_t Div[2] = {1,1};  //x,y
        Int_t size[2] = {600,600}; //x,y
        Int_t cPos[2] = {2000,0}; //x,y 

        //====================================== Load root file
        TFile *f0 = new TFile (rootfile, "read"); 
        TTree *tree = (TTree*)f0->Get("tree");
        printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
        gStyle->SetOptStat(1112211);
        gStyle->SetOptFit(1);
   
        //======================================================== Cut/Gate
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


        //--------------- Simple Cut
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

        //-----------------stack neutron gate
        Double_t gnmin = -12.0, gnmax = -8.0; //ns

        gateStr.Form("(%f < sta1hTOF && sta1hTOF < %f) || (%f + 98.5 < sta1hTOF && sta1hTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu1h = gateStr;
        gateStr.Form("(%f < sta2hTOF && sta2hTOF < %f) || (%f + 98.5 < sta2hTOF && sta2hTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu2h = gateStr;
        gateStr.Form("(%f < sta1vTOF && sta1vTOF < %f) || (%f + 98.5 < sta1vTOF && sta1vTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu1v = gateStr;
        gateStr.Form("(%f < sta2vTOF && sta2vTOF < %f) || (%f + 98.5 < sta2vTOF && sta2vTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu2v = gateStr;
        gateStr.Form("(%f < sta3vTOF && sta3vTOF < %f) || (%f + 98.5 < sta3vTOF && sta3vTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu3v = gateStr;
        gateStr.Form("(%f < sta4vTOF && sta4vTOF < %f) || (%f + 98.5 < sta4vTOF && sta4vTOF < %f + 98.5)", gnmin, gnmax, gnmin, gnmax); TCut gateStaNeu4v = gateStr;

        TCut gateStaNeu = (gateStaNeu1h || gateStaNeu2h || gateStaNeu1v || gateStaNeu2v || gateStaNeu3v || gateStaNeu4v); 
       
        //-------------- other

        TCut gateGRLAS = "coinReg == 16";
        TCut gateFinite = "TMath::Finite(grXC)";

        TCut gate1 = "adc[0]<20 && adc[1]<20";
        TCut gate2 = "adc[2]<20 && adc[3]<20";
        TCut gate3 = "adc[4]<20 && adc[5]<20";
        TCut gate4 = "adc[6]<20 && adc[7]<20";
        TCut gate5 = "adc[8]<20 && adc[9]<20";
        TCut gate6 = "adc[10]<20 && adc[11]<20";

        TCut gateStaPed = gate1 + gate2 + gate3 + gate4 + gate5 + gate6;
        
        TCut gateXC_cent = "-171.79 < grXC && grXC < 107.67";
        TCut gateXC_side1 = "-331.53 < grXC && grXC < -171.79";
        TCut gateXC_side2 = "107.67 < grXC && grXC < 267.41";

        //----------- Liquid discrimination
        
        TCut gateL1 =  "liqlf < 3.75 *liqld +  25.    && liqld <= 20";
        TCut gateL2 =  "liqlf < 5./3.*liqld + 200./3. && liqld > 20";

        TCut gateR1 =  "liqrf < 4.   *liqrd + 20.     && liqrd <= 20";
        TCut gateR2 =  "liqrf < 5./3.*liqrd + 200./3. && liqrd > 20";

        TCut gateL = gateL1 || gateL2;
        TCut gateR = gateR1 || gateR2;

        TCut gateGammaL = "-85<liqlTOF && liqlTOF<-78";
        TCut gateGammaR = "-85<liqrTOF && liqrTOF<-78";

        
        //------- complex gate
        TCut gate3He = "cut3He_a || cut3He_b";
        TCut gateBloTri  = (gateBlo1  || gateBlo2  || gateBlo3  || gateBlo4);
        TCut gateStaTri  = (gateSta1h || gateSta2h || gateSta1v || gateSta1v || gateSta2v || gateSta3v || gateSta4v);
        TCut gateTem  = "vetogate == 1" && gateFinite + gateGRLAS + gate3He + gateStaTri;

        TCut gateBloAll  = gateBlo1  || gateBlo2  || gateBlo3  || gateBlo4;
        TCut gateStaAll  = gateSta1h || gateSta2h || gateSta1v || gateSta1v || gateSta2v || gateSta3v || gateSta4v;
        TCut gateSta  = gateStaAll + "vetogate == 1 && (cut3He_a || cut3He_b)" + gateGRLAS;
        
        printf("........ loaded gates\n"); 
        //====================================================== Test code

        //--------- Compiling Selector
        //gROOT->ProcessLine(".L Selector_disc.C+");
        //tree->Process("Selector_disc");

        //h1px->Draw();

        
        
        /**/
        //======================================================== Browser or Canvas
        if( analysis == 0) {
          printf("............. end of ana.C\n");
          return;
        }   
        printf(".......... start analysis \n");

        TCanvas * cAna = new TCanvas("cAna", "cAna", cPos[0], cPos[1] , size[0]*Div[0], size[1]*Div[1]);
        cAna->Divide(Div[0],Div[1]);
        cAna->cd(1);
        
        /////======================================================== analysis

        TH2F * h1 = new TH2F("h1", "h1", 200, 0, 400, 200, 0, 700);

        tree->Draw("liqlf:liqld>>h1", gateL, "colz");

        TProfile * h1px = new TProfile("h1px", "h1px", 600, 0, 400);

        h1->ProfileX("h1px")->Draw("same");

       
        /* //------------------------------------------- Fitting of Liqid discrimination
        tree->Draw("ratio1>>h1(200,0.9,1.2)", gate3He + gateL+ "liqld>20" , "colz");

        TH1F* k1;
        k1 = Fit_2Gauss_sub(h1, 250, 1.01, 0.02, 300, 1.05, 0.02); 
        cAna->cd(2);
        
        k1->Draw();
        
        cAna->cd(3);

        tree->Draw("ratio2>>h2(200,0.9,1.2)", gate3He + gateR+ "liqrd>20" , "colz");
        TH1F* k2;
        k2 = Fit_2Gauss_sub(h2, 40, 1.02, 0.02, 80, 1.07, 0.02); 
        
        cAna->cd(4);
        k2->Draw();
        //=================================================================/**/
        
		//tree->Draw("adc[0]:adc[1]>>h1h(100, -50, 150, 100, -50, 100)",	 gateGRLAS, "colz");
		//tree->Draw("adc[2]:adc[2]>>h2h(100, -50, 150, 100, -50, 100)",   gateGRLAS, "colz");
		//tree->Draw("adc[4]:adc[4]>>h1v(100, -50, 150, 100, -50, 100)",   gateGRLAS, "colz");
		//tree->Draw("adc[6]:adc[6]>>h2v(100, -50, 150, 100, -50, 100)",   gateGRLAS, "colz");
		//tree->Draw("adc[8]:adc[9]>>h3v(100, -50, 150, 100, -50, 100)",   gateGRLAS, "colz");
		//tree->Draw("adc[10]:adc[11]>>h4v(100, -50, 150, 100, -50, 100)", gateGRLAS, "colz");
 
//		tree->Draw("sqrt((adc[0]+50)*(adc[1]+50)):sta1h*120>>h1(100, -50, 150, 100, 0, 150)", gateSta, "colz");
//
//
//		tree->Draw("sqrt((adc[0]+50)*(adc[1]+50))+sqrt((adc[2]+50)*(adc[3]+50))+sqrt((adc[4]+50)*(adc[5]+50))+sqrt((adc[6]+50)*(adc[7]+50))+sqrt((adc[8]+50)*(adc[9]+50))+sqrt((adc[10]+50)*(adc[11]+50)):sta_sum>>h2(100, -0, 20, 100, 300, 600)", gateSta, "colz");
//
//		
//		tree->Draw("sqrt((adc[0]+50)*(adc[1]+50))+sqrt((adc[2]+50)*(adc[3]+50))+sqrt((adc[4]+50)*(adc[5]+50))+sqrt((adc[6]+50)*(adc[7]+50))+sqrt((adc[8]+50)*(adc[9]+50))+sqrt((adc[10]+50)*(adc[11]+50)):sta_ratio>>h3(100, -1, 1, 100, 300, 600)", gateSta, "colz");
//		tree->Draw("sta_sum:sta_ratio>>h4(100, -1, 1, 100, 0, 20)", gateSta, "colz");
//
		//tree->Draw("sta_sum:sta_ratio>>h4a(100, -1, 1, 100, 0, 20)", gateSta + !gateStaPed, "colz");


        //        tree->Draw("adc[0]:adc[1]:adc[2]:adc[3]:adc[4]:adc[5]:adc[6]:adc[7]:adc[8]:adc[9]:adc[10]:adc[11]:sta_sum:sta_ratio", gateSta, "para");

        /**/
}

{

        gROOT->Reset();
        gROOT->ProcessLine(".!date");
        gROOT->ProcessLine(".L ~/ana/yulei/rootmacro/constant.h");
        
        //======================================================== InPut setting
        char * rootfile = "run1035.root";
        Int_t Div[2] = {1,1};  //x,y
        Int_t size[2] = {400,400}; //x,y
        
        Double_t BGscale = 1.0;

	//====================================== Load root file
        TFile *f0 = new TFile (rootfile, "read"); 
        TTree *tree = (TTree*)f0->Get("tree");
        printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
        gStyle->SetOptStat(1112211);
        gStyle->SetOptFit(1);
        
        //======================================================== Browser or Canvas
        TCanvas * cAna = new TCanvas("cAna", "cAna", 2000, 0 , size[0]*Div[0], size[1]*Div[1]);
        cAna->Divide(Div[0],Div[1]);
        cAna->cd(1);
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

        TCut gateGRLAS = "coinReg == 16";

        TCut gate1 = "adc[0]<20 && adc[1]<20";
        TCut gate2 = "adc[2]<20 && adc[3]<20";
        TCut gate3 = "adc[4]<20 && adc[5]<20";
        TCut gate4 = "adc[6]<20 && adc[7]<20";
        TCut gate5 = "adc[8]<20 && adc[9]<20";
        TCut gate6 = "adc[10]<20 && adc[11]<20";

        TCut gateStaPed = gate1 + gate2 + gate3 + gate4 + gate5 + gate6;
        
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
        TCut gateBloAll  = gateBlo1  || gateBlo2  || gateBlo3  || gateBlo4;
        TCut gateStaAll  = gateSta1h || gateSta2h || gateSta1v || gateSta1v || gateSta2v || gateSta3v || gateSta4v;
        TCut gateSta  = gateStaAll + "vetogate == 1 && (cut3He_a || cut3He_b)" + gateGRLAS;
        
        
        /////======================================================== analysis
		tree->Draw("adc[0]:adc[1]>>h1h(100, -50, 150, 100, -50, 100)",	 gateGRLAS, "colz");
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
}

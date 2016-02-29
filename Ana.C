
{
        gROOT->Reset();
        //gROOT->ProcessLine(".!date");
        gROOT->ProcessLine(".L constant.h");
        gROOT->ProcessLine(".L nuclei_mass.h");
        gROOT->ProcessLine(".L Fit_2Gauss.c");
        //======================================================== InPut setting
        char * rootfile = "run1223.root";
        Int_t Div[2] = {1,1};  //x,y
        Int_t size[2] = {400,400}; //x,y
        
        Bool_t analysis = 0; // 0 = no analysis, only load root file and gate; 1 = analysis.
        
        Double_t BGscale = 1.0;
        
	    //====================================== Load root file
        TFile *f0 = new TFile (rootfile, "read"); 
        TTree *tree = (TTree*)f0->Get("tree");
        printf("=====> /// %15s //// is loaded. Total #Entry: %10d \n", rootfile,  tree->GetEntries());
        gStyle->SetOptStat(1112211);
        gStyle->SetOptFit(1);

        //======================================================== Cut/Gate      
        //------- Graphical Cut
        TCutG * gate3He1_a = new TCutG("cut3He1_a", 5);
        gate3He1_a->SetVarX("grTOF1");
        gate3He1_a->SetVarY("grdE1");
        gate3He1_a->SetPoint(0,193.772,274.187); //#1035
		gate3He1_a->SetPoint(1,193.772,245.466);
		gate3He1_a->SetPoint(2,209.298,226.228);
		gate3He1_a->SetPoint(3,209.357,253.053);
		gate3He1_a->SetPoint(4,193.772,274.187);

		//gate3He1_a->SetPoint(0,186.604,257.051); //#1034
		//gate3He1_a->SetPoint(1,186.604,232.051);
		//gate3He1_a->SetPoint(2,200.154,211.538);
		//gate3He1_a->SetPoint(3,200.154,234.615);
		//gate3He1_a->SetPoint(4,186.604,257.051);
	
        TCutG * gate3He1_b = new TCutG("cut3He1_b", 5);
        gate3He1_b->SetVarX("grTOF1");
        gate3He1_b->SetVarY("grdE1");
        gate3He1_b->SetPoint(0,193.772+98.85,274.187); //#1035
		gate3He1_b->SetPoint(1,193.772+98.85,245.466);
		gate3He1_b->SetPoint(2,209.298+98.85,226.228);
		gate3He1_b->SetPoint(3,209.357+98.85,253.053);
		gate3He1_b->SetPoint(4,193.772+98.85,274.187);

		//gate3He1_b->SetPoint(0,186.604+98.85,257.051); //#1034
		//gate3He1_b->SetPoint(1,186.604+98.85,232.051);
		//gate3He1_b->SetPoint(2,200.154+98.85,211.538);
		//gate3He1_b->SetPoint(3,200.154+98.85,234.615);
		//gate3He1_b->SetPoint(4,186.604+98.85,257.051);
        //   #1035
        TCutG *gate3He2_a = new TCutG("cut3He2_a",8);
        gate3He2_a->SetVarX("grTOF2");
        gate3He2_a->SetVarY("grdE2");
        gate3He2_a->SetPoint(0,195.517,448.418);
        gate3He2_a->SetPoint(1,195.517,412.816);
        gate3He2_a->SetPoint(2,202.037,428.639);
        gate3He2_a->SetPoint(3,207.317,433.386);
        gate3He2_a->SetPoint(4,210.819,420.728);
        gate3He2_a->SetPoint(5,210.873,465.032);
        gate3He2_a->SetPoint(6,207.425,471.361);
        gate3He2_a->SetPoint(7,195.517,448.418);

        TCutG *gate3He2_b = new TCutG("cut3He2_b",8);
        gate3He2_b->SetVarX("grTOF2");
		gate3He2_b->SetVarY("grdE2");
		gate3He2_b->SetPoint(0,195.517+98.85,448.418);
		gate3He2_b->SetPoint(1,195.517+98.85,412.816);
		gate3He2_b->SetPoint(2,202.037+98.85,428.639);
		gate3He2_b->SetPoint(3,207.317+98.85,433.386);
		gate3He2_b->SetPoint(4,210.819+98.85,420.728);
		gate3He2_b->SetPoint(5,210.873+98.85,465.032);
		gate3He2_b->SetPoint(6,207.425+98.85,471.361);
		gate3He2_b->SetPoint(7,195.517+98.85,448.418);
        /*    // #1034
        TCutG *gate3He2_a = new TCutG("cut3He2_a",5);
        gate3He2_a->SetVarX("grTOF2");
        gate3He2_a->SetVarY("grdE2");
        gate3He2_a->SetPoint(0,188.366,467.379);
        gate3He2_a->SetPoint(1,188.366,407.893);
        gate3He2_a->SetPoint(2,201.588,360.07);
        gate3He2_a->SetPoint(3,201.588,402.061);
        gate3He2_a->SetPoint(4,188.366,467.379);

        TCutG *gate3He2_b = new TCutG("cut3He2_b",5);
        gate3He2_b->SetVarX("grTOF2");
		gate3He2_b->SetVarY("grdE2");
		gate3He2_b->SetPoint(0,188.366+98.85,467.379);
        gate3He2_b->SetPoint(1,188.366+98.85,407.893);
        gate3He2_b->SetPoint(2,201.588+98.85,360.07);
        gate3He2_b->SetPoint(3,201.588+98.85,402.061);
        gate3He2_b->SetPoint(4,188.366+98.85,467.379);*/

		TCutG *gateN3He = new TCutG("cutN3He",5);//narrow 3He in #1034
        gateN3He->SetVarX("grx");
        gateN3He->SetVarY("grthC*TMath::RadToDeg()");
        gateN3He->SetPoint(0,53.0702,0.5);
        gateN3He->SetPoint(1,37.2102,-1.16049);
        gateN3He->SetPoint(2,106.133,-0.209877);
        gateN3He->SetPoint(3,120.23,1.19136);
        gateN3He->SetPoint(4,53.0702,0.5);

        TCutG * gateStaGam = new TCutG("cutsta_g",9);
        gateStaGam->SetVarX("sta_ratio");
        gateStaGam->SetVarY("sta_sum");
        gateStaGam->SetPoint(0,-0.667692,1.29108);
        gateStaGam->SetPoint(1,-0.52544,5.80139);
        gateStaGam->SetPoint(2,-0.265832,7.68068);
        gateStaGam->SetPoint(3,0.00088907,8.52636);
        gateStaGam->SetPoint(4,0.271166,7.68068);
        gateStaGam->SetPoint(5,0.427643,5.80139);
        gateStaGam->SetPoint(6,0.530775,2.13676);
        gateStaGam->SetPoint(7,0.548556,0.633322);
        gateStaGam->SetPoint(8,-0.667692,1.29108);

        //--------- Simple cut
        TString gateStr; 
        
        gateStr.Form("0<blo1Tavg && blo1Tavg<250*%f",LAS_CH2NS); TCut gateBlo1 = gateStr;
        gateStr.Form("0<blo2Tavg && blo2Tavg<250*%f",LAS_CH2NS); TCut gateBlo2 = gateStr;
        gateStr.Form("0<blo3Tavg && blo3Tavg<250*%f",LAS_CH2NS); TCut gateBlo3 = gateStr;
        gateStr.Form("0<blo4Tavg && blo4Tavg<250*%f",LAS_CH2NS); TCut gateBlo4 = gateStr;
        
        gateStr.Form("-5<sta1hTavg && sta1hTavg<400*%f",LAS_CH2NS); TCut gateSta1h = gateStr;
        gateStr.Form("-5<sta2hTavg && sta2hTavg<400*%f",LAS_CH2NS); TCut gateSta2h = gateStr;
        gateStr.Form("-5<sta1vTavg && sta1vTavg<400*%f",LAS_CH2NS); TCut gateSta1v = gateStr;
        gateStr.Form("-5<sta2vTavg && sta2vTavg<400*%f",LAS_CH2NS); TCut gateSta2v = gateStr;
        gateStr.Form("-5<sta3vTavg && sta3vTavg<400*%f",LAS_CH2NS); TCut gateSta3v = gateStr;
        gateStr.Form("-5<sta4vTavg && sta4vTavg<400*%f",LAS_CH2NS); TCut gateSta4v = gateStr;

        //-----------------stack neutron gate
        Double_t gnmin = 9.0, gnmax = 12.0; //ns #1035
        //Double_t gnmin = 13.0, gnmax = 16.0; //ns #1034

        gateStr.Form("(%f < sta1hTOFC && sta1hTOFC < %f) || (%f + 98.85 < sta1hTOFC && sta1hTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu1h = gateStr;
        gateStr.Form("(%f < sta2hTOFC && sta2hTOFC < %f) || (%f + 98.85 < sta2hTOFC && sta2hTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu2h = gateStr;
        gateStr.Form("(%f < sta1vTOFC && sta1vTOFC < %f) || (%f + 98.85 < sta1vTOFC && sta1vTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu1v = gateStr;
        gateStr.Form("(%f < sta2vTOFC && sta2vTOFC < %f) || (%f + 98.85 < sta2vTOFC && sta2vTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu2v = gateStr;
        gateStr.Form("(%f < sta3vTOFC && sta3vTOFC < %f) || (%f + 98.85 < sta3vTOFC && sta3vTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu3v = gateStr;
        gateStr.Form("(%f < sta4vTOFC && sta4vTOFC < %f) || (%f + 98.85 < sta4vTOFC && sta4vTOFC < %f + 98.85)", gnmin, gnmax, gnmin, gnmax);
        TCut gateStaNeu4v = gateStr;
        //TCut gateStaNeu2h = "(gnmin < sta2hTOF && sta2hTOF < gnmax) || (gnmin + 98.5 < sta2hTOF && sta2hTOF < gnmax + 98.5)";
        //TCut gateStaNeu1v = "(gnmin < sta1vTOF && sta1vTOF < gnmax) || (gnmin + 98.5 < sta1vTOF && sta1vTOF < gnmax + 98.5)";
        //TCut gateStaNeu2v = "(gnmin < sta2vTOF && sta2vTOF < gnmax) || (gnmin + 98.5 < sta2vTOF && sta2vTOF < gnmax + 98.5)";
        //TCut gateStaNeu3v = "(gnmin < sta3vTOF && sta3vTOF < gnmax) || (gnmin + 98.5 < sta3vTOF && sta3vTOF < gnmax + 98.5)";
        //TCut gateStaNeu4v = "(gnmin < sta4vTOF && sta4vTOF < gnmax) || (gnmin + 98.5 < sta4vTOF && sta4vTOF < gnmax + 98.5)";

        TCut gateStaNeu = (gateStaNeu1h || gateStaNeu2h || gateStaNeu1v || gateStaNeu2v || gateStaNeu3v || gateStaNeu4v); 
       

        //TCut gateXC_cent = "-171.79 < grXC && grXC < 107.67";
        TCut gateXC_cent = "-160.59 < grXC && grXC < 111.03";
        //TCut gateXC_side1 = "-331.53 < grXC && grXC < -171.79";
        TCut gateXC_side1 = "-319.83 < grXC && grXC < = -160.59";
        //TCut gateXC_side2 = "107.67 < grXC && grXC < 267.41";
        TCut gateXC_side2 = "111.03 < grXC && grXC < 270.27";
       

        TCut gateGRLAS = "coinReg == 16";
        TCut gateGRFinite = "TMath::Finite(grXC) && TMath::Finite(grthC) && TMath::Finite(gry) && TMath::Finite(grph)";
        TCut gateLASFinite = "TMath::Finite(sta_sum) && TMath::Finite(sta_ratio)";

        /* //----------- Liquid discrimination
        
        TCut gateL1 =  "liqlf < 3.75 *liqld +  25.    && liqld <= 20";
        TCut gateL2 =  "liqlf < 5./3.*liqld + 200./3. && liqld > 20";

        TCut gateR1 =  "liqrf < 4.   *liqrd + 20.     && liqrd <= 20";
        TCut gateR2 =  "liqrf < 5./3.*liqrd + 200./3. && liqrd > 20";

        TCut gateL = gateL1 || gateL2;
        TCut gateR = gateR1 || gateR2;
        */
         //------- complex gate
        TCut gate3He = "((cut3He1_a || cut3He1_b) && (cut3He2_a || cut3He2_b))";
        TCut gateVeto =  "vetogate == 1";
        TCut gateBloTri  = (gateBlo1  || gateBlo2  || gateBlo3  || gateBlo4);
        TCut gateStaTri  = (gateSta1h || gateSta2h || gateSta1v || gateSta2v || gateSta3v || gateSta4v);
 
        TCut gateAcc = "lastgr > 1000";
        TCut gateTrue = "760 < lastgr && lastgr < 820";
        TCut gateEve = "eventID < 3.8e6";
        TCut gatePhi = "-0.4 < grph*TMath::RadToDeg() && grph*TMath::RadToDeg() < 0.4";

        TCut gateTem1  =  gateEve + gatePhi + gateVeto + gateGRLAS + gateGRFinite + gateLASFinite + gateStaTri + gateStaNeu + gate3He + "!cutsta_g";
        TCut gateTem2  =  gateGRFinite+gateEve+gatePhi+gate3He;

        //------------stack position gate

        TCut gateSta1h1 = "sta1hTdif < -2";
        TCut gateSta1h2 = "-2 <= sta1hTdif && sta1hTdif < 0";
        TCut gateSta1h3 = " 0 <= sta1hTdif && sta1hTdif < 2";
        TCut gateSta1h4 = "sta1hTdif >= 2";
        TCut gateSta2h1 = "sta2hTdif < -2";
        TCut gateSta2h2 = "-2 <= sta2hTdif && sta2hTdif < 0";
        TCut gateSta2h3 = "0 <= sta2hTdif && sta2hTdif < 2";
        TCut gateSta2h4 = "sta2hTdif >= 2";

        TCut gateSta1v1 = "sta1vTdif >=0";
        TCut gateSta1v2 = "sta1vTdif < 0";
        TCut gateSta2v1 = "sta2vTdif >=0";
        TCut gateSta2v2 = "sta2vTdif < 0";
        TCut gateSta3v1 = "sta3vTdif >=0";
        TCut gateSta3v2 = "sta3vTdif < 0";
        TCut gateSta4v1 = "sta4vTdif >=0";
        TCut gateSta4v2 = "sta4vTdif < 0";



        TCut gateSta11 = gateSta1h1 || gateSta1v1;
        TCut gateSta12 = gateSta1h2 || gateSta2v1;
        TCut gateSta13 = gateSta1h3 || gateSta3v1;
        TCut gateSta14 = gateSta1h4 || gateSta4v1;
        TCut gateSta21 = gateSta2h1 || gateSta1v2;
        TCut gateSta22 = gateSta2h2 || gateSta2v2;
        TCut gateSta23 = gateSta2h3 || gateSta3v2;
        TCut gateSta24 = gateSta2h4 || gateSta4v2;
                                               

        printf("........ loaded gates\n"); 
   
        if( analysis == 0) {
          printf("............. end of Ana.C\n");
          return;
        }   

        printf(".......... start analysis \n");
        //======================================================== Browser or Canvas
        TCanvas * cAna = new TCanvas("cAna", "cAna", 200, 0 , size[0]*Div[0], size[1]*Div[1]);
        cAna->Divide(Div[0],Div[1]);
        cAna->cd(1);

        /////======================================================== analysis
        /*
        tree->Draw("ratio1>>h1(200,0.9,1.2)", gate3He + gateL+ "liqld>20" , "colz");

        Fit_2Gauss(h1, 250, 1.01, 0.02, 300, 1.05, 0.02);

        TH1F* k1 = new TH1F("k1", "k1", 200, 0.9, 1.2);  
        
        for( Int_t i = 1; i <=200; i++){
           double x = h1->GetBinCenter(i); 
           k1->Fill(x, Gauss1->Eval(x));
        }
        
        TH1F* k2 = new TH1F(*h1 - *k1);

        k2->Draw();
        

        
        tree->Draw("liqlf:liqld>>h1(500,0,140,500,0,200", gate3He, "colz");
        //tree->Draw("liqrf:liqrd>>h2(500,0,140,500,0,200", gate3He, "colz");

        cAna->cd(2);
        tree->Draw("liqlf:liqld>>h1g(500,0,140,500,0,200", gate3He + gateL, "colz");
        //tree->Draw("liqrf:liqrd>>h2g(500,0,140,500,0,200", gate3He, "colz");

        
        //TF1 * fL = new TF1("fL", 3.75*x + 25, 0, 20); 

        */
        
        //tree->Draw("grdE1:grTOF1>>h1(500, 100, 350, 500, 0, 500)", "", "colz");
        //tree->Draw("grdE1:grTOF1>>h1g(500, 100, 350, 500, 0, 500)", gate3He, "colz");
        /* tree->Draw("grthC*TMath::RadToDeg():grXC>>h2(1200,-600,600,600,-3,3)", gate3He, "colz");
        
        cAna->cd(2);
        tree->Draw("grph*TMath::RadToDeg():gry>>h3(600,-100,100,600,-3,-3)", "cut3He_a || cut3He_b", "colz");

        cAna->cd(3);
        tree->Draw("grph*TMath::RadToDeg():grth*TMath::RadToDeg()>>h4(600,-5,5,600,-5,5)", "cut3He_a || cut3He_b", "colz");

        cAna->cd(4);
        tree->Draw("gry:grx>>h5(600,-1000,1000,600,-100,100)", "cut3He_a || cut3He_b", "colz");
        */
        //tree->Draw("grthC*TMath::RadToDeg():grXC>>h3(1200,-600,600,600,-1.5,1.5)", gate3He, "colz");
        tree->Draw("grXC>>h1(1200,-600,600)",  gateGRFinite + gateEve + gatePhi + gate3He);

        Fit_2Gauss(h1, 1000, -20, 50, 10, 300, 400);
        //tree->Draw("sta_sum:sta_ratio>>s1(200,-1.1,1.1, 200,0.,35.)", "cutsta_g" + gateTem, "colz");
        //cAna->cd(2);
        //tree->Draw("sta_sum:sta_ratio>>s2(200,-1.1,1.1, 200,0.,35.)", "!cutsta_g" + gateTem + gateXC_cent + gateStaNeu, "colz");
        //tree->Draw("grXC>>h2(1200,-600,600)",gate3He + "-0.5 < grthC*TMath::RadToDeg() && grthC*TMath::RadToDeg() < 0");
        //tree->Draw("grXC>>h1(1200,-600,600)",gateTem);
        //tree->Draw("sta2hTOF>>s2(200,-50,150)", "cutsta_g" + gateTem);
        //cAna->cd(3);
        //tree->Draw("grthC*TMath::RadToDeg():grXC>>h4(1200,-600,600,600,-1.5,1.5)", gateTem + gateStaNeu,"colz");
        //tree->Draw("sta_sum:sta_ratio>>s3(200,-1.1,1.1, 200,0.,35.)","!cutsta_g"+ gateTem + (gateXC_side1 || gateXC_side2) + gateStaNeu, "colz");
        //tree->Draw("sta_sum:sta_ratio>>s3(200,-1.1,1.1, 200,0.,35.)", "!cutsta_g" + gateTem + gateStaNeu, "colz");
       
        //tree->Draw("grXC>>h3(1200,-600,600)",gate3He + "0 < grthC*TMath::RadToDeg() && grthC*TMath::RadToDeg() < 0.5");
        //cAna->cd(4);
        //tree->Draw("sta_sum:sta_ratio>>s4(200,-1.1,1.1, 200,0.,30.)", gateTem +gateStaNeu, "colz");
        //tree->Draw("grXC>>h1(1200,-600,600)",gate3He + gateXC_cent);
        //tree->Draw("sta1hTOF>>s4(200,-50,150)", gateTem + gateXC_cent + gateStaNeu);
        //tree->Draw("grthC*TMath::RadToDeg():grXC>>h4(1200,-600,600,600,-1.5,1.5)",gate3He + gateXC_cent,"colz");
        //tree->Draw("grthC*TMath::RadToDeg()>>h4a(600,-1.5,1.5)","grthC*TMath::RadToDeg() > 0.5","colz");
        //tree->Draw("sta2hTOF>>s4(200,-50,150)", "!cutsta_g" + gateTem + gateStaNeu);

	
        /*
        THStack *mS = new THStack("mS", "Stack of Ex for 22O - 18O");
        mS->Add(m22);
        mS->Add(m21);
        mS->Add(m20);
        //mS->Add(m19);
        //mS->Add(m18);
        
        leg = new TLegend(0.1,0.6,0.3,0.9);
        //leg->SetHeader("");
        leg->SetTextSize(0.05);
        leg->AddEntry(fit22o, "(^{23}F,^{22}O)", "l");
        leg->AddEntry(fit21o, "(^{23}F,^{21}O)", "l");
        leg->AddEntry(fit20o, "(^{23}F,^{20}O)", "l"); // latex OK!
        leg->Draw();
        

    */
//=================================================== Count 3He by 3 Gauss
/*
	    tree->Draw("grXAux>>h2px(1200,-600,600)", "cut3He_a || cut3He_b", "colz");

        //h2->ProjectionX("h2px")->Draw();
        
        TF1* fit = new TF1("fit", "gaus(0)+gaus(3)+gaus(6)", -600, 600);
        Double_t para[9]={1100, -20, 50, 150, 110, 120, 340, -2, 340};
        fit->SetParameters(para);
        //fit->FixParameter(1, 0);
        fit->SetLineColor(1);
        fit->SetLineStyle(1);
        fit->SetLineWidth(2);
        h2px->Fit("fit", "R");
        printf("reduced chi-squared = %f \n", fit->GetChisquare()/fit->GetNDF());

        fit->GetParameters(para);
        TF1* g1 = new TF1("g1", "gaus(0)", -600, 600); g1->SetParameters(&para[0]); g1->SetLineColor(4); g1->Draw("same");
        TF1* g2 = new TF1("g2", "gaus(0)", -600, 600); g2->SetParameters(&para[3]); g2->SetLineColor(3); g2->Draw("same");
        TF1* g3 = new TF1("g3", "gaus(0)", -600, 600); g3->SetParameters(&para[6]); g3->SetLineColor(2); g3->Draw("same");

        TH1F* k1 = new TH1F("k1", "k1", 1200, -600, 600);  
        
        for( Int_t i = 1; i <=1200; i++){
           double x = h2px->GetBinCenter(i); 
           k1->Fill(x, g2->Eval(x)+g3->Eval(x));
        }
        
        cAna->cd(2);
        k1->Draw();
        
        cAna->cd(3);
        TH1F* k2 = new TH1F(*h2px - *k1);
        k2->Draw();
*/      
/* =================================================== discrimination 
    Double_t paraL1[7] = {-23.7487397849174,12.4173924736471,-0.828505201407496,3.21874524544915e-2,-6.68727859698445e-4,7.02202234122151e-6,-2.927396710109e-8};   
   TF1* polL1 = new TF1 ("polL1", "pol6(0)", 0, 40);
   polL1->SetParameters(paraL1);

   Double_t paraL2[7] = {0.172266977459559,5.9637634482046,-0.304035887742127,1.07119395996909e-2,-2.06969268428427e-4,2.0368291152803e-6,-7.95796852547104e-9};
   TF1* polL2 = new TF1 ("polL2", "pol6(0)", 0, 40);
   polL2->SetParameters(paraL2);

   Double_t X1 = 10.;
   Double_t dpolL1 = polL1->Derivative(X1); 
   Double_t dpolL2 = polL2->Derivative(X1);

   TF1 * L1 = new TF1("L1", "pol1", -100, 40);
   L1->SetParameter(0,polL1->Eval(X1)-dpolL1*X1);
   L1->SetParameter(1,dpolL1);
   
   TF1 * L2 = new TF1("L2", "pol1", -100, 40);
   L2->SetParameter(0,polL2->Eval(X1)-dpolL2*X1);
   L2->SetParameter(1,dpolL2);

   L1->SetLineColor(4);L1->Draw("");
   L2->SetLineColor(6);L2->Draw("same");

   polL1->Draw("same");
   polL2->Draw("same");
   
   
   Double_t xp = -(polL1->Eval(X1)-dpolL1*X1-polL2->Eval(X1)+dpolL2*X1)/(dpolL1-dpolL2);
   Double_t yp = dpolL1*(xp-X1) + polL1->Eval(X1);

   printf("xp, yp: %f, %f\n", xp, yp);

*/      

/*++++++++++++++++++++++++++++++++++++++++++++++ Momentum*/
	/*	
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p22(  9, -300, 300)", gate  + cut22o +"TMath::Abs(Ex)<20") ; //o22->SetTitle(gateStr + " X 22O");
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p21(  18, -300, 300)", gate  + cut21o +"TMath::Abs(Ex-8)<20") ; //o21->SetTitle(gateStr + " X 21O");
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p20(  18, -300, 300)", gate  + cut20o +"TMath::Abs(Ex-15)<20") ; //o20->SetTitle(gateStr + " X 20O");
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p22c( 9, -300, 300)", gatec + cut22o +"TMath::Abs(Ex)<20") ; //o22c->SetTitle(gateStr + " X 22O");
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p21c( 18, -300, 300)", gatec + cut21o +"TMath::Abs(Ex-8)<20") ; //o21c->SetTitle(gateStr + " X 21O");
        //tree->Draw("(0.002*kE + kp)/sqrt(1-0.002*0.002)>>p20c( 18, -300, 300)", gatec + cut20o +"TMath::Abs(Ex-15)<20") ; //o20c->SetTitle(gateStr + " X 20O");

        tree->Draw("kMomt>>p22(  40, 0, 400)", gate  + "cut22o" +"TMath::Abs(Ex)<20") ; //o22->SetTitle(gateStr + " X 22O");
        tree->Draw("kMomt>>p21(  40, 0, 400)", gate  + "cut21o" +"TMath::Abs(Ex-10)<10") ; //o21->SetTitle(gateStr + " X 21O");
        tree->Draw("kMomt>>p20(  40, 0, 400)", gate  + "cut20o" +"TMath::Abs(Ex-20)<20") ; //o20->SetTitle(gateStr + " X 20O");
        tree->Draw("kMomt>>p22c( 40, 0, 400)", gatec + "cut22o" +"TMath::Abs(Ex)<20") ; //o22c->SetTitle(gateStr + " X 22O");
        tree->Draw("kMomt>>p21c( 40, 0, 400)", gatec + "cut21o" +"TMath::Abs(Ex-10)<10") ; //o21c->SetTitle(gateStr + " X 21O");
        tree->Draw("kMomt>>p20c( 40, 0, 400)", gatec + "cut20o" +"TMath::Abs(Ex-20)<20") ; //o20c->SetTitle(gateStr + " X 20O");

        p22c->Scale(BGscale); TH1F* q22 = new TH1F(*p22 - *p22c); q22->SetLineColor(1); q22->SetName("q22"); q22->SetTitle("23F(p,2p)22Ogs");
        p21c->Scale(BGscale); TH1F* q21 = new TH1F(*p21 - *p21c); q21->SetLineColor(1); q21->SetName("q21"); q21->SetTitle("23F(p,2p)22O*->21O*+n, 6.86");
        p20c->Scale(BGscale); TH1F* q20 = new TH1F(*p20 - *p20c); q20->SetLineColor(1); q20->SetName("q20"); q20->SetTitle("23F(p,2p)22O*->20O*+2n, 10.65");

        q22->SetXTitle("k [MeV/c]"); q22->SetYTitle("count / 10 MeV"); q22->SetTitle("Ex = 0 MeV");
        q21->SetXTitle("k [MeV/c]"); q21->SetYTitle("count / 10 MeV"); q21->SetTitle("Ex = 10 MeV,");
        q20->SetXTitle("k [MeV/c]"); q20->SetYTitle("count / 10 MeV"); q20->SetTitle("Ex = 20 MeV");
        
        Ex00->Draw("k>>hk0a(40, 0, 400)", "xsec1d5", "E"); hk0a->SetLineColor(2); 
        Ex00->Draw("k>>hk0b(40, 0, 400)", "xsec1p1", "E"); hk0b->SetLineColor(3); 
        Ex00->Draw("k>>hk0c(40, 0, 400)", "xsec1p3", "E"); hk0c->SetLineColor(4); 
        Ex00->Draw("k>>hk0d(40, 0, 400)", "xsec2s1", "E"); hk0d->SetLineColor(6);
        Ex10->Draw("k>>hk1a(40, 0, 400)", "xsec1d5", "E"); hk1a->SetLineColor(2); 
        Ex10->Draw("k>>hk1b(40, 0, 400)", "xsec1p1", "E"); hk1b->SetLineColor(3); 
        Ex10->Draw("k>>hk1c(40, 0, 400)", "xsec1p3", "E"); hk1c->SetLineColor(4); 
        Ex10->Draw("k>>hk1d(40, 0, 400)", "xsec2s1", "E"); hk1d->SetLineColor(6); 
        Ex20->Draw("k>>hk2a(40, 0, 400)", "xsec1d5", "E"); hk2a->SetLineColor(2); 
        Ex20->Draw("k>>hk2b(40, 0, 400)", "xsec1p1", "E"); hk2b->SetLineColor(3); 
        Ex20->Draw("k>>hk2c(40, 0, 400)", "xsec1p3", "E"); hk2c->SetLineColor(4); 
        Ex20->Draw("k>>hk2d(40, 0, 400)", "xsec2s1", "E"); hk2d->SetLineColor(6); 
        
        Double_t temp = hk0a->GetMaximum(); 
        hk0a->Scale(q22->GetMaximum()/temp);
        hk0b->Scale(q22->GetMaximum()/temp);
        hk0c->Scale(q22->GetMaximum()/temp);
        hk0d->Scale(q22->GetMaximum()/temp);
        Double_t temp = hk1b->GetMaximum(); 
        hk1a->Scale(q21->GetMaximum()/temp);
        hk1b->Scale(q21->GetMaximum()/temp);
        hk1c->Scale(q21->GetMaximum()/temp);
        hk1d->Scale(q21->GetMaximum()/temp);
        Double_t temp = hk2c->GetMaximum(); 
        hk2a->Scale(q20->GetMaximum()/temp);  
        hk2b->Scale(q20->GetMaximum()/temp);  
        hk2c->Scale(q20->GetMaximum()/temp);
        hk2d->Scale(q20->GetMaximum()/temp);
        
        q22->SetMaximum(hk0a->GetMaximum()*1.1);
        q21->SetMaximum(hk1a->GetMaximum()*1.1);
        q20->SetMaximum(hk2a->GetMaximum()*1.1);
        
        
        q22->Draw();hk0a->Draw("same");hk0b->Draw("same");hk0c->Draw("same"); hk0d->Draw("same");
//        q21->Draw();hk1a->Draw("same");hk1b->Draw("same");hk1c->Draw("same"); hk1d->Draw("same");
*/

/*++++++++++++++++++++++++++++++++++++++++++++++ Asymmetry*//*
        TCutG* cutAsy = gate20o;
        TString plotTitle = "(23F,20O)"; 
        TCut energy = "abs(Ex-18)<10";
                
        cutAsy->SetName("cutAsy");

        tree->Draw("theta1>>h2aL(5,20,70)", gate + "cutAsy" + energy + "runNum<=36", "colz");
        tree->Draw("theta2>>h2aR(5,20,70)", gate + "cutAsy" + energy + "runNum<=36", "colz");
        
        tree->Draw("theta1>>h2aLg(5,20,70)", gatec + "cutAsy" + energy + "runNum<=36", "colz");
        tree->Draw("theta2>>h2aRg(5,20,70)", gatec + "cutAsy" + energy + "runNum<=36", "colz");
        
        TH1F* h2aLs = new TH1F(*h2aL - *h2aLg); h2aLs->SetName("h2aLs"); h2aLs->SetLineColor(4);
        TH1F* h2aRs = new TH1F(*h2aR - *h2aRg); h2aRs->SetName("h2aRs"); h2aRs->SetLineColor(2);
        
        h2aLs->SetXTitle("theta [deg]"); h2aLs->SetYTitle("count / 10 deg"); h2aLs->SetTitle(plotTitle + "Spin up");
        h2aRs->SetXTitle("theta [deg]"); h2aRs->SetYTitle("count / 10 deg"); h2aRs->SetTitle(plotTitle + "Spin up");
        
        cAna->cd(1);
        h2aRs->Draw("");
        h2aLs->Draw("same");
        
        
        tree->Draw("theta1>>h2bL(5,20,70)", gate + "cutAsy" + energy + "runNum>36", "colz");
        tree->Draw("theta2>>h2bR(5,20,70)", gate + "cutAsy" + energy + "runNum>36", "colz");
        
        tree->Draw("theta1>>h2bLg(5,20,70)", gatec + "cutAsy" + energy + "runNum>36", "colz");
        tree->Draw("theta2>>h2bRg(5,20,70)", gatec + "cutAsy" + energy + "runNum>36", "colz");
        
        TH1F* h2bLs = new TH1F(*h2bL - *h2bLg); h2bLs->SetName("h2bLs"); h2bLs->SetLineColor(4);
        TH1F* h2bRs = new TH1F(*h2bR - *h2bRg); h2bRs->SetName("h2bRs"); h2bRs->SetLineColor(2);
        
        h2bLs->SetXTitle("theta [deg]"); h2bLs->SetYTitle("count / 10 deg"); h2bLs->SetTitle(plotTitle + "Spin down");
        h2bRs->SetXTitle("theta [deg]"); h2bRs->SetYTitle("count / 10 deg"); h2bRs->SetTitle(plotTitle + "Spin down");
        
        Ex10->Draw("theta1>>g1(5,20,70)", "xsec1p1*(1+asym1p1)","");
        Ex10->Draw("theta2>>g2(5,20,70)", "xsec1p1*(1+asym1p1)","");
        
        cAna->cd(2);
        h2bLs->Draw("E"); h2bLs->Draw("same");
        h2bRs->Draw("E same"); h2bRs->Draw("same");
        
        cAna->cd(1);
        h2aLs->Draw("E"); h2aLs->Draw("same");
        h2aRs->Draw("E same"); h2aRs->Draw("same");
        
        //--------------------calculate AyP
        Int_t nBin = h2aLs->GetNbinsX();
        
        Double_t yLu[nBin], yRu[nBin], yLd[nBin], yRd[nBin];
        Double_t xBin[nBin], AyP[nBin];
        
        for ( Int_t l = 1; l <= nBin; l++){
                Int_t i = l-1;
                xBin[i] = h2aLs->GetBinCenter(l);
                yLu[i]=h2aLs->GetBinContent(l);
                yRu[i]=h2aRs->GetBinContent(l);
                yLd[i]=h2bLs->GetBinContent(l);
                yRd[i]=h2bRs->GetBinContent(l);
                
                Double_t tYL = g1->GetBinContent(l);
                Double_t tYR = g2->GetBinContent(l);
                
                Double_t tAy = (tYL-tYR)/(tYL+tYR);
                
                Double_t YL = TMath::Sqrt(yLu[i]*yRd[i]);
                Double_t YR = TMath::Sqrt(yRu[i]*yLd[i]);
                
                Double_t dYL = TMath::Sqrt(yLu[i] + yRd[i])/2;
                Double_t dYR = TMath::Sqrt(yLd[i] + yRu[i])/2;
                
                AyP[i] = (YL - YR)/(YL + YR);
                
                Double_t dAyP = TMath::Sqrt(YL*YL*dYR*dYR + YR*YR*dYL*dYL)/TMath::Power(YL + YR,2);
                
                Double_t Ay = AyP[i]/0.3;
                Double_t dAy = Ay * TMath::Sqrt( TMath::Power(dAyP/AyP[i],2) + TMath::Power(dAyP/AyP[i],2));
                
                printf("x:%4.0f, yLu:%4.0f, yRu:%4.0f, yLd:%4.0f, yRd:%4.0f, AyP:%6.4f, dAyp:%10.6f, Ay:%10.6f, dAy:%10.6f, tAy:%10.6f\n", xBin[i], yLu[i], yRu[i], yLd[i], yRd[i], AyP[i], dAyP, Ay, dAy, tAy); 
        }

        */

        printf("............. end of Ana.C\n");
}

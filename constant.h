//_______________________________Exp Parameters
  const Double_t GR_CH2NS[2] = {0.123,0.123};
  const Double_t LAS_CH2NS[10] = {0.1437,0.1451,0.1448,0.1448,0.1448,0.1434,0.1441,0.1437,0.1448,0.1444};
//const Double_t LAS_CH2NS[10] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  
  const Double_t GR_TOF1_OFFSET = 241.545; //ns
  const Double_t BLOCK_TOF_OFFSET[4] = {103.87, 104.37, 102.87, 100.87}; //ns
  const Double_t STACK_TOF_OFFSET[6] = {146.0+2.37, 144.4+2.37, 139.0+2.37, 135.6+2.37, 139.0+2.37, 136.8+2.37}; //ns
  const Double_t GR_LENGTH = 20000.; // mm 
  const Double_t GR_RHO = 3000.; //mm 
  const Double_t BRHO_3HE = GR_RHO*0.548506; // T.m
  
//_____________________________ Matrix Element

  const Double_t X_X = -0.41682; //mm/mm 
  const Double_t X_A =  0;       //mm/rad 
  const Double_t X_D = 15457.38; //mm/1% 

  const Double_t A_X = -1.34006; //rad/mm 
  const Double_t A_A = -2.3991;  //rad/rad 
  const Double_t A_D =  1.13316; //rad/1% 


//_____________________________ Universial constant  
  const double deg2rad = 3.141592654/180;
  const Double_t cVAC = 299.792458; // mm/ns
  const Double_t mp = 938.272; // MeV/c^2
  const Double_t mn = 939.565; // MeV/c^2
  const Double_t md = 1875.613;
  const Double_t m3He = 2808.391;
  const Double_t mt = 2808.921;
  const Double_t m4He = 3727.379;
  const Double_t m10B = 9324.436;
  const Double_t m10C = 9327.573;
  const Double_t m12C = 11174.863;
  const Double_t m14N = 13040.203;
  const Double_t m14O = 13044.836;
  const Double_t m16O = 14895.08;
  const Double_t amu = 931.494;   

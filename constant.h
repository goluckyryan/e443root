#ifndef constant
#define constant

//_______________________________Exp Parameters
const double GR_CH2NS[2] = {0.123,0.123};
const double LAS_CH2NS[10] = {0.1437,0.1451,0.1448,0.1448,0.1448,0.1434,0.1441,0.1437,0.1448,0.1444};
//const double LAS_CH2NS[10] = {1., 1., 1., 1., 1., 1., 1., 1., 1., 1.};
  
const double GR_TOF1_OFFSET = 241.545; //ns
const double BLOCK_TOF_OFFSET[4] = {103.87, 104.37, 102.87, 100.87}; //ns
const double STACK_TOF_OFFSET[6] = {146.0+2.37, 144.4+2.37, 139.0+2.37, 135.6+2.37, 139.0+2.37, 136.8+2.37}; //ns
//  const double BLOCK_TOF_OFFSET[4] = {0., 0., 0., 0.}; //ns
//  const double STACK_TOF_OFFSET[6] = {0., 0., 0., 0., 0., 0.}; //ns
  
const double GR_LENGTH = 20000.; // mm 
const double GR_RHO = 3000.; //mm 
const double BRHO_3HE = GR_RHO*0.548506; // T.m
  
//_____________________________ Matrix Element

const double X_X = -0.41682; //mm/mm 
const double X_A =  0;       //mm/rad 
const double X_D = 15457.38; //mm/1% 

const double A_X = -1.34006; //rad/mm 
const double A_A = -2.3991;  //rad/rad 
const double A_D =  1.13316; //rad/1% 


//_____________________________ Universial constant  
const double deg2rad = 3.141592654/180;
const double rad2deg = 1/deg2rad;
const double cVAC = 299.792458; // mm/ns

//-------- or using Nucleus_Mass(Z,A) in nucleus_mass.h
const double mp = 938.272; // MeV/c^2
const double mn = 939.565; // MeV/c^2
const double md = 1875.613;
const double m3He = 2808.391;
const double mt = 2808.921;
const double m4He = 3727.379;
const double m10B = 9324.436;
const double m10C = 9327.573;
const double m12C = 11174.863;
const double m14N = 13040.203;
const double m14O = 13044.836;
const double m16O = 14895.08;
const double amu = 931.494;



#endif

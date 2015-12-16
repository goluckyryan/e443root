
#ifndef NUCLEIMASS_H
#define NUCLEIMASS_H

#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>  //atoi
#include "constant.h"

using namespace std;

const string dataPath="/home/e443/ana/yulei/rootmacro/mass12.txt";

// about the mass12.txt
// Mass Excess = (ATOMIC MASS - A)*amu | e.g. n : (1.088664.91585E-6-1)*amu
//  mass excess uncertaintly
// BEA = (Z*M(1H) + N*M(1n) - Me(A,Z))/A , Me is the mass with electrons
// BEA = (Z*mp + N*mn - M(A,Z))/A , M is the mass without electrons

double Nucleus_Mass(int Z, int A)
{
  string line;
  int    lineNum=0;
  int    list_A, list_Z;
  double list_BEA, mass;

  ifstream myfile;
  int    flag=0;

  int numLineStart = 40;
  int numLineEnd  = 3392;  

  if ( A >= 50 && A < 100) numLineStart = 447; //545;
  if ( A >=100 && A < 150) numLineStart = 1072;//1100;
  if ( A >=150 && A < 200) numLineStart = 1833;//1899;
  if ( A >=200 ) numLineStart = 2534;//2622;

  myfile.open(dataPath.c_str());

  if (myfile.is_open())
  {
    while (/*! myfile.eof() &&*/ flag == 0 && lineNum <numLineEnd)
    {
      lineNum ++ ;
      //printf("%3d  ",lineNum);
      getline (myfile,line);

      if (lineNum >= numLineStart ){
        list_Z     = atoi((line.substr(10,5)).c_str());
      	list_A     = atoi((line.substr(15,5)).c_str());
      	list_BEA = atof((line.substr(55,11)).c_str());
      	if ( A == list_A && Z == list_Z) {
      		mass = list_Z*mp + (list_A-list_Z)*mn - list_BEA/1000*list_A;
     		flag = 1;
      	}else if ( list_A > A) {
          break;
        }

      }
    }
    myfile.close();
  }
  else printf("Unable to open %s\n", dataPath.c_str());

  if (flag == 1){
    return mass;
  }else{
  	return -404;
  }

  return 0;
}

string Nucleus_Name(int Z, int A){

  string line, name;
  int    lineNum=0;
  int    list_A, list_Z;
  ifstream myfile;
  int    flag=0;

  int numLineStart = 40;
  int numLineEnd  = 3392;

  if ( A >= 50 && A < 100) numLineStart = 545;
  if ( A >=100 && A < 150) numLineStart = 1100;
  if ( A >=150 && A < 200) numLineStart = 1899;
  if ( A >=200 ) numLineStart = 2622;

  myfile.open(dataPath.c_str());

  if (myfile.is_open())
  {
    while ( flag == 0 && lineNum <numLineEnd)
    {
      lineNum ++ ;
      getline (myfile,line);

      if (lineNum >= numLineStart ){
      	list_Z     = atoi((line.substr(10,5)).c_str());
      	list_A     = atoi((line.substr(15,5)).c_str());
      	name       = line.substr(19,3);
      	if ( A == list_A && Z == list_Z) {
     		   flag = 1;
      	}else if ( list_A > A) {
          break;
        }

      }
    }
    myfile.close();
  }
  else printf("Unable to open %s\n", dataPath.c_str());

  if (flag == 1){
      return name;
  }else{
  	   return "---";
  }

}

//Sp is seperation energy for breaking in to a nucleus with Nn neutron and Np proton. 
double Sp (int Z, int A, int Nn, int Np){

  if ( A == 1 ) return 0;

  double massP, massD1, massD2;

  massP = Nucleus_Mass(Z,A);
  massD1 = Nucleus_Mass(Z-Np,A-Np-Nn);
  massD2 = Nucleus_Mass(Np, Np+Nn);

  if ( massP != -404 && massD1 != -404 && massD2 != -404) {
    return massD1 + massD2 - massP;
  }else{
    return -0.000;
  }

}

// SpBreak is seperation energy for breaking in to Nn neutron and Np proton. 
double SpBreak (int Z, int A, int Nn, int Np){
  double mass, massD;

  mass  = Nucleus_Mass(Z,A);
  massD = Nucleus_Mass(Z-Np, A-Np-Nn);

  if( mass != -404 && massD != -404){
    return massD + Nn*mn + Np*mp - mass;
  }else{
    return -0.000;
  }
}

#endif

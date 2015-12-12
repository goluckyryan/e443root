#include "nuclei_mass.h"
using namespace std;

int main (int argc, char *argv[])
{

   int    Z, A, Za, Aa;
   double mass;
   string name = "--";

   if (argc == 3){
      A= atoi(argv[1]);
      Z= atoi(argv[2]);
   }else if ( argc == 5){
      A= atoi(argv[1]);
      Z= atoi(argv[2]);
      Aa= atoi(argv[3]);
      Za= atoi(argv[4]);      
   }else{
      cout << "./nuclei_mass.o A Z" << endl;
      cout << "./nuclei_mass.o A Z A' Z'" << endl;
      exit(0);
   }

   mass = Nucleus_Mass(Z,A);
   name = Nucleus_Name(Z,A);

   double BE = (Z*mp + (A-Z)*mn - mass)/A;

   if (mass > 0){
     printf(" mass of \e[31m%3d%s\e[m nucleus (Z,A)=(%3d,%3d) is \e[31m%12.5f\e[m MeV, BE/A=%7.5f\n",A,name.c_str(), Z,A, mass, BE);     
     printf(" mass in amu : %12.5f u\n",mass/amu);
     printf(" S1p: %8.4f| S1n: %8.4f| S(2H ): %8.4f| S1p1n : %8.4f\n", SpBreak(Z,A, 0, 1), SpBreak(Z,A, 1, 0), Sp(Z,A, 1, 1), SpBreak(Z,A,1,1));
     printf(" S2p: %8.4f| S2n: %8.4f| S(3He): %8.4f| S(3H) : %8.4f\n", SpBreak(Z,A, 0, 2), SpBreak(Z,A, 2, 0), Sp(Z,A, 1, 2), Sp(Z,A, 2, 1));
     printf(" S3p: %8.4f| S3n: %8.4f| S(4He): %8.4f|\n", SpBreak(Z,A, 0, 3), SpBreak(Z,A, 3, 0), Sp(Z,A, 2, 2));
     printf(" S4p: %8.4f| S4n: %8.4f| \n", SpBreak(Z,A, 0, 4), SpBreak(Z,A, 4, 0));
     
   }else{
      printf("Error %6.0f, no nucleus with (Z,A) = (%3d,%3d). \n", mass, Z, A);
   }

   
   if( argc == 5 ){
      string name2 = Nucleus_Name(Za, Aa);
      string name3 = Nucleus_Name(Z-Za, A-Aa);
      printf(" S(%2d%2s) = S(%2d%2s) : %8.4f MeV\n", Aa, name2.c_str(), A-Aa, name3.c_str(), Sp(Z,A,Aa-Za,Za));
      printf(" Q-factor = m(%2d%2s) - m(%2d%2s) - m(%2d%2s) = %8.4f MeV\n",A, name.c_str(), Aa, name2.c_str(), A-Aa, name3.c_str(),Nucleus_Mass(Z,A)-Nucleus_Mass(Za,Aa)-Nucleus_Mass(Z-Za,A-Aa));

            
   }
   return 0;
}


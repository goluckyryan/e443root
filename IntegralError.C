//____________________________________________________________________
double IntegralError(int npar, double * c, double * errPar, 
   double * covMatrix = 0) {   
// calculate the error on the integral given the parameter error and 
// the integrals of the gradient functions c[] 

   double err2 = 0; 
   for (int i = 0; i < npar; ++i) { 
      if (covMatrix == 0) { // assume error are uncorrelated
         err2 += c[i] * c[i] * errPar[i] * errPar[i]; 
      } else {
         double s = 0; 
         for (int j = 0; j < npar; ++j) {
            s += covMatrix[i*npar + j] * c[j]; 
         }
         err2 += c[i] * s; 
      }
   }

   return TMath::Sqrt(err2);
}

//____________________________________________________________________
double df_dPar(double * x, double * p) { 
   // derivative of the function w.r..t parameters
   // use calculated derivatives from TF1::GradientPar

   double grad[NPAR]; 
   // p is used to specify for which parameter the derivative is computed 
   int ipar = int(p[0] ); 
   //printf("p[0] : %f, ipar: %d \n",p[0], ipar);
   assert (ipar >=0 && ipar < NPAR );

   assert(total);
   total->GradientPar(x, grad);

   return grad[ipar]; 
}





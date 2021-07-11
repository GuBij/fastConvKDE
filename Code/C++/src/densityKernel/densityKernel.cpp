 #include <iostream>
 #include <vector>
 #include <algorithm>
 #include <mathimf.h>
 #include "densityKernel.h"

namespace densityKernel
{

 void bandWidth( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc, bool useProdKernel )
 {
    if ( useProdKernel )
	bandWidthProdEpanechnikov(bandWidthVec, distrSigmaSlaves, noProcs, dataSizeProc);
    else
	bandWidth3DKernel(bandWidthVec, distrSigmaSlaves, noProcs, dataSizeProc);
 }

 void bandWidth3DKernel( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
        varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] /= nop; 
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    double A = (varPos[0]*varPos[0]*varPos[1]*varPos[1]+varPos[0]*varPos[0]*varPos[2]*varPos[2]+varPos[2]*varPos[2]*varPos[1]*varPos[1])/(varPos[0]*varPos[0]*varPos[1]*varPos[1]*varPos[2]*varPos[2])+0.5*(varPos[0]*varPos[1]+varPos[1]*varPos[2]+varPos[0]*varPos[2])*(varPos[0]*varPos[1]+varPos[1]*varPos[2]+varPos[0]*varPos[2])/(varPos[0]*varPos[1]*varPos[2]*varPos[0]*varPos[1]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[1] = sqrt(varPos[1]);
    varPos[2] = sqrt(varPos[2]); 
    A = 8*varPos[0]*varPos[1]*varPos[2]*2/A;
    A *= 3*15/14*sqrt(acos(-1))*7*7;
    *bandWidthVec[0] = pow(A/nop,1.0/7.0);
    *bandWidthVec[1] = *bandWidthVec[0];
 }

 void bandWidthProdKernel( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[2] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1];
	varPos[1] += distrSigmaSlaves[i*dataSizeProc+2];
    }
    varPos[0] = sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop);

    vector<double> distrSigma(2*nop,0);
    vector<double> distrSigmaVec(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    unsigned insertAt = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
        unsigned ii = i*dataSizeProc, nopi = long(distrSigmaVec[ii]);
        move(distrSigmaVec.begin()+ii+3,distrSigmaVec.begin()+ii+2+nopi,distrSigma.begin()+insertAt);
        move(distrSigmaVec.begin()+ii+3+(dataSizeProc-3)/2,distrSigmaVec.begin()+ii+2+(dataSizeProc-3)/2+nopi,distrSigma.begin()+nop+insertAt);
        insertAt += distrSigmaVec[ii];
    }
    distrSigmaVec.clear();
    sort(distrSigma.begin(),distrSigma.begin()+nop-1);
    sort(distrSigma.begin()+nop,distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop);
    double IQR = (distrSigma[m75]-distrSigma[m25])/1.34;
    if ( nop > 100 && IQR < varPos[0] && IQR > 0)
      *bandWidthVec[0] = IQR;
    else
      *bandWidthVec[0] = varPos[0];

    IQR = (distrSigma[nop+m75]-distrSigma[nop+m25]); //1.34;

      double norm2ndDeriv = 0;

    double A = 4*20*20.0/7*4;
    *bandWidthVec[0] *= 0.85*pow(A/nop,1.0/6);

    double beta = 0, alpha = 0, sign = -1.0;
    for ( unsigned k = 0; k < 7; k++)
    {
        sign *= -1.0;
        beta += sign/((2*k+1)*fac(k)*fac(6-k));
        if (k < 4)
          alpha += sign/((2*k+3)*fac(k)*fac(3-k));
    }
    beta *= 2*fac(6); alpha *= 2*6;
    if ( norm2ndDeriv == 0 )
    {
        A = 8*sqrt(acos(-1))/3*beta/(alpha*alpha);
        *bandWidthVec[1] = 0.85*varPos[1]*pow(A/nop,0.2);
    }
    else
    {
        A = beta/(alpha*alpha*norm2ndDeriv);
        *bandWidthVec[1] = pow(A/nop,0.2);
    }
 }

 void bandWidthProdEpanechnikov( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
	varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop); 
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    double A = (varPos[0]*varPos[0]+varPos[2]*varPos[2])/(varPos[0]*varPos[0]*varPos[2]*varPos[2])+0.5*(varPos[0]+varPos[2])*(varPos[0]+varPos[2])/(varPos[0]*varPos[2]*varPos[0]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[2] = sqrt(varPos[2]); 
    A = 4*varPos[0]*varPos[2]*2/A; 
    A *= 2*4/3*6*6; 
    *bandWidthVec[0] = pow(A/nop,1.0/6);
      *bandWidthVec[1] = varPos[1]; 

   A = 8*sqrt(acos(-1))*5; // 8*sqrt(acos(-1))/3*3/5*5*5
   *bandWidthVec[1] *= pow(A/nop,0.2);
//   *bandWidthVec[1] = pow(15.0/(<bandWidth>*nop),0.2); //Replace <bandWidth> by the appropriate value according to Table 3 for the FD simulations
 }

 void bandWidthGaussian( double** bandWidthVec, double* distrSigmaSlaves, int noProcs, unsigned long dataSizeProc )
 {
    unsigned long nop = 0;
    double varPos[3] = {};
    for (unsigned i = 0; i < noProcs; i++ )
    {
        nop += distrSigmaSlaves[i*dataSizeProc];
        varPos[0] += distrSigmaSlaves[i*dataSizeProc+1]; // x-dir
        varPos[1] += distrSigmaSlaves[i*dataSizeProc+2]; // z-dir
        varPos[2] += distrSigmaSlaves[i*dataSizeProc+3]; // y-dir
    }
    varPos[0] /= nop; // sqrt(varPos[0]/nop);
    varPos[1] = sqrt(varPos[1]/nop);
    varPos[2] /= nop; // sqrt(varPos[2]/nop);

    vector<double> distrSigma(nop,0);
    vector<double> distrSigmaVec(distrSigmaSlaves, distrSigmaSlaves + noProcs*dataSizeProc);
    unsigned insertAt = 0;
    for (unsigned i = 0; i < noProcs; i++ )
    {
        unsigned ii = i*dataSizeProc, nopi = long(distrSigmaVec[ii]);
        move(distrSigmaVec.begin()+ii+4,distrSigmaVec.begin()+ii+3+nopi,distrSigma.begin()+insertAt);
        insertAt += distrSigmaVec[ii];
    }
    distrSigmaVec.clear();
    sort(distrSigma.begin(),distrSigma.end());

    unsigned long m25=floor(0.25*nop), m75=floor(0.75*nop);
    double IQR = (distrSigma[m75]-distrSigma[m25])/1.34;

    double A = (varPos[0]*varPos[0]+varPos[2]*varPos[2])/(varPos[0]*varPos[0]*varPos[2]*varPos[2])+0.5*(varPos[0]+varPos[2])*(varPos[0]+varPos[2])/(varPos[0]*varPos[2]*varPos[0]*varPos[2]);
    varPos[0] = sqrt(varPos[0]);
    varPos[2] = sqrt(varPos[2]);
    A = 2*varPos[0]*varPos[2]*2/A;
    *bandWidthVec[0] = pow(A/nop,1.0/6);

    if ( nop > 100 && IQR < varPos[1] && IQR > 0)
      *bandWidthVec[1] = IQR;
    else
      *bandWidthVec[1] = varPos[1];

   A = 4.0/3.0; 
   *bandWidthVec[1] *= pow(A/nop,0.2);
 }

 double weight(string type)
 {
   double w;

   if (type == "bw")
   {
     w = 720.0/675675; //beta value
     w *= 3*12006225/(4*acos(-1)*36); //inverse alpha squared value
     w *= 32*pow(acos(-1),1.5)/(3*(2+3));
     w = 0.85*pow(3*w,1.0/7.0);
//     w *= 0.85*pow(nop,-1/7);
   }
   else if (type == "kernel")
   {
     w = 0.25*3*315/(48*acos(-1));
   }

   return w;
 }

 unsigned fac(unsigned n)
 {
   unsigned nfac=1;

   if (n==0)
	return nfac;

   for (unsigned i = 1; i < (n+1); i++)
	nfac *= i;

   return nfac;
 }

 double integralKernel( bool useProdKernel )
 {
   double intK;

   if ( useProdKernel )
     intK = integralProdEpanechnikov();
   else 
     intK = integral3DEpanechnikov();

   return intK;
 }

 double integral3DKernel ( ) //double zc, double h, double zMax)
 {
   unsigned a = 3; double zLow, zUp;
/*
   if (zc < h)
	zLow = -zc/h;
   else
*/
	zLow = -1;
/*
   if (zc+h > zMax)
	zUp = (zMax-zc)/h;
   else
*/
	zUp = 1;

   double intK = (zUp-zLow)/fac(a+1), zUpPow = zUp, zLowPow = zLow;
   int sign = 1;
   for (unsigned k = 1; k < (a+2); k++)
   {
     sign *= -1;
//     zUpPow *= zUp*zUp; zLowPow *= zLow*zLow;
     intK += sign*(zUpPow-zLowPow)/((2*k+1)*fac(k)*fac(a+1-k));
   }

   intK *= fac(a)*acos(-1);

   return intK;
 }

 double integralProdKernel ( ) // double zc, double h, double zMax)
 {
   const unsigned a = 3; double zLow, zUp;
/*
   if (zc < h)
        zLow = -zc/h;
   else
*/
        zLow = -1;
/*
   if (zc+h > zMax)
        zUp = (zMax-zc)/h;
   else
*/
        zUp = 1;

   double intK = (zUp-zLow)/fac(a), zUpPow = zUp, zLowPow = zLow; 
   int sign = 1; 
   for (int k = 1; k < (a+1); k++)
   {
     sign *= -1;
//     zUpPow *= zUp*zUp; zLowPow *= zLow*zLow;
     intK += sign*(zUpPow-zLowPow)/((2*k+1)*fac(k)*fac(a-k));
   }

//   intK = acos(-1)/(a+1);

   intK *= fac(a)*acos(-1)/(a+1);

   return intK;
 }

 double integralProdEpanechnikov ( ) // double zc, double h, double zMax)
 {
   return 2.0*acos(-1)/3.0;
 }

 double integral3DEpanechnikov ( )
 {
   return 8.0*acos(-1)/15.0;
 }

 double integralGaussian ( )
 {
   return pow(2*acos(-1),1.5);
 }

}

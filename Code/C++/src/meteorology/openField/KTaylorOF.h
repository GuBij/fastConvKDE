#ifndef KTaylorOF_H
#define	KTaylorOF_H

 #include <iostream>
 #include "openField.h"

 class KTaylorOF
 :
 public openField
 {
   const string name_;
   double* eddyParams;
   double* w2bar;
   double* w2barGrad;
   double* tauL_w;
   double* tauL_u = NULL;
   double* corr;
   double* varianceRate;
   double* timeStep_;
   double* timeStep_u = NULL;
//   double* timeStepHor_;
//   int timeStep_;
   const double dtByTauL_;
   double minTimeStep_;
//   double tau_u;
//   double tau_v;
//  double tau_w;
//   const double minTauL_ = 4.0;
//   double sigma_u;
//   double sigma_v;
//   double sigma_w;
   double hABL;
   double heff_;
   double propConst;
   double L;
   double uStar;
   unsigned long metId;
   int mfreq;
   bool unstable;
   bool stable;

  protected:
   double sigma_u;
   double sigma_v;
   double sigma_w;

  public:
   KTaylorOF( const dictionary&, bool );
   void eddyDiff( double*, const double& );
   void eddyDiffRough ( double*, const double*, const double& );
   void stdevWind ( double*, const double& );
   void autoCorr( double*, const double& );
   void autoCorr( double*, const double&, const double& );
   void tauL ( double*, const double& );
   void updateParams( unsigned ); //ifstream* );
   void timeStep ( double*, const double& );
   void nonHomTurbCorr( double&, const double&, const double& );
//   double timeStep( const double& );

   inline string name()
   {
       return name_;
   }

   inline double minTimeStep()
   {
	return minTimeStep_;
   }

   inline double dtByTauL()
   {
	return dtByTauL_;
   }

   inline double sigmaV( unsigned meteoId )
   {
	return eddyParams[meteoId*8+1];
   }

   inline double sigmaU( unsigned meteoId )
   {
	return eddyParams[meteoId*8];
   }

   inline bool boolUnstable()
   {
	if ( L < 0 && L > -500 )
	  return true;

	return false;
   }

 };

#endif

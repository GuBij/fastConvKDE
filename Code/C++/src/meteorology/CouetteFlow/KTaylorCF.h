#ifndef KTaylorCF_H
#define	KTaylorCF_H

 #include <iostream>
 #include "CouetteFlow.h"

 class KTaylorCF
 :
 public CouetteFlow
 {
   const string name_;
   double* eddyParams;
   double* w2bar;
   double* tauL_w;
   double* corr;
   double* timeStep_;
   const double dtByTauL_;
   double minTimeStep_;
   double hABL;
   double heff_;
   double propConst;
   double L;
   double uStar;
   unsigned long metId;
   int mfreq;

  protected:
   double sigma_u;
   double sigma_v;
   double sigma_w;

  public:
   KTaylorCF( const dictionary&, bool );
   void eddyDiff( double*, const double*, const double& );
   void stdevWind ( double*, const double& );
   void autoCorr( double*, const double& );
   void autoCorr( double*, const double&, const double& );
   void tauL ( double*, const double& );
   void updateParams( unsigned ); 
   void timeStep ( double*, const double& );
   double nonHomTurbCorr( const double& , const double& );

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

   inline bool boolUnstable()
   {
        if ( L < 0 && L > -200 )
          return true;

        return false;
   }

 };

#endif

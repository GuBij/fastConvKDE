#ifndef UCouette_H
#define UCouette_H

#include "CouetteFlow.h"

 class UCouette
 :
   public CouetteFlow
 {
    const string name_;
    unsigned long metId;
    double* windParams;
    double* Umag_;
    double uStar;
    double Uref;
    double heff_;
    double windDirection;
    double sine;
    double cosine;
    double Psi0;

   public:
    UCouette(const dictionary&, bool);
    void U ( double*, const double& );
    void updateParams ( unsigned ); 
    inline string name()
    {
	return name_;
    }
    inline double wdir()
    {
	return windDirection;
    }
    inline double heff()
    {
	return heff_;
    }
    inline double* Umag()
    {
        return Umag_;
    }
 };

#endif

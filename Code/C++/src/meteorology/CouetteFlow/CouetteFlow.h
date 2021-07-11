#ifndef CouetteFlow_H
#define	CouetteFlow_H

 #include <fstream>
 #include <mathimf.h>
 #include "dictionary.h"
 #include "meteorology.h"

 class CouetteFlow
 :
 public dictionary, public meteorology
 {
   const string name_;
   const string meteofile_;
   const string humidityType_;
   const string windVarType_;
   const double g;
   const double zref_T;
   const double zref_T0;
   const double zref_wvar;
   const double pressure;
   const double rho;
   const double R;
   const double cp;
   const double d;
   const double flowRate;
   const double Ts;
   const bool useLangevin;
//   const unsigned latitude;
   bool humidityUniform_;
   bool windVarIsVelocity_;
   double humidityValue_;
   void picard(double, double, double, double&, double&, double&);
   void plumeRise(double, double, double, double, double&);

   protected:
    const double zref;
    const double kappa;
    const double z0;
    const double stackHeight;
    const double PI = acos(-1);
    const double fCori = 2*7.2921e-5*sin(58*PI/180);
    const unsigned noGridPointsZ;
    const double resZ;
    unsigned noRecords;
    double stabKernel(const double&, const string&, const bool&);
 
   public:
    CouetteFlow(const dictionary&, bool);
    CouetteFlow(const dictionary&, bool, bool);
    void calcParams(int,double*);
    inline void updateParams(ifstream*)
    {
	// not implemented;
    }

    inline virtual string name()
    {
       return name_;
    }

    inline string meteofile()
    {
       return meteofile_;
    }

    inline unsigned noRecords_()
    {
	return noRecords;
    }

    inline bool humidityUniform()
    {
	return humidityUniform_;
    }
 };

#endif

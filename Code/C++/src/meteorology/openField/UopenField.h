#ifndef UopenField_H
#define UopenField_H

#include "openField.h"

 class UopenField
 :
   public openField
 {
    const string name_;
    unsigned long metId;
    double* windParams;
    double* Umag_;
    double* windDirFluc_;
    double uStar_;
    double L;
    double heff_;
    double windDirection;
    double windDirNoise;
    double sigmaAzim;
    double sine;
    double cosine;
    double Psi0;
//    double elevation;

   public:
    UopenField(const dictionary&, bool);
    void U ( double*, const double& );
    void updateParams ( unsigned ); // ifstream* );
    inline void updateWindDirNoise( double noise )
    {
	 windDirNoise = noise;
    }
    inline void updateWindDirFluc()
    {
//	if ( sigmaAzim > minMeso )
//	{
	  windDirection = windParams[metId*5+3] + windDirNoise*sigmaAzim;
	  sine = sin(windDirection);
          cosine = cos(windDirection);
	  windDirFluc_[metId] = windDirNoise*sigmaAzim;
//	}
    }
    inline string name()
    {
	return name_;
    }
    inline double wdir( unsigned meteoId )
    {
	return windParams[meteoId*5+3] + windDirFluc_[meteoId];
    }
    inline double windDirFluc( unsigned meteoId )
    {
	return windDirFluc_[meteoId];
    }
    inline double mesoFlucSD( unsigned meteoId )
    {
	return windParams[meteoId*5+4];
    }
    inline void windDirFlucReset()
    {
	fill(windDirFluc_,windDirFluc_+noRecords,'\0');
    }
    inline double heff()
    {
	return heff_;
    }
    inline double Umag( const double& z )
    {
	unsigned cellIndex = z/resZ;
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;

        return Umag_[cellIndex];
    }
    inline double uStarSqr()
    {
	return uStar_*uStar_;
    }
 };

#endif

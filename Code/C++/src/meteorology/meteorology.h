#ifndef meteorology_H
#define meteorology_H

#include <iostream>
#include "dictionary.h"

 class meteorology
 {
     string type;
   public:
     meteorology ();
     static meteorology* New(bool, bool);
     static meteorology* New(string, bool, bool);
     virtual string name() =0;
     virtual string meteofile() =0;
     virtual unsigned noRecords_() =0;
     virtual void U(double* u, const double& z)
     {
        cerr << "Error: function 'U ( double* , const double& )' not implemented in class " << this->name() << endl;
     }

     virtual inline double Umag( const double& z )
     {
        cerr << "Error: function 'Umag ( const double& )' not implemented in class " << this->name() << endl;
     }

     virtual inline double uStarSqr()
     {
	cerr << "Error: function 'uStarSqr ( )' not implemented in class " << this->name() << endl;
     }

     virtual void eddyDiff(double* eddyVec, const double& z)
     {
        cerr << "Error: function 'eddyDiff ( double* , const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void eddyDiffRough ( double* eddyVec, const double* tauLVec, const double& z)
     {
        cerr << "Error: function 'eddyDiffRough ( double* , const double*, const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void stdevWind ( double* stdevVec, const double& z)
     {
	cerr << "Error: function 'stdevWind ( double* , const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void GradK(double* k, const double& z)
     {
	k[0]=0; k[1]=0; k[2]=0;
     }

     virtual void updateTimeStep(double newTimeStep, double newTimeStepVert)
     {
	cerr << "Error: function 'updateTimeStep ( double, double )' not implemented in class " << this->name() << endl;
     }

     virtual void nonHomTurbCorr(double& a, const double& z, const double& zdot )
     {
        cerr << "Error: function 'nonHomTurbCorr ( double&, const double&, const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void updateParams( unsigned ) //ifstream* ifile)
     {
	cerr << "Error: function 'updateParams ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual void updateParams( unsigned, double ) //ifstream* ifile)
     {
        cerr << "Error: function 'updateParams ( unsigned, double )' not implemented in class " << this->name() << endl;
     }

     virtual inline void updateWindDirNoise ( double )
     {
        cerr << "Error: function 'updateWindDirNoise ( double )' not implemented in class " << this->name() << endl;
     }

     virtual inline void updateWindDirFluc ()
     {
	cerr << "Error: function 'updateWindDirFluc ( double )' not implemented in class " << this->name() << endl;
     }

     virtual void calcParams(int timeIndex, double* data)
     {
	cerr << "Error: function 'calcParams ( int, double* )' not implemented in class " << this->name() << endl;
     }

     virtual void tauL (double* tauVec, const double& z)
     {
	cerr << "Error: function 'tauL ( double* , const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void autoCorr( double* corrVec, const double& z )
     {
        cerr << "Error: function 'autoCorr ( double* , const double& )' not implemented in class " << this->name() << endl;
     }

     virtual void autoCorr( double* corrVec, const double& z, const double& timeStepVert )
     {
	cerr << "Error: function 'autoCorr ( double* , const double&, const double& )' not implemented in class " << this->name() << endl;
     }

     virtual double wdir( unsigned meteoId )
     {
	cerr << "Error: function 'wdir ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual double windDirFluc( unsigned meteoId )
     {
        cerr << "Error: function 'windDirFluc ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual void windDirFlucReset()
     {
        cerr << "Error: function 'windDirFluc ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual double mesoFlucSD( unsigned meteoId )
     {
        cerr << "Error: function 'mesoFlucSD ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual double heff()
     {
        cerr << "Error: function 'heff ( )' not implemented in class " << this->name() << endl;
     }

     virtual void timeStep( double* timeStepVec, const double& z )
     {
        cerr << "Error: function 'timeStep ( double*, const double& )' not implemented in class " << this->name() << endl;
     }

     virtual inline double minTimeStep()
     {
        cerr << "Error: function 'minTimeStep ( )' not implemented in class " << this->name() << endl;
     }

     virtual inline bool boolUnstable()
     {
	cerr << "Error: function 'boolUnstable ( )' not implemented in class " << this->name() << endl;
     }

     virtual void reset()
     {
	cerr << "Error: function 'reset ( )' not implemented in class " << this->name() << endl;
     }

     virtual inline bool humidityUniform()
     {
	cerr << "Error: function 'humidityUniform ( )' not implemented in class " << this->name() << endl;
     }

     virtual inline double sigmaV( unsigned meteoId )
     {
	cerr << "Error: function 'sigmaV ( unsigned )' not implemented in class " << this->name() << endl;
     }

     virtual inline double sigmaU( unsigned meteoId )
     {
        cerr << "Error: function 'sigmaU ( unsigned )' not implemented in class " << this->name() << endl;
     }
 };

#endif

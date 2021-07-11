#include "UopenField.h"

 UopenField::UopenField(const dictionary& dict, bool master_proc )
 :
 openField(dict, master_proc),
 windDirNoise(0),
 name_("UopenField")
 {
    windParams = new double[noRecords*5];
    ifstream windParamFile("output/PARAMS_" + name_ + "_" + dict.get(0));
    for (unsigned i = 0; i < noRecords; i++)
    {
	unsigned timeIndex, ii = 5*i;
	windParamFile >> timeIndex >> windParams[ii] >> windParams[ii+1] >> windParams[ii+2] >> windParams[ii+3] >> windParams[ii+4];
	windParams[ii+3] *= PI/180; windParams[ii+3] += PI;
	windParams[ii+4] *= PI/180;
    }
    windParamFile.close();

    double Psi0;
    Umag_ = new double[noRecords*noGridPointsZ];
    for ( unsigned long i = 0; i < (noRecords*noGridPointsZ); i++ )
    {
	unsigned long ii = i % noGridPointsZ, r = i/noGridPointsZ;
	double zi = resZ*(ii+0.5);

	if ( ii == 0 )
	{
	  updateParams(r); Psi0 = openField::stabKernel(z0/L,"M",true); 
	}

	if (zi < 30*z0)
	  Umag_[i] = uStar_/kappa*(log(30)-openField::stabKernel(30*z0/L,"M",true)+Psi0)/(30*z0)*zi;
	else
	  Umag_[i] = uStar_/kappa*(log(zi/z0)-openField::stabKernel(zi/L,"M",true)+Psi0);
    }
    windDirFluc_ = new double[noRecords]();
    updateParams(0);
 }

/*
 double UopenField::U( const double& z )
 {
     double umag;

     if (z < 2*z0)
       umag = uStar_/kappa*(log(2)-openField::stabKernel(2*z0/L,"M",true)+Psi0)/(2*z0)*z;
     else
        umag = uStar_/kappa*(log(z/z0)-openField::stabKernel(z/L,"M",true)+Psi0);

     return umag;
 }
*/

 void UopenField::U(double* Uvec, const double& z)
 {

	unsigned cellIndex = z/resZ;
	if (cellIndex > (noGridPointsZ-1))
	   cellIndex = noGridPointsZ-1;

	cellIndex += metId*noGridPointsZ;
        *Uvec = Umag_[cellIndex]*cosine; //(windDirection); //*cos(elevation);
        *(Uvec + 1) = Umag_[cellIndex]*sine; //(windDirection); //*cos(elevation);
        *(Uvec + 2) = 0.0;
 }

 void UopenField::updateParams( unsigned meteo_id ) //ifstream* ifile)
 {
	if (meteo_id < noRecords)
	{
	   uStar_ = windParams[meteo_id*5];
	   L = windParams[meteo_id*5+1];
	   heff_ = windParams[meteo_id*5+2];
	   windDirection = windParams[meteo_id*5+3];
	   sigmaAzim = windParams[meteo_id*5+4];
/*
	   if ( noise != 1000 )
	   {
		windDirNoise = noise; cout << "\nnoise wdir: " << noise << ", sigmaAzim: " << sigmaAzim << ", wdir fluc: " << windDirNoise*sigmaAzim << endl;
	   }
*/
	   windDirection += windDirNoise*sigmaAzim;
	   sine = sin(windDirection); 
	   cosine = cos(windDirection);
	   metId = meteo_id;
	}
	else
	{
	   cout << "\n\t ERROR: meteo record with id " << meteo_id << " does not exist.\n" << endl;
	   exit(EXIT_FAILURE);	   
	}
 }


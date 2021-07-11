#include "UCouette.h"

 UCouette::UCouette(const dictionary& dict, bool master_proc )
 :
 CouetteFlow(dict, master_proc),
 name_("UCouette")
 {
    windParams = new double[noRecords*4];
    ifstream windParamFile("output/PARAMS_" + name_ + "_" + dict.get(0));
    for (unsigned i = 0; i < noRecords; i++)
    {
	unsigned timeIndex, ii = 4*i;
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
	  updateParams(r); 
	}

	Umag_[i] = Uref; // *zi/zref; 
    }

    updateParams(0);
 }

 void UCouette::U(double* Uvec, const double& z)
 {

	unsigned cellIndex = floor(z/resZ);
	if (cellIndex > (noGridPointsZ-1))
	   cellIndex = noGridPointsZ-1;

	cellIndex += metId*noGridPointsZ;
        *Uvec = Umag_[cellIndex]*cosine; //(windDirection); //*cos(elevation);
        *(Uvec + 1) = Umag_[cellIndex]*sine; //(windDirection); //*cos(elevation);
        *(Uvec + 2) = 0.0;
 }

 void UCouette::updateParams( unsigned meteo_id ) //ifstream* ifile)
 {
	if (meteo_id < noRecords)
	{
	   uStar = windParams[meteo_id*4];
	   Uref = windParams[meteo_id*4+1];
	   heff_ = windParams[meteo_id*4+2];
	   windDirection = windParams[meteo_id*4+3];
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


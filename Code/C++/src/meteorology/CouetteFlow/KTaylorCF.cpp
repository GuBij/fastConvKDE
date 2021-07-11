#include "KTaylorCF.h"

 KTaylorCF::KTaylorCF(const dictionary& dict, bool master_proc )
 :
 CouetteFlow(dict,true, master_proc),
 name_("KTaylorCF"),
 dtByTauL_(dict.getfl(11)),
// minTimeStep_(dict.getfl(12)),
 mfreq(dict.getfl(7))
 {
    eddyParams = new double[noRecords*8];
    ifstream eddyParamFile("output/PARAMS_" + name_ + "_" + dict.get(0));
    for (unsigned i = 0; i < noRecords; i++)
    {
        unsigned timeIndex, ii = 8*i;
        eddyParamFile >> timeIndex >> eddyParams[ii] >> eddyParams[ii+1] >> eddyParams[ii+2] >> eddyParams[ii+3] >> eddyParams[ii+4] >> eddyParams[ii+5] >> eddyParams[ii+6] >> eddyParams[ii+7];
    }
    eddyParamFile.close();

    w2bar = new double[noRecords*noGridPointsZ];
    tauL_w = new double[noRecords*noGridPointsZ]();
    corr = new double[2*noRecords*noGridPointsZ]();
    timeStep_ = new double[noRecords*noGridPointsZ](); 
    double pct, lengthScale, tauL_uv = 0.0;
    bool unstable, neutral, updateTimeStep, update = true;
    unsigned heffCell;
    unsigned long i = 0;
    while ( i < (noRecords*noGridPointsZ) )
    {
        unsigned long ii = i % noGridPointsZ, r = i/noGridPointsZ;

        if ( update )
        {
          updateParams(r); updateTimeStep = true; update = false; pct = dtByTauL_; heffCell = heff_/resZ;
	  if ( heffCell > (noGridPointsZ-1) )
	     heffCell = noGridPointsZ-1;

	  i += heffCell; ii = heffCell;
	  minTimeStep_ = pct*0.5*30*z0/sigma_w/(1+15*fCori*30*z0/uStar);
        }
	double zi = resZ*(ii+0.5);

	w2bar[i] = sigma_w; // *pow(zi/zref,1.0/3.0);
        tauL_w[i] = (0.5*heff_/sigma_w)/(1+15*fCori*heff_/uStar); // /(30*fCori*sigma_w); 

	double timeStepi; 
	double tauL_wi;
	tauL_wi = tauL_w[i];
	timeStepi = pct*tauL_wi;

	if ( timeStepi < minTimeStep_ )
	{
           timeStep_[i] = minTimeStep_;
//	   tauL_wi = minTimeStep_/pct;
//	   tauL_w[i] = tauL_wi;
	}
	else
	   timeStep_[i] = timeStepi;

	if ( ii == heffCell )
        {
	  timeStep_[i] = mfreq/ceil(mfreq/timeStep_[i]); pct = timeStep_[i]/tauL_wi; timeStepi = timeStep_[i];
        }

	corr[i] = exp(-timeStep_[i]/tauL_wi);
	corr[i+noRecords*noGridPointsZ] = corr[i];

	i++;

	if ( ii == heffCell && heffCell == 0 )
        {
	  updateTimeStep = false;
	}
	else if ( ii == heffCell )
        {
	  i = r*noGridPointsZ;
	}
	else if ( heffCell > 0 && ii == heffCell - 1 ) 
        {
          i++; updateTimeStep = false;
        }

	if ( i == (r+1)*noGridPointsZ )
	{
	  update = true;
	}
    } 

    updateParams(0);
 }

 void KTaylorCF::tauL ( double* tauVec, const double& z )
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;
        *tauVec = tauL_w[cellIndex];
        *(tauVec+1) = tauL_w[cellIndex];
	*(tauVec+2) = tauL_w[cellIndex];
 }

 void KTaylorCF::eddyDiff ( double* eddyVec, const double* tauLVec, const double& z)
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

	cellIndex += metId*noGridPointsZ;

        *eddyVec=sqrt(2/tauLVec[0])*sigma_u;
        *(eddyVec+1)=sqrt(2/tauLVec[1])*sigma_v;
        *(eddyVec+2)=sqrt(2/tauLVec[2])*w2bar[cellIndex];
 }

 void KTaylorCF::stdevWind ( double* stdevVec, const double& z)
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;

	*stdevVec=sigma_u;
	*(stdevVec+1)=sigma_v;
	*(stdevVec+2)=w2bar[cellIndex];
 }

 double KTaylorCF::nonHomTurbCorr ( const double& z, const double& zdot )
 {
	return 0.0; 
 }

 void KTaylorCF::timeStep( double* timeStepVec, const double& z )
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

	cellIndex += metId*noGridPointsZ;

	*(timeStepVec+1) = timeStep_[cellIndex];
	*timeStepVec = timeStep_[cellIndex];
 }

 void KTaylorCF::autoCorr( double* corrVec, const double& z )
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;
	corrVec[0] = corr[cellIndex];
	corrVec[1] = corr[cellIndex+noRecords*noGridPointsZ];
 }

 void KTaylorCF::autoCorr( double* corrVec, const double& z, const double& timeStepVert ) //doesn't work yet for stable
 {
        unsigned cellIndex = floor(z/resZ);
        if (cellIndex > (noGridPointsZ-1))
           cellIndex = noGridPointsZ-1;

        cellIndex += metId*noGridPointsZ;

        corrVec[0] = corr[cellIndex];
        corrVec[1] = exp(-timeStepVert/tauL_w[cellIndex]);
 }

 void KTaylorCF::updateParams ( unsigned meteo_id ) 
 {
        if (meteo_id < noRecords)
        {
           sigma_u = eddyParams[meteo_id*8];
           sigma_v = eddyParams[meteo_id*8+1];
           sigma_w = eddyParams[meteo_id*8+2];
	   propConst = eddyParams[meteo_id*8+3];
           uStar = eddyParams[meteo_id*8+4];
           L = eddyParams[meteo_id*8+5];
	   heff_ = eddyParams[meteo_id*8+6];
	   hABL = eddyParams[meteo_id*8+7];
	   metId = meteo_id; 
       }
        else
        {
           cout << "\n\t ERROR: meteo record with id " << meteo_id << " does not exist.\n" << endl; 
           exit(EXIT_FAILURE);
        }	
 }

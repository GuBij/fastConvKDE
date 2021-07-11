#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sys/stat.h>
#include <mpi.h>
#include <mathimf.h>
#include "mkl_vsl.h"
#include "meteorology.h"
#include "densityKernel.h"

 int main(int argc, char **argv)
 {
   bool LangevinModel = false;
   bool stationary = false;
   bool prodKernel = true;
   unsigned long noStations;
   unsigned seedNr = 1;
   double* stations;  
   {
     int i=1;
     while (  argc > 1 && i < argc )
     {
        if ( string(argv[i]) == "-Langevin" )
	{
	    LangevinModel = true;
	}
	else if ( string(argv[i]) == "-stationary" )
	{
	    stationary = true;
	}
	else if ( string(argv[i]) == "-seed" )
	{
	    i++;
	    stringstream ss(argv[i]);
	    ss >> seedNr;
	}
	else
	{
	    cout << "\n\t ERROR: option '" << argv[i] << "' is not implemented. Valid options are: \n" << endl;
	    cout << "\t-Langevin\n\t-stationary" << endl;
	    exit(EXIT_FAILURE);
	}

	i++;
     }

     ifstream ifile("stationCoord.txt");
     ifile.unsetf(ios_base::skipws);
     if (ifile.is_open())
       noStations = count(istream_iterator<char>(ifile),istream_iterator<char>(),'\n');
     else
     {
       cerr << "\nfile 'stationsCoord.txt' not found.\n" << endl;
       exit(EXIT_FAILURE);
     }
     ifile.setf(ios_base::skipws);
     ifile.close();

     stations = new double[3*noStations]();
     ifile.open("stationCoord.txt", ifstream::in);
     for (unsigned long j=0; j<noStations; j++)
     {
       ifile >> stations[3*j] >> stations[3*j+1] >> stations[3*j+2];
     }
     ifile.close();
   }

   int noProcs, procid;
   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &noProcs);
   MPI_Comm_rank(MPI_COMM_WORLD, &procid);

   int master = noProcs - 1; 
   bool isMaster = false;
   if ( procid == master )
	isMaster = true;

   meteorology* windModel(meteorology::New("windModel", LangevinModel, isMaster));
   meteorology* eddyDiffModel(meteorology::New("eddyDiffModel", LangevinModel, isMaster));

   dictionary dict({"stackHeight","releaseDuration","measureFreq","noPartPerRel","xMax","yMax","zMax","sourceStrength"});

   unsigned long ppr = dict.getfl(3); 
   int mfreq = dict.getfl(2);
   double xMax = dict.getfl(4), yMax = dict.getfl(5), zMax = dict.getfl(6);
   ppr = ceil(double(ppr)/noProcs); 

   if ( isMaster ) 
   {
	if ( LangevinModel ) { cout << "\nLangevin model selected." << endl; }
	if ( stationary ) { cout << "\nStationary release assumed." << endl; }
	cout << "\n - domain: [-" << xMax  << " m, " << xMax << " m ] x [-" << yMax << " m, " << yMax << " m ] x [0 m, " << zMax << " m ]" << endl;
	cout << "\n - source strength: " << dict.getfl(7) << " Bq/h " << endl; 
	cout << "\n - stack height: " << dict.getfl(0) << " m" << endl;
	cout << "\n - " << ppr*noProcs << " particles per time step are released" << endl;
	cout << "\n - concentration will be calculated for " << noStations << " stations" << endl;
   }

   const int nor = windModel->noRecords_();
   double endTime;    
   bool warmingUp = false;
   double dt[2]; eddyDiffModel->timeStep(dt,windModel->heff());
   if ( stationary )
   {
	endTime = nor*mfreq;
	warmingUp = true;
   }
   else if ( dict.getfl(1) == 0 ) 
	endTime = dt[0];
   else
	endTime = dict.getfl(1);

   const double Pi = acos(-1);
   const double Q = dict.getfl(7)*1000; // /3600;
   unsigned Nt = ceil(double(mfreq)/eddyDiffModel->minTimeStep());
   unsigned long pi = 0;
   int warmingUpPeriod = 0;
   bool unstableTrue = eddyDiffModel->boolUnstable();
   long tm = mfreq;
   double releaseTime = 0;
   double zt[2] = {windModel->heff(), 0.0}, zdot[2]={};
   const double x0[2] = {}, xdot_eject[3] = {}; 
   double tL[3], beta[3], corr[2], pathIntU[2], u0[3] = {}, sigma_u[3]; 
   double var_xy[2]={}, integralCorr[2]={}, var_xy_dummy[9]={}; var_xy_dummy[8] = 1.0;

   unsigned seed_(noProcs*seedNr + procid);
   VSLStreamStatePtr stream;
   vslNewStream( &stream, VSL_BRNG_MT19937, seed_);
   {
     double xsi[4];
     const int status = vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, 4, xsi, 0.0, 1.0);
     if (status != VSL_ERROR_OK)
     {
       cout << "\n\t ERROR: sampling by function 'vdRngGaussian' failed; status flag: " << status << "\n" << endl;
       exit(EXIT_FAILURE);
     }
     eddyDiffModel->stdevWind(sigma_u, zt[0]);

     u0[0] = sigma_u[0]*xsi[0];
     u0[1] = sigma_u[1]*xsi[1];
     u0[2] = sigma_u[2]*xsi[2]; zdot[0] = u0[2] + xdot_eject[2];
     windModel->updateWindDirNoise(xsi[3]);
     windModel->updateWindDirFluc();
   }

   double Uvec0[3]; windModel->U(Uvec0,zt[0]);
   double tp = dt[0], tp_z = 0; // dt[1];
   double dQ = Q*dt[1]/(ppr*noProcs); //*dt[0]
   pathIntU[0] = Uvec0[0]*dt[0]; //dt[1]
   pathIntU[1] = Uvec0[1]*dt[0]; //dt[1]
   eddyDiffModel->tauL(tL,zt[0]);
   eddyDiffModel->eddyDiff(beta,zt[0]);
   eddyDiffModel->autoCorr(corr,zt[0]); // ,dt[0]); // dt[1]);
   tL[0] *= -1; tL[1] *= -1; tL[2] *= -1;

   double tLbetaSqr, corrSqr;

     tLbetaSqr=beta[0]*beta[0]*tL[0], corrSqr = corr[0]*corr[0];
     var_xy[0] += (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[0]+dt[0])*tLbetaSqr*tL[0];
     var_xy[0] += 2*(corr[0]-1)*tL[0]*var_xy_dummy[1]+(corr[0]-1)*(corr[0]-1)*tL[0]*tL[0]*var_xy_dummy[0];
     integralCorr[0] += (corr[0]-1)*tL[0]*var_xy_dummy[8];

     var_xy_dummy[0] = corrSqr*(var_xy_dummy[0]+0.5*tLbetaSqr)-0.5*tLbetaSqr; //Itilde
     var_xy_dummy[2] = corrSqr*(var_xy_dummy[2]+var_xy_dummy[3]); //Ihat3
     var_xy_dummy[3] = -0.5*tLbetaSqr*(1-corrSqr);
     var_xy_dummy[1] = corr[0]*var_xy_dummy[1]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[0]*0.5+(1-1/corr[0])*var_xy_dummy[2]*tL[0]; //Ihat1

     tLbetaSqr=beta[1]*beta[1]*tL[1];
     var_xy[1] += (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[1]+dt[0])*tLbetaSqr*tL[1];
     var_xy[1] += 2*(corr[0]-1)*tL[1]*var_xy_dummy[5]+(corr[0]-1)*(corr[0]-1)*tL[1]*tL[1]*var_xy_dummy[4];
     integralCorr[1] += (corr[0]-1)*tL[1]*var_xy_dummy[8];

     var_xy_dummy[4] = corrSqr*(var_xy_dummy[4]+0.5*tLbetaSqr)-0.5*tLbetaSqr;
     var_xy_dummy[6] = corrSqr*(var_xy_dummy[6]+var_xy_dummy[7]);
     var_xy_dummy[7] = -0.5*tLbetaSqr*(1-corrSqr);
     var_xy_dummy[5] = corr[0]*var_xy_dummy[5]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[1]*0.5+(1-1/corr[0])*var_xy_dummy[6]*tL[1];
     var_xy_dummy[8] *= corr[0];


   double* meanDrift = new double[4*nor]();
   double* positions = new double[3*ppr*nor]();
   double* pathVar = new double[2*ppr*nor]();
   double* C = new double[nor*noStations]();

   do {
	double* xsi = new double[Nt+5]();
	const int status = vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, Nt+5, xsi, 0.0, 1.0);
	if (status != VSL_ERROR_OK)
	{
	  cout << "\n\t ERROR: sampling by function 'vdRngGaussian' failed; status flag: " << status << "\n" << endl;
	  exit(EXIT_FAILURE);
	}

	double A0Scheme[2]={}, tL0[3] = {}, sigma_u0[3] = {}, beta0[3] = {};
	unsigned kk = 0;
	bool updateTimeStep = true; 
	eddyDiffModel->tauL(tL0, 0); eddyDiffModel->stdevWind(sigma_u0, 0); eddyDiffModel->eddyDiff(beta0,0);
	while ( ceil(tp*1000)/1000 < tm || tp_z < tm ) 
	{
	    double A0; eddyDiffModel->nonHomTurbCorr(A0,zt[0],zdot[0]); 
	    if ( unstableTrue )
	    {
		zt[1] = zt[0] + zdot[0]*dt[1];
		zdot[1] = corr[1]*(zdot[0] + A0*dt[1] + beta[2]*sqrt(dt[1])*xsi[kk]); // (1+dt[1]/tL[2])
	    }
	    else
	    {
              const double diff_z = sqrt(0.5*((corr[1]-2)*(corr[1]-2)-1)*tL[2]+dt[1])*(-tL[2]);
              const double diff_zdot = sqrt(0.5*(corr[1]*corr[1]-1)*tL[2]);

              zt[1] = zt[0] + A0Scheme[0]*A0 - (1-corr[1])*tL[2]*(zdot[0] + A0*tL[2]) - dt[1]*A0*tL[2]+ diff_z*beta[2]*xsi[kk]; 
              zdot[1] = corr[1]*(zdot[0] + A0Scheme[1]*A0) - A0*(1-corr[1])*tL[2] + beta[2]*diff_zdot*xsi[kk]; 
	    }

            if ( ceil(tp*1000)/1000 < tm && ((( tp_z <= tp ) && ( tp_z+dt[1] >= tp )) || ( tp < tp_z ))) 
            {
                int ii = 1;
                if ( ( abs(tp - tp_z) < tp_z + dt[1] - tp ) || ( tp < tp_z ))
                  ii--;

                double zz = zt[ii] + zdot[0]*(tp - tp_z - ii*dt[1]); 
            	if ( zz < 0 && unstableTrue )
            	{
                  const double t0 = -zt[0]/zdot[0], varGround = sigma_u0[2]*sigma_u0[2]*(1-exp(-2*(dt[1]-t0)/tL0[2]));
                  double zzdot = -2*varGround*log(1-exp(-0.5*zdot[0]*zdot[0]/varGround)); zzdot = sqrt(zzdot);
                  zz = (dt[1]-t0)*zzdot;
            	}
            	else if ( zz < 0)
            	  zz *= -1;

		eddyDiffModel->timeStep(dt,zz);
		eddyDiffModel->autoCorr(corr,zz);
		eddyDiffModel->tauL(tL,zz);
		eddyDiffModel->eddyDiff(beta,zz);
		tL[0] *= -1; tL[1] *= -1; tL[2] *= -1;

		tLbetaSqr=beta[0]*beta[0]*tL[0], corrSqr = corr[0]*corr[0];
		var_xy[0] += (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[0]+dt[0])*tLbetaSqr*tL[0];
		var_xy[0] += 2*(corr[0]-1)*tL[0]*var_xy_dummy[1]+(corr[0]-1)*(corr[0]-1)*tL[0]*tL[0]*var_xy_dummy[0];
		integralCorr[0] += (corr[0]-1)*tL[0]*var_xy_dummy[8];

		var_xy_dummy[0] = corrSqr*(var_xy_dummy[0]+0.5*tLbetaSqr)-0.5*tLbetaSqr; //Itilde
		var_xy_dummy[2] = corrSqr*(var_xy_dummy[2]+var_xy_dummy[3]); //Ihat3
		var_xy_dummy[3] = -0.5*tLbetaSqr*(1-corrSqr);
		var_xy_dummy[1] = corr[0]*var_xy_dummy[1]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[0]*0.5+(1-1/corr[0])*var_xy_dummy[2]*tL[0]; //Ihat1

		tLbetaSqr=beta[1]*beta[1]*tL[1];
		var_xy[1] += (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[1]+dt[0])*tLbetaSqr*tL[1];
		var_xy[1] += 2*(corr[0]-1)*tL[1]*var_xy_dummy[5]+(corr[0]-1)*(corr[0]-1)*tL[1]*tL[1]*var_xy_dummy[4];
		integralCorr[1] += (corr[0]-1)*tL[1]*var_xy_dummy[8];

		var_xy_dummy[4] = corrSqr*(var_xy_dummy[4]+0.5*tLbetaSqr)-0.5*tLbetaSqr;
		var_xy_dummy[6] = corrSqr*(var_xy_dummy[6]+var_xy_dummy[7]);
		var_xy_dummy[7] = -0.5*tLbetaSqr*(1-corrSqr);
		var_xy_dummy[5] = corr[0]*var_xy_dummy[5]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[1]*0.5+(1-1/corr[0])*var_xy_dummy[6]*tL[1];
		var_xy_dummy[8] *= corr[0];

		tp += dt[0];

		double Uvec[3];
            	windModel->U(Uvec,zz);
            	pathIntU[0] += Uvec[0]*dt[0];
            	pathIntU[1] += Uvec[1]*dt[0];
            }

	    const double zdot1 = zdot[0];
	    if ( unstableTrue && zt[1] < 0 )
	    {
                const double t0 = -zt[0]/zdot[0], varGround = sigma_u0[2]*sigma_u0[2]*(1-exp(-2*(dt[1]-t0)/tL0[2]));
                zdot[0] = -2*varGround*log(1-exp(-0.5*zdot[0]*zdot[0]/varGround)); zdot[0] = sqrt(zdot[0]);
                zt[1] = (dt[1]-t0)*zdot[0];
		zdot[1] = exp(-(dt[1]-t0)/tL0[2])*(zdot[0] + beta0[2]*sqrt(dt[1]-t0)*xsi[kk]); 
	    }
            else if (zt[1] < 0)
            {
                zt[1] *= -1; zdot[0] *= -1; zdot[1] *= -1; A0Scheme[0] *= -1; A0Scheme[1] *= -1;
            }
	    const double zt0 = zt[1], zdot0 = zdot[0];
            zt[0] = zt[1];
            zdot[0] = zdot[1];
	    tp_z += dt[1];

	    if ( unstableTrue && tp_z > tm )
	    {
		zt[0] += zdot0*(tm - tp_z); // zdot0 is reflected velocity 
		if ( zt[0] < 0 )
		{
		  const double t0 = tp_z-zt0/zdot0; 
                  zt[0] = (tm-t0)*zdot1; // zdot1 is inicident velocity
		}
		tp_z = tm;
	    }
	    else if ( tp_z > tm )
	    {
                zt[0] += zdot0*(tm - tp_z); 
                tp_z = tm;
		if ( zt[0] < 0 )
		{
		  zt[0] *= -1; zdot[0] *= -1;
		}
	    }
	    else
	    {
	      eddyDiffModel->autoCorr(corr,zt[0]);
	      eddyDiffModel->tauL(tL,zt[0]);
              eddyDiffModel->eddyDiff(beta,zt[0]);
	      tL[2] *= -1;
	    }
	    eddyDiffModel->timeStep(dt,zt[0]);
	    kk++;
	}

       double argExp[2] = {};
       {
	 double windAngle = windModel->wdir(tm/mfreq-1), cosine = cos(windAngle), sine = sin(windAngle);
         argExp[0] = -pathIntU[0]-x0[0] - integralCorr[0]*(cosine*(u0[0] + xdot_eject[0]) - sine*(u0[1] + xdot_eject[1])); 
         argExp[1] = -pathIntU[1]-x0[1] - integralCorr[1]*(sine*(u0[0] + xdot_eject[0]) + cosine*(u0[1] + xdot_eject[1])); 
       }

	if ( !warmingUp )
	{
           unsigned long i = 4*(tm/mfreq-1), ii = ((tm/mfreq-1)*ppr + (pi % ppr))*3;
           meanDrift[i+2] += zt[0]; meanDrift[i+3] += 1;
           positions[ii] = argExp[0]; positions[ii+1] = argExp[1]; positions[ii+2] = zt[0];
	   pathVar[(2*ii)/3] = var_xy[0]; 
	   pathVar[(2*ii)/3+1] = var_xy[1]; 
	}

	const bool death = (tm/mfreq == nor) || (abs(argExp[0]) > xMax) || (abs(argExp[1]) > yMax ) || (abs(argExp[2]) > zMax);
	bool newPuff = false;

	if ( (warmingUp || (releaseTime < 0)) && death )
	{
          if ( ( pi < ppr ) && ( tm > warmingUpPeriod ) )
          {
            warmingUpPeriod = tm;
          }

	  tm = releaseTime/mfreq;
	  pi++;
	  if ( pi == ppr )
	  {
	    releaseTime -= warmingUpPeriod;
	    if ( isMaster ) { cout << "\nWarming up period: " << warmingUpPeriod << " s" << endl; }
	  }

	  windModel->updateParams(0); eddyDiffModel->updateParams(0);
	  eddyDiffModel->timeStep(dt, windModel->heff());
	  if (pi % ppr == 0)
	  {
	    newPuff = true;
	    releaseTime += dt[0];
	  }
	  if ((pi >= ppr) && (releaseTime < 0))
	     warmingUp = true;
	  else if ((pi >= ppr) && (releaseTime == 0))
	  {
	     warmingUp = false; tm++; 
	  }
	  else if (pi < ppr)
          {
            tm++;
          }
	  tm *= mfreq;
	  tp = releaseTime;
	  tp_z = releaseTime;
	  zt[0] = windModel->heff();

          windModel->U(Uvec0,zt[0]);
          pathIntU[0] = Uvec0[0]*dt[0];
          pathIntU[1] = Uvec0[1]*dt[0];
          eddyDiffModel->tauL(tL,zt[0]);
          eddyDiffModel->eddyDiff(beta,zt[0]);
          eddyDiffModel->autoCorr(corr,zt[0]); //,dt[1]);
          tL[0] *= -1; tL[1] *= -1; tL[2] *= -1;

	  unstableTrue = eddyDiffModel->boolUnstable();
	  eddyDiffModel->stdevWind(sigma_u,zt[0]);
	  u0[0] = sigma_u[0]*xsi[kk]; 
	  u0[1] = sigma_u[1]*xsi[kk+1]; 
	  u0[2] = sigma_u[2]*xsi[kk+2]; zdot[0] = u0[2] + xdot_eject[2]; 
	  eddyDiffModel->tauL(tL0, 0); 
	  eddyDiffModel->stdevWind(sigma_u0, 0);
	  eddyDiffModel->eddyDiff(beta0,0);

          for (unsigned i=0; i < 8; i++)
          {
            var_xy_dummy[i]=0;
          }
          var_xy_dummy[8] = 1.0;

            tLbetaSqr=beta[0]*beta[0]*tL[0], corrSqr = corr[0]*corr[0];
            var_xy[0] = (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[0]+dt[0])*tLbetaSqr*tL[0];
            var_xy[0] += 2*(corr[0]-1)*tL[0]*var_xy_dummy[1]+(corr[0]-1)*(corr[0]-1)*tL[0]*tL[0]*var_xy_dummy[0];
            integralCorr[0] = (corr[0]-1)*tL[0]*var_xy_dummy[8];

            var_xy_dummy[0] = corrSqr*(var_xy_dummy[0]+0.5*tLbetaSqr)-0.5*tLbetaSqr; //Itilde
            var_xy_dummy[2] = corrSqr*(var_xy_dummy[2]+var_xy_dummy[3]); //Ihat3
            var_xy_dummy[3] = -0.5*tLbetaSqr*(1-corrSqr);
            var_xy_dummy[1] = corr[0]*var_xy_dummy[1]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[0]*0.5+(1-1/corr[0])*var_xy_dummy[2]*tL[0]; //Ihat1

            tLbetaSqr=beta[1]*beta[1]*tL[1];
            var_xy[1] = (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[1]+dt[0])*tLbetaSqr*tL[1];
            var_xy[1] += 2*(corr[0]-1)*tL[1]*var_xy_dummy[5]+(corr[0]-1)*(corr[0]-1)*tL[1]*tL[1]*var_xy_dummy[4];
            integralCorr[1] = (corr[0]-1)*tL[1]*var_xy_dummy[8];

            var_xy_dummy[4] = corrSqr*(var_xy_dummy[4]+0.5*tLbetaSqr)-0.5*tLbetaSqr;
            var_xy_dummy[6] = corrSqr*(var_xy_dummy[6]+var_xy_dummy[7]);
            var_xy_dummy[7] = -0.5*tLbetaSqr*(1-corrSqr);
            var_xy_dummy[5] = corr[0]*var_xy_dummy[5]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[1]*0.5+(1-1/corr[0])*var_xy_dummy[6]*tL[1];
            var_xy_dummy[8] *= corr[0];

	    tp += dt[0];
	}
	else if (death)
	{
          pi++;
          tm = releaseTime/mfreq+1;
	  windModel->updateParams(tm-1); eddyDiffModel->updateParams(tm-1);

	  unsigned I = 0;
          eddyDiffModel->timeStep(dt, windModel->heff());
          if ( (releaseTime > 0) && ( (releaseTime + dt[I] - floor((releaseTime + dt[I])/mfreq)*mfreq) == 0 ) && (tm + 1 <= nor) && (pi % ppr == 0) )
          {
	    newPuff = true;
            releaseTime += dt[I];
            tm++;
            eddyDiffModel->updateParams(tm-1); windModel->updateParams(tm-1); 
	    eddyDiffModel->timeStep(dt, windModel->heff()); 
	    if ( isMaster ) { cout << "\ntime: " << releaseTime << " s" << endl; }
          }
          else if (pi % ppr == 0)
          {
	     newPuff = true;
             releaseTime += dt[I];
	     if ( isMaster ) { cout << "\ntime: " << releaseTime << " s" << endl; }
          }

	  Nt = ceil((tm*mfreq-releaseTime)/eddyDiffModel->minTimeStep());
	  if ( newPuff )
	  {
	    windModel->updateWindDirNoise(xsi[kk]); kk++;
	  }
          tm *= mfreq;
	  tp = releaseTime; //+dt[0];
	  tp_z = releaseTime;
          zt[0] = windModel->heff();

	  windModel->U(Uvec0,zt[0]);
   	  pathIntU[0] = Uvec0[0]*dt[0]; //dt[1]
   	  pathIntU[1] = Uvec0[1]*dt[0]; //dt[1]
   	  eddyDiffModel->tauL(tL,zt[0]);
   	  eddyDiffModel->eddyDiff(beta,zt[0]);
   	  eddyDiffModel->autoCorr(corr,zt[0]); //,dt[1]);
   	  tL[0] *= -1; tL[1] *= -1; tL[2] *= -1;

	  unstableTrue = eddyDiffModel->boolUnstable();
          eddyDiffModel->stdevWind(sigma_u,zt[0]);
          u0[0] = sigma_u[0]*xsi[kk]; 
          u0[1] = sigma_u[1]*xsi[kk+1]; 
          u0[2] = sigma_u[2]*xsi[kk+2]; zdot[0] = u0[2] + xdot_eject[2]; 
	  eddyDiffModel->tauL(tL0, 0); 
	  eddyDiffModel->stdevWind(sigma_u0, 0);
	  eddyDiffModel->eddyDiff(beta0,0);

          for (unsigned i=0; i < 8; i++)
          {
            var_xy_dummy[i]=0;
          }
	  var_xy_dummy[8] = 1.0;

   	    tLbetaSqr=beta[0]*beta[0]*tL[0], corrSqr = corr[0]*corr[0];
   	    var_xy[0] = (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[0]+dt[0])*tLbetaSqr*tL[0];
   	    var_xy[0] += 2*(corr[0]-1)*tL[0]*var_xy_dummy[1]+(corr[0]-1)*(corr[0]-1)*tL[0]*tL[0]*var_xy_dummy[0];
   	    integralCorr[0] = (corr[0]-1)*tL[0]*var_xy_dummy[8];

   	    var_xy_dummy[0] = corrSqr*(var_xy_dummy[0]+0.5*tLbetaSqr)-0.5*tLbetaSqr; //Itilde
   	    var_xy_dummy[2] = corrSqr*(var_xy_dummy[2]+var_xy_dummy[3]); //Ihat3
   	    var_xy_dummy[3] = -0.5*tLbetaSqr*(1-corrSqr);
   	    var_xy_dummy[1] = corr[0]*var_xy_dummy[1]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[0]*0.5+(1-1/corr[0])*var_xy_dummy[2]*tL[0]; //Ihat1

   	    tLbetaSqr=beta[1]*beta[1]*tL[1];
   	    var_xy[1] = (0.5*((corr[0]-2)*(corr[0]-2)-1)*tL[1]+dt[0])*tLbetaSqr*tL[1];
   	    var_xy[1] += 2*(corr[0]-1)*tL[1]*var_xy_dummy[5]+(corr[0]-1)*(corr[0]-1)*tL[1]*tL[1]*var_xy_dummy[4];
   	    integralCorr[1] = (corr[0]-1)*tL[1]*var_xy_dummy[8];

   	    var_xy_dummy[4] = corrSqr*(var_xy_dummy[4]+0.5*tLbetaSqr)-0.5*tLbetaSqr;
   	    var_xy_dummy[6] = corrSqr*(var_xy_dummy[6]+var_xy_dummy[7]);
   	    var_xy_dummy[7] = -0.5*tLbetaSqr*(1-corrSqr);
   	    var_xy_dummy[5] = corr[0]*var_xy_dummy[5]+(1-corr[0])*(1-corr[0])*tLbetaSqr*tL[1]*0.5+(1-1/corr[0])*var_xy_dummy[6]*tL[1];
   	    var_xy_dummy[8] *= corr[0];

	    tp += dt[0];

	}
	else 
	{
	  tm += mfreq;
	  windModel->updateParams(tm/mfreq-1); eddyDiffModel->updateParams(tm/mfreq-1);
	  if ( windModel->windDirFluc(tm/mfreq-1) == 0 )
	  {
	    windModel->updateWindDirNoise(xsi[kk]); kk++;
	    windModel->updateWindDirFluc();
	  }
	}

        delete[] xsi;
        unsigned long* nocp = NULL;
        double *bandWidths = NULL, *distrSigmaVec = NULL;
        if ( newPuff && noProcs > 1 )
        {
           nocp = new unsigned long[nor]; //number of ocntributing particles
           for (unsigned i = 0; i < nor; i++)
                nocp[i] = meanDrift[4*i+3];

           double* meanDriftContainer = NULL;
           if ( isMaster )
             meanDriftContainer = new double[4*nor*noProcs];

           MPI_Gather(meanDrift, 4*nor, MPI_DOUBLE, meanDriftContainer, 4*nor, MPI_DOUBLE, master, MPI_COMM_WORLD);
           for ( unsigned i = 0; isMaster && i < 4*nor*noProcs; i++ )
           {
              if ( i/(4*nor) == procid )
                continue;

              unsigned ii = i % (4*nor);
              meanDrift[ii] += meanDriftContainer[i];
           }
           if ( isMaster )
             delete[] meanDriftContainer;

           MPI_Bcast(meanDrift, 4*nor, MPI_DOUBLE, master, MPI_COMM_WORLD);
           bandWidths = new double[2*nor];
        }
        else if ( newPuff )
        {
           nocp = new unsigned long[nor]; //number of ocntributing particles
           for (unsigned i = 0; i < nor; i++)
                nocp[i] = meanDrift[4*i+3];

           bandWidths = new double[2*nor];
        }

        for (unsigned long i = 0; newPuff && i < nor*ppr; i++)
        {
           unsigned long r = i/ppr, ii = i % ppr, j = (r*ppr + ii)*3;
           if ( ii >= nocp[r] )
                continue;
           else if ( ii == 0 )
           {
             distrSigmaVec = new double[ppr+4]();
             distrSigmaVec[0] = nocp[r];
             meanDrift[4*r+2] /= meanDrift[4*r+3]; 
	   }

           if ( ii < nocp[r] )
           {
                distrSigmaVec[2] += (positions[j+2]-meanDrift[4*r+2])*(positions[j+2]-meanDrift[4*r+2]);
                distrSigmaVec[ii + 4] = positions[j+2];
           }

           if ( ii == nocp[r] - 1 && isMaster && noProcs > 1 )
           {
             double* h[2] = { &bandWidths[r], &bandWidths[nor+r]};
             double* distrSigmaContainer = new double[(4+ppr)*noProcs]();
             MPI_Gather(distrSigmaVec, 4+ppr, MPI_DOUBLE, distrSigmaContainer, 4+ppr, MPI_DOUBLE, master, MPI_COMM_WORLD);
             densityKernel::bandWidth(h, distrSigmaContainer, noProcs, 4+ppr, prodKernel);
             double h_cpy[2] = {*h[0],*h[1]};
             MPI_Bcast(h_cpy, 2, MPI_DOUBLE, master, MPI_COMM_WORLD);
             meanDrift[4*r] = 0; meanDrift[4*r+1] = 0; meanDrift[4*r+2] = 0; meanDrift[4*r+3] = 0;
             delete[] distrSigmaVec; delete[] distrSigmaContainer;
	   }
           else if ( ii == nocp[r] - 1 && noProcs > 1 )
           {
             double h[2]={};
             MPI_Gather(distrSigmaVec, 4+ppr, MPI_DOUBLE, NULL, 4+ppr, MPI_DOUBLE, master, MPI_COMM_WORLD);
             MPI_Bcast(h, 2, MPI_DOUBLE, master, MPI_COMM_WORLD);
             bandWidths[r] = h[0]; bandWidths[nor+r] = h[1];
             meanDrift[4*r] = 0; meanDrift[4*r+1] = 0; meanDrift[4*r+2] = 0; meanDrift[4*r+3] = 0;
             delete[] distrSigmaVec;
           }
           else if ( ii == nocp[r] - 1 )
           {
             double* h[2] = { &bandWidths[r], &bandWidths[nor+r]};
             densityKernel::bandWidth(h, distrSigmaVec, 1, 4+ppr, prodKernel);
             meanDrift[4*r] = 0; meanDrift[4*r+1] = 0; meanDrift[4*r+2] = 0; meanDrift[4*r+3] = 0;
             delete[] distrSigmaVec;
           }
        }

        {
         unsigned long i = 0, r = 0, ii = 0, jj = 0, k = 0;
         double h[2], Cda, sine, cosine; 
         if ( newPuff )
         {
           h[1] = bandWidths[nor]; Cda = 4.0/3.0; sine = sin(windModel->wdir(0)); cosine = cos(windModel->wdir(0)); 
         }
         while ( newPuff && i < nor*noStations )
         {
            if ( nocp[r] == 0 && i + noStations < nor*noStations )
            {
                i += noStations; r++; jj += noStations;
                h[1] = bandWidths[nor+r]; Cda = 4.0/3.0; sine = sin(windModel->wdir(r)); cosine = cos(windModel->wdir(r)); 
                continue;
            }
            else if ( nocp[r] == 0 )
            {
                i += noStations; continue;
            }

           unsigned long j = 3*(r*ppr+k);

           double dCxy = -0.5*(stations[ii]+positions[j])*(stations[ii]+positions[j])*(cosine*cosine/pathVar[(j*2)/3] + sine*sine/pathVar[(j*2)/3+1]);
           dCxy -= 0.5*(stations[ii+1]+positions[j+1])*(stations[ii+1]+positions[j+1])*(sine*sine/pathVar[(j*2)/3] + cosine*cosine/pathVar[(j*2)/3+1]);
	   dCxy -= (stations[ii]+positions[j])*(stations[ii+1]+positions[j+1])*sine*cosine*(1/pathVar[(j*2)/3] - 1/pathVar[(j*2)/3+1]);
           dCxy = exp(dCxy)/(2*Pi*sqrt(pathVar[(j*2)/3]*pathVar[(j*2)/3+1]));
           double dist2 = (stations[ii+2] - positions[j+2])*(stations[ii+2] - positions[j+2])/(h[1]*h[1]);

            if ( dist2 < 1 )
            {
              double z = 1-dist2;
              C[jj] += dQ*dCxy*z/(h[1]*Cda);
            }

            dist2 = (stations[ii+2] + positions[j+2])*(stations[ii+2] + positions[j+2])/(h[1]*h[1]);

            if ( dist2 < 1 )
            {
              double z = 1-dist2; 
              C[jj] += dQ*dCxy*z/(h[1]*Cda);
            }

            k++;
            if ( k == nocp[r] && i+1 < nor*noStations )
            {
                i++; k = 0; r = i/noStations; ii = i % noStations; jj = r*noStations+ii; ii = 3*i;
                h[1] = bandWidths[nor+r]; Cda = 4.0/3.0; sine = sin(windModel->wdir(r)); cosine = cos(windModel->wdir(r)); // densityKernel::integralKernel( prodKernel ); 
            }
            else if ( k == nocp[r] )
               i++;
         }
        }

	double* C_slaves = NULL;
        if ( isMaster && newPuff && noProcs > 1 ) 
	   C_slaves = new double[3*nor*noStations*noProcs];

	if ( newPuff && noProcs > 1 )
	{
	  MPI_Gather(C, 3*nor*noStations, MPI_DOUBLE, C_slaves, 3*nor*noStations, MPI_DOUBLE, master, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);
	}

        if ( newPuff )
        {
          delete[] nocp; delete[] bandWidths;
          fill(meanDrift,meanDrift+4*nor,'\0');
	  fill(positions,positions+3*ppr*nor,'\0');
	  fill(pathVar,pathVar+2*ppr*nor,'\0');
	  windModel->windDirFlucReset();
          eddyDiffModel->timeStep(dt, zt[0]);
	  windModel->updateWindDirFluc();
          dQ = Q*dt[1]/(ppr*noProcs); //*dt[0]
          MPI_Barrier(MPI_COMM_WORLD);
        }

	if (warmingUp && pi >= ppr && tm > 0)
	{
	   warmingUp = false; 
	}

   } while( ceil(releaseTime*1000)/1000 < endTime );

   vslDeleteStream( &stream);

   if ( isMaster && noProcs > 1 )
   {
     double* CContainer = new double[nor*noStations*noProcs];
     MPI_Gather(C, nor*noStations, MPI_DOUBLE, CContainer, nor*noStations, MPI_DOUBLE, master, MPI_COMM_WORLD);
     unsigned long i = 0;
     while ( i < nor*noStations*noProcs )
     {
        if ( i/(nor*noStations) == master )
        {
          i += nor*noStations; continue;
        }

        unsigned long ii = i % (nor*noStations);
        C[ii] += CContainer[i];
        i++;
     }
     delete[] CContainer;
   }
   else if ( noProcs > 1 )
   {
     MPI_Gather(C, nor*noStations, MPI_DOUBLE, NULL, nor*noStations, MPI_DOUBLE, master, MPI_COMM_WORLD);
   }

  MPI_Finalize();

   if ( isMaster )
   {
     cout << "\nWriting output" << endl;
     ofstream ofile("output/CONC_" + windModel->meteofile());
     for (unsigned i = 0; i < noStations*nor; i++)
     {
        ofile << i/noStations << setw(15) << i % noStations << setw(15) << C[i] << endl; 
     }
     ofile.close();
   }

   delete windModel, delete eddyDiffModel, delete[] C, delete[] stations, delete[] meanDrift, delete[] positions, delete[] pathVar; 

   return 0;
 }

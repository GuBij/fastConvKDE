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
   bool prodKernel = false;
   bool refSolution = false;
   unsigned long noStations;
   unsigned seedNr = 1;
   double *stations, *cellVol;  
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
	else if ( string(argv[i]) == "-prodKernel" )
	{
	    prodKernel = true;
	}
        else if ( string(argv[i]) == "-seed" )
        {
            i++;
            stringstream ss(argv[i]);
            ss >> seedNr;
        }
        else if ( string(argv[i]) == "-refSolution" )
        {
            refSolution = true;
        }
	else
	{
	    cout << "\n\t ERROR: option '" << argv[i] << "' is not implemented. Valid options are: \n" << endl;
	    cout << "\t-Langevin\n\t-stationary\n\t-prodKernel\n\t-seed\n\t-refSolution" << endl;
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

     cellVol = new double[noStations]();
     ifile.open("cellVolumes.txt", ifstream::in);
     if (ifile.is_open())
     {
        for (unsigned long j=0; j<noStations; j++)
        {
           ifile >> cellVol[j];
        }
     }
     else
     {
        cerr << "\nWARNING: file 'cellVolumes.txt' not found. All Cells are assumed to have unit volume.\n" << endl;

        for (unsigned long j=0; j<noStations; j++)
        {
           cellVol[j] = 1.0;
        }
     }
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

   int mfreq = dict.getfl(2);
   double xMax = dict.getfl(4), yMax = dict.getfl(5), zMax = dict.getfl(6);
   const unsigned long ppr = ceil(dict.getfl(3)/noProcs); //dict.getfl(8));

   if ( isMaster ) 
   {
	if ( LangevinModel ) { cout << "\nLangevin model selected." << endl; }
	if ( stationary ) { cout << "\nStationary release assumed." << endl; }
	if ( prodKernel ) { cout << "\nProduct kernel will be used for concentration estimation." << endl; }
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

   const double Q = dict.getfl(7)*1000; // /3600;
   const unsigned Nt = ceil(double(mfreq)/eddyDiffModel->minTimeStep());
   unsigned long pi = 0;
   int warmingUpPeriod = 0;
   bool unstableTrue = eddyDiffModel->boolUnstable();
   long tm = mfreq;
   double releaseTime = 0;
   double xt[3] = {}; xt[2] = windModel->heff();
   const double xdot[3] = {};
   double ut[3] = {};
   double tL[3]; 
   double* meanDrift = new double[4*nor]();
   double* positions = new double[3*ppr*nor]();

   double MISE = 0;
   double* Cref = new double[nor*noStations]();
/*
   if ( !refSolution )
   {
     ifstream ifile("output/puffReleases_" + windModel->meteofile() + "/CONC_PUFF0.txt", ifstream::in);
     double skip;
     for (unsigned long i=0; i<noStations*nor; i++)
     {
       ifile >> skip >> skip >> Cref[i];
     }
     ifile.close();
   }
*/
   unsigned seed_(noProcs*seedNr + procid); // *(nor+3)*mfreq*ppr); //+3 to estimate the warming up period
   VSLStreamStatePtr stream;
   vslNewStream( &stream, VSL_BRNG_MT19937, seed_);

   {
//     double Uvec0[3]; 
//     windModel->U(Uvec0,xt[2]);
     const int status = vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, 3, ut, 0.0, 1.0);
     if (status != VSL_ERROR_OK)
     {
        cout << "\n\t ERROR: sampling by function 'vdRngGaussian' failed; status flag: " << status << "\n" << endl;
        exit(EXIT_FAILURE);
     }
     double sigma_u[3];
     eddyDiffModel->stdevWind(sigma_u,xt[2]);
     ut[0] *= sigma_u[0]; // ut[0] -= Uvec0[0]; // xdot[0];
     ut[1] *= sigma_u[1]; // ut[1] -= Uvec0[1]; // xdot[1];
     ut[2] *= sigma_u[2]; // ut[2] -= Uvec0[2]; // xdot[2];
//     cout << "\n U0_x: " << Uvec0[0] << ", dt: " << dt[0] << ", sigma_u: " << sigma_u[0] << endl;
   }

   double tp = 0;
   double dQ = Q; // *dt[0]; 
   unsigned puffId = 0;
   do {
	double* xsi = new double[3*(Nt+1)]();
	const int status = vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, 3*(Nt+1), xsi, 0.0, 1.0);
	if (status != VSL_ERROR_OK)
	{
	  cout << "\n\t ERROR: sampling by function 'vdRngGaussian' failed; status flag: " << status << "\n" << endl;
	  exit(EXIT_FAILURE);
	}

	unsigned kk = 0;
	double Uvec[3], beta[3], dt_hor = 0.0, dt_hor_max = 0.0, tL0[3] = {}, sigma_u0[3] = {}, beta0[3] = {}; 
        windModel->U(Uvec,xt[2]);
        eddyDiffModel->tauL(tL,xt[2]);
        eddyDiffModel->eddyDiff(beta,xt[2]);
        tL[0] = -1/tL[0]; tL[1] = -1/tL[1]; tL[2] = -1/tL[2];
	eddyDiffModel->tauL(tL0, 0); eddyDiffModel->stdevWind(sigma_u0, 0); eddyDiffModel->eddyDiff(beta0,0);
	while ( ceil(tp*1000)/1000 < tm || dt_hor < dt_hor_max )
	{
          double A0; eddyDiffModel->nonHomTurbCorr(A0,xt[2],ut[2]);
	  double zz[2] = {xt[2],Uvec[2]+ut[2]};
          xt[2] += (Uvec[2]+ut[2])*dt[1];
          if ( !unstableTrue )
          {
             double A0_; eddyDiffModel->nonHomTurbCorr(A0_,xt[2],ut[2]); A0 *= 0.5; A0 += 0.5*A0_;
          }
          ut[2] = (1+dt[1]*tL[2])*ut[2] + A0*dt[1] + beta[2]*sqrt(dt[1])*xsi[kk];
	  kk++;

          if ( ceil(tp*1000)/1000 < tm && ((( dt_hor <= dt_hor_max ) && ( dt_hor+dt[1] >= dt_hor_max )) || dt_hor_max < dt_hor )) // dt_hor >= dt_hor_max
          {
	    double zzz;
            if ( (abs(dt_hor - dt_hor_max) < dt_hor - dt_hor_max + dt[1]) || dt_hor_max < dt_hor  )
                  zzz = zz[0] + zz[1]*(dt_hor_max - dt_hor);
            else
                  zzz = xt[2] + zz[1]*(dt_hor_max - dt_hor - dt[1]);

           if ( zzz < 0 && unstableTrue )
            {
               const double t0 = -zz[0]/zz[1],  varGround = sigma_u0[2]*sigma_u0[2]*(1-exp(-2*(dt[1]-t0)/tL0[2]));
               double utr = -2*varGround*log(1-exp(-0.5*zz[1]*zz[1]/varGround)); utr = sqrt(utr); // ut reflected
               zzz = (dt[1]-t0)*utr;
            }
            else if ( zzz < 0 )
               zzz *= -1;

	    windModel->U(Uvec,zzz);
	    dt_hor = dt_hor - dt_hor_max + dt[1];
	    eddyDiffModel->timeStep(dt,zzz);
	    eddyDiffModel->tauL(tL,zzz);
	    eddyDiffModel->eddyDiff(beta,zzz);
	    tL[0] = -1/tL[0]; tL[1] = -1/tL[1];

            xt[0] += (Uvec[0]+ut[0])*dt[0]; 
            xt[1] += (Uvec[1]+ut[1])*dt[0]; 
	    ut[0] = (1+dt[0]*tL[0])*ut[0] + beta[0]*sqrt(dt[0])*xsi[kk]; 
	    ut[1] = (1+dt[0]*tL[1])*ut[1] + beta[1]*sqrt(dt[0])*xsi[kk+1]; 
	    kk += 2;
	    tp += dt[0];
	    dt_hor_max = dt[0];
	  }
	  else
		dt_hor += dt[1];

            const double uti = zz[1];
            if ( unstableTrue && xt[2] < 0 )
            {
                const double t0 = -zz[0]/zz[1], varGround = sigma_u0[2]*sigma_u0[2]*(1-exp(-2*(dt[1]-t0)/tL0[2]));
                zz[1] = -2*varGround*log(1-exp(-0.5*zz[1]*zz[1]/varGround)); zz[1] = sqrt(zz[1]);
                xt[2] = (dt[1]-t0)*zz[1];
                ut[2] = (1-(dt[1]-t0)/tL0[2])*zz[1] + beta0[2]*sqrt(dt[1]-t0)*xsi[kk];
            }
            else if ( xt[2] < 0 )
            {
              xt[2] *= -1; ut[2] *= -1; // dt_prev *= -1;
            }

            if ( unstableTrue && tp >= tm && dt_hor >= dt_hor_max )
            {
                const double zt0 = xt[2];
                xt[2] += zz[1]*(tm - tp + dt_hor_max - dt_hor);
                if ( xt[2] < 0 )
                {
                  const double t0 = tp - dt_hor_max + dt_hor - zt0/zz[1];
                  xt[0] = (tm-t0)*uti;
                }
                tp = tm; dt_hor = dt_hor_max;
            }
            else if ( tp >= tm && dt_hor >= dt_hor_max )
            {
		xt[2] += zz[1]*(tm - tp + dt_hor_max - dt_hor); xt[2] = abs(xt[2]);
                tp = tm; dt_hor = dt_hor_max;
            }
            else
            {
                eddyDiffModel->tauL(tL,xt[2]);
                eddyDiffModel->eddyDiff(beta,xt[2]);
                tL[2] = -1/tL[2];
            }
            eddyDiffModel->timeStep(dt,xt[2]);
	}
//	cout << "\ntp: " << tp << ", dt/tauL_w: " << -dt[1]*tL[2] << ", kk: " << kk << ", tp: " << tp << ", tp_z: " << tp - dt_hor_max + dt_hor << endl;

	if ( !warmingUp )
	{
//	   cout << "x: " << xt[0] << ", y: " << xt[1] << ", z: " << xt[2] << "\n" << endl;
	   unsigned long i = 4*(tm/mfreq-1), ii = ((tm/mfreq-1)*ppr + (pi % ppr))*3;
	   meanDrift[i] += xt[0]; meanDrift[i+1] += xt[1]; meanDrift[i+2] += xt[2]; meanDrift[i+3] += 1;
	   positions[ii] = xt[0]; positions[ii+1] = xt[1]; positions[ii+2] = xt[2];
//	   pathIntU[1] += pathIntU[0];
	}

	const bool death = (tm/mfreq == nor) || (abs(xt[0]) > xMax) || (abs(xt[1]) > yMax ) || (abs(xt[2]) > zMax);
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
	  for (unsigned i = 0; i < 3; i++)
	  {
		xt[i] = 0; ut[i] = 0;
	  }
	  xt[2] = windModel->heff();

	  {
	    double Uvec0[3]; 
	    windModel->U(Uvec0,xt[2]);
            unstableTrue = eddyDiffModel->boolUnstable();
	    double sigma_u[3];
	    eddyDiffModel->stdevWind(sigma_u,xt[2]);
	    ut[0] = - Uvec0[0] + sigma_u[0]*xsi[kk]; // xdot[0] + sigma_u[0]*xsi[kk]; // - Uvec0[0];
	    ut[1] = - Uvec0[1] + sigma_u[1]*xsi[kk+1]; // xdot[1] + sigma_u[1]*xsi[kk+1]; // - Uvec0[1];
	    ut[2] = - Uvec0[2] + sigma_u[2]*xsi[kk+2]; // xdot[2] + sigma_u[2]*xsi[kk+2]; //  - Uvec0[2];
            eddyDiffModel->tauL(tL0, 0);
            eddyDiffModel->stdevWind(sigma_u0, 0);
            eddyDiffModel->eddyDiff(beta0,0);
	  }
	}
	else if (death)
	{
          pi++;
          tm = releaseTime/mfreq+1;
	  windModel->updateParams(tm-1); eddyDiffModel->updateParams(tm-1);

          eddyDiffModel->timeStep(dt, windModel->heff());
          if ( (releaseTime > 0) && ( (releaseTime + dt[0] - floor((releaseTime + dt[0])/mfreq)*mfreq) == 0 ) && (tm + 1 <= nor) && (pi % ppr == 0) )
          {
	    newPuff = true;
            releaseTime += dt[0];
            tm++;
            eddyDiffModel->updateParams(tm-1); windModel->updateParams(tm-1); //windParamPtr);
	    if ( isMaster ) { cout << "\ntime: " << releaseTime << " s" << endl; }
          }
          else if (pi % ppr == 0)
          {
	     newPuff = true;
             releaseTime += dt[0];
	     if ( isMaster ) { cout << "\ntime: " << releaseTime << " s" << endl; }
          }
          tm *= mfreq;
	  tp = releaseTime;
          for (unsigned i = 0; i < 3; i++)
          {
                xt[i] = 0; ut[i] = 0; 
          }
          xt[2] = windModel->heff();

          {
//            double Uvec0[3]; 
//            windModel->U(Uvec0,xt[2]);
	    unstableTrue = eddyDiffModel->boolUnstable();
            double sigma_u[3];
            eddyDiffModel->stdevWind(sigma_u,xt[2]);
            ut[0] = sigma_u[0]*xsi[kk]; // xdot[0] + sigma_u[0]*xsi[kk]; // - Uvec0[0];
            ut[1] = sigma_u[1]*xsi[kk+1]; // xdot[1] + sigma_u[1]*xsi[kk+1]; // - Uvec0[1];
            ut[2] = sigma_u[2]*xsi[kk+2]; // xdot[2] + sigma_u[2]*xsi[kk+2]; //  - Uvec0[2];
            eddyDiffModel->tauL(tL0, 0);
            eddyDiffModel->stdevWind(sigma_u0, 0);
            eddyDiffModel->eddyDiff(beta0,0);
//	    cout << "\n U0_x: " << Uvec0[0] << ", dt: " << dt[0] << ", sigma_u*u: " << sigma_u[0]*xsi[kk] << endl;
          }
	}
	else 
	{
	  tm += mfreq;
	  windModel->updateParams(tm/mfreq-1); eddyDiffModel->updateParams(tm/mfreq-1);
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

//	double minMaxVal[6] = {2*xMax,-2*xMax,2*yMax,-2*yMax,0,0};
	for (unsigned long i = 0; newPuff && i < nor*ppr; i++)
	{
	   unsigned long r = i/ppr, ii = i % ppr, j = (r*ppr + ii)*3;
	   if ( ii >= nocp[r] )
		continue;
	   else if ( ii == 0 )
	   {
	     distrSigmaVec = new double[ppr+4]();
	     distrSigmaVec[0] = nocp[r];
	     meanDrift[4*r] /= meanDrift[4*r+3]; meanDrift[4*r+1] /= meanDrift[4*r+3]; meanDrift[4*r+2] /= meanDrift[4*r+3];
//	     cout << "mean z: " << pathIntU[1]/meanDrift[4*r+3] << "\n" << endl;
	   }

	   if ( ii < nocp[r] )
	   {
		distrSigmaVec[1] += (positions[j]-meanDrift[4*r])*(positions[j]-meanDrift[4*r]);
		distrSigmaVec[2] += (positions[j+2]-meanDrift[4*r+2])*(positions[j+2]-meanDrift[4*r+2]);
		distrSigmaVec[3] += (positions[j+1]-meanDrift[4*r+1])*(positions[j+1]-meanDrift[4*r+1]);
		distrSigmaVec[ii + 4] = positions[j+2];
	   }
/*
	   else if ( ii < nocp[r] )
	   {
	     distrSigmaVec[1] += (positions[j]-meanDrift[4*r])*(positions[j]-meanDrift[4*r])+(positions[j+1]-meanDrift[4*r+1])*(positions[j+1]-meanDrift[4*r+1])+(positions[j+2]-meanDrift[4*r+2])*(positions[j+2]-meanDrift[4*r+2]);
	     distrSigmaVec[ii+4] = sqrt((positions[j]-meanDrift[4*r])*(positions[j]-meanDrift[4*r])+(positions[j+1]-meanDrift[4*r+1])*(positions[j+1]-meanDrift[4*r+1])+(positions[j+2]-meanDrift[4*r+2])*(positions[j+2]-meanDrift[4*r+2]));
//sqrt(positions[j]*positions[j] + positions[j+1]*positions[j+1] + positions[j+2]*positions[j+2]);
	   }

	     for ( unsigned m = 0; m < 3; m++ )
	     {
	       if ( positions[j+m] < minMaxVal[2*m] )
		  minMaxVal[2*m] = positions[j+m];

	       if ( positions[j+m] > minMaxVal[2*m+1] )
		  minMaxVal[2*m+1] = positions[j+m];
	     }
//	     sd = positions[j+2]-meanDrift[4*r+2]; 
//	     distrSigmaVec[2] += sd*sd;
//	     distrSigmaVec[ppr+ii+3] = sd;
*/
	   if ( ii == nocp[r] - 1 && isMaster && noProcs > 1 )
	   {
	     double* h[2] = { &bandWidths[r], &bandWidths[nor+r]};
	     double* distrSigmaContainer = new double[(4+ppr)*noProcs]();
	     MPI_Gather(distrSigmaVec, 4+ppr, MPI_DOUBLE, distrSigmaContainer, 4+ppr, MPI_DOUBLE, master, MPI_COMM_WORLD);
	     densityKernel::bandWidth(h, distrSigmaContainer, noProcs, 4+ppr, prodKernel);
	     double h_cpy[2] = {*h[0],*h[1]};
	     MPI_Bcast(h_cpy, 2, MPI_DOUBLE, master, MPI_COMM_WORLD);
//	     bandWidths[r] = h; //bandWidths[nor+r] = h[1];
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
//	     densityKernel::bandWidthHist(h, distrSigmaVec, 1, 1+2*ppr, minMaxVal);
//	     bandWidths[3*r] = h[0]; bandWidths[3*r+1] = h[1]; bandWidths[3*r+2] = h[2];
	     densityKernel::bandWidth(h, distrSigmaVec, 1, 4+ppr, prodKernel);
//	     cout << "bandWidth: " << h[0] << ", " << h[1] << endl;
	     meanDrift[4*r] = 0; meanDrift[4*r+1] = 0; meanDrift[4*r+2] = 0; meanDrift[4*r+3] = 0;
//	     minMaxVal[0] = 2*xMax; minMaxVal[1] = -2*xMax; minMaxVal[2] = 2*yMax; minMaxVal[3] = -2*yMax; minMaxVal[3] = 0; minMaxVal[3] = 0;
	     delete[] distrSigmaVec;
	   }
	}

	{
	 unsigned long i = 0, r = 0, ii = 0, jj = 0, k = 0;
	 double h[2], h2, Cda;
	 if ( newPuff )
	 {
	   h[0] = bandWidths[0]; h[1] = bandWidths[nor]; h2 = h[0]*h[0]; Cda = densityKernel::integralKernel( prodKernel ); //stations[2],h,zMax,prodKernel);
//	   h[0] = bandWidths[0]; h[1] = bandWidths[1]; h[2] = bandWidths[2]; 
	 }
	 while ( newPuff && i < nor*noStations ) 
	 {
	    if ( nocp[r] == 0 && i + noStations < nor*noStations )
	    {
		i += noStations; r++; jj += noStations; 
//		h[0] = bandWidths[3*r]; h[1] = bandWidths[3*r+1]; h[2] = bandWidths[3*r+2]; 
		h[0] = bandWidths[r]; h[1] = bandWidths[nor+r]; h2 = h[0]*h[0]; Cda = densityKernel::integralKernel( prodKernel ); //stations[2],h,zMax,prodKernel);
		continue;
	    }
	    else if ( nocp[r] == 0 )
	    {
		i += noStations; continue;
	    }

	    unsigned long j = 3*(r*ppr+k);

	    double dist2 = (stations[ii] - positions[j])*(stations[ii] - positions[j])/h2;
	    dist2 += (stations[ii+1] - positions[j+1])*(stations[ii+1] - positions[j+1])/h2;
	    dist2 += (stations[ii+2] - positions[j+2])*(stations[ii+2] - positions[j+2])/h2;

            double dC = 0;
	    if ( prodKernel )
	    {
		double distVert2 = (stations[ii+2] - positions[j+2])*(stations[ii+2] - positions[j+2]);
		dist2 -= distVert2/h2;
		distVert2 /= h[1]*h[1];
		if ( dist2 < 1 && distVert2 < 1 )
		{
		  double x = 1-dist2, z = 1-distVert2;
		  dC += dQ*x*z/(h2*h[1]*Cda); 
//		  double z = (1-dist2)*(1-4*dist2);
//                  dC += -6*z/(h[1]*h[1]*h[1]*96.0/105.0);
//		  C[jj] += dQ*x*x*x*z*z*z/(h2*h[1]*Cda);
		}
	    }
	    else if ( dist2 < 1 ) //distx < h[0] && disty < h[1] && abs(stations[ii+2] - (floor(positions[j+2]/h[2])+0.5)*h[2]) < h[2] ) // && distVert2 < 1 )
	    {
		double x = 1-dist2; // z = 1-distVert2;
		dC += dQ*x/(h2*h[0]*Cda);
	    }

            dist2 = (stations[ii] - positions[j])*(stations[ii] - positions[j])/h2;
            dist2 += (stations[ii+1] - positions[j+1])*(stations[ii+1] - positions[j+1])/h2;
            dist2 += (stations[ii+2] + positions[j+2])*(stations[ii+2] + positions[j+2])/h2;

	    if ( prodKernel )
            {
		double distVert2 = (stations[ii+2] + positions[j+2])*(stations[ii+2] + positions[j+2]);
                dist2 -= distVert2/h2;
                distVert2 /= h[1]*h[1];
                if ( dist2 < 1 && distVert2 < 1 )
                {
                  double x = 1-dist2, z = 1-distVert2;
		  dC += dQ*x*z/(h2*h[1]*Cda); 
//		  double z = (1-dist2)*(1-4*dist2);
//                  dC += -6*z/(h[1]*h[1]*h[1]*96.0/105.0);
//                  C[jj] += dQ*x*x*x*z*z*z/(h2*h[1]*Cda);
                }
	    }
            else if ( dist2 < 1 ) // && distVert2 < 1 )
            {
                double x = 1-dist2; // z = 1-distVert2;
                dC += dQ*x/(h2*h[0]*Cda);
            }

	    Cref[jj] += dC/(ppr*noProcs);
/*
            if ( refSolution )
	    {
                Cref[jj] += dC/(ppr*noProcs);
	    }
            else
	    {
                MISE += (dC-Cref[jj])*(dC-Cref[jj])/(ppr*noProcs*dQ*dQ);
	    }
*/
	    k++;
	    if ( k == nocp[r] && i+1 < nor*noStations )
            {
                i++; k = 0; r = i/noStations; ii = i % noStations; jj = r*noStations+ii; ii = 3*i;
//		h[0] = bandWidths[3*r];  h[1] = bandWidths[3*r+1]; h[2] = bandWidths[3*r+2];
                h[0] = bandWidths[r]; h[1] = bandWidths[nor+r]; h2 = h[0]*h[0]; Cda = densityKernel::integralKernel( prodKernel ); //stations[ii+2],h,zMax,prodKernel);
            }
	    else if ( k == nocp[r] )
	       i++;
	 }
	}

	if ( newPuff && isMaster && noProcs > 1 )
	{
	  delete[] nocp; delete[] bandWidths;
          double* CContainer = new double[nor*noStations*noProcs];
          MPI_Gather(Cref, nor*noStations, MPI_DOUBLE, CContainer, nor*noStations, MPI_DOUBLE, master, MPI_COMM_WORLD);
          unsigned long i = 0;
          while ( i < nor*noStations*noProcs )
          {
            if ( i/(nor*noStations) == master )
            {
                i += nor*noStations; continue;
            }

            unsigned long ii = i % (nor*noStations);
            Cref[ii] += CContainer[i];
            i++;
          }
          delete[] CContainer;
	}
	else if ( newPuff && noProcs > 1 )
	{
	  delete[] nocp; delete[] bandWidths;
          MPI_Gather(Cref, nor*noStations, MPI_DOUBLE, NULL, nor*noStations, MPI_DOUBLE, master, MPI_COMM_WORLD);
          fill(Cref,Cref+noStations*nor,'\0');
          eddyDiffModel->timeStep(dt, xt[2]);
          dQ = Q; // *dt[0];
          MPI_Barrier(MPI_COMM_WORLD);
	}
	else if ( newPuff )
	{
	  delete[] nocp; delete[] bandWidths;
	}

        if ( refSolution && newPuff && isMaster )
        {
          eddyDiffModel->timeStep(dt, xt[2]);
          dQ = Q; // *dt[0];

          string odir = "output/puffReleases_" + windModel->meteofile();
          if ( puffId == 0 && mkdir(odir.c_str(),0777) == -1 )
          {
             cout << "\n\t ERROR: failed to create directory " << odir << endl;
             exit(EXIT_FAILURE);
          }

          ofstream ofile(odir + "/CONC_PUFF" + to_string(puffId) + ".txt");
          for (unsigned long i = 0; i < noStations*nor; i++)
          {
            ofile << i/noStations << setw(15) << i % noStations << setw(15) << Cref[i] << endl;
	    Cref[i] = 0;
          }
          ofile.close();
          puffId++;
	  MPI_Barrier(MPI_COMM_WORLD);
        }
        else if ( newPuff && isMaster )
        {
          eddyDiffModel->timeStep(dt, xt[2]);
          dQ = Q; // *dt[0]; 

          ifstream ifile("output/puffReleases_" + windModel->meteofile() + "/CONC_PUFF" + to_string(puffId) + ".txt", ifstream::in);
	  if ( !ifile.is_open() )
	  {
	    cerr << "\nfile 'output/puffReleases_" + windModel->meteofile() + "/CONC_PUFF" + to_string(puffId) + ".txt' not found.\n" << endl;
	    exit(EXIT_FAILURE);
	  }
          double skip, C;
          for (unsigned long i=0; i<noStations*nor; i++)
          {
	    unsigned long j = i % noStations;
            ifile >> skip >> skip >> C; // Cref[i];
	    MISE += (Cref[i]-C)*(Cref[i]-C)*cellVol[j];
	    Cref[i] = 0;
          }
          ifile.close();

	  puffId++;
          MPI_Barrier(MPI_COMM_WORLD);
        }

	if (warmingUp && pi >= ppr && tm > 0)
	{
	   warmingUp = false; 
	}

   } while( ceil(releaseTime*1000)/1000 < endTime );

   vslDeleteStream( &stream);
/*
   if ( !refSolution && isMaster && noProcs > 1 )
   {
     double* MISEContainer = new double[noProcs];
     MPI_Gather(&MISE, 1, MPI_DOUBLE, MISEContainer, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
     unsigned long i = 0;
     while ( i < noProcs )
     {
        if ( i == master )
        {
          i++; continue;
        }

        MISE += MISEContainer[i];
        i++;
     }
     delete[] MISEContainer;
   }
   else if ( !refSolution && noProcs > 1 )
   {
     MPI_Gather(&MISE, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
   }
*/
   MPI_Finalize();

   if ( !refSolution && isMaster )
   {
     ofstream ofile("output/MISE" + to_string(ppr*noProcs) + "_" + windModel->meteofile(), ofstream::app);
     ofile << seedNr << "\t" << MISE << endl;
     ofile.close();
   }

   delete windModel, delete eddyDiffModel, delete[] Cref, delete[] stations, delete[] meanDrift, delete[] positions;

   return 0;
 }

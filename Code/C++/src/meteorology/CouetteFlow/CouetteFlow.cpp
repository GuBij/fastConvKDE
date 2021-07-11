#include <iomanip>
#include <algorithm>
#include <iterator>
#include "CouetteFlow.h"

 CouetteFlow::CouetteFlow(const dictionary& dict, bool master_proc )
 :
 CouetteFlow(dict,false, master_proc)
 {}
 
 CouetteFlow::CouetteFlow(const dictionary& dict, bool LangevinModel, bool master_proc )
 :
 dictionary(dict,"CouetteFlow",{"z0","zref","zrefTemp","zrefTemp0","zrefWindVar"}),
 name_("CouetteFlow"),
 kappa(0.4),
 g(9.81),
 z0(getfl(0)),
 zref(getfl(1)),
 zref_T(getfl(2)),
 zref_T0(getfl(3)),
 zref_wvar(getfl(4)),
 pressure(101300.0),
 rho(1.225),
 R(287.0),
 cp(1004.0),
 fCori(2*7.2921e-5*sin(dict.getfl(9)*PI/180)),
 meteofile_(dict.get(0)),
 stackHeight(dict.getfl(2)),
 d(dict.getfl(3)),
 flowRate(dict.getfl(4)),
 Ts(dict.getfl(5)+273.15),
 noGridPointsZ(ceil(4*dict.getfl(6))),
 resZ(dict.getfl(6)/ceil(4*dict.getfl(6))),
 humidityType_(dict.get(8)),
 windVarType_(dict.get(10)),
 humidityUniform_(false),
 windVarIsVelocity_(false),
 useLangevin(LangevinModel)
 {
   ifstream ifile(dict.get(0));
   ifile.unsetf(ios_base::skipws);
   if (ifile.is_open())
     noRecords = count(istream_iterator<char>(ifile),istream_iterator<char>(),'\n');
   else
     cerr << "file " << dict.get(0) << "not found. Command executed from class CouetteFlow." << endl;

   ifile.setf(ios_base::skipws);
   ifile.close();

   dictionary humiDict(dict,humidityType_,{"Value"});
   if ( !( humidityType_ == "RH" || humidityType_ == "vaporPressure" ))
   {
      cout << "\n\t ERROR: entry " << humidityType_ << " is not defined. Possible values are: \n" << endl;
      cout << "\n\tRH \n\tvaporPressure" << endl;
      exit(EXIT_FAILURE);
   }
   else if ( humiDict.subDictExist() )
   {
      humidityUniform_ = true;
      humidityValue_ = humiDict.getfl(0);
   }

   if ( !( windVarType_ == "velocity" || windVarType_ == "angle" ) )
   {
      cout << "\n\t ERROR: entry " << windVarType_ << " is not defined. Possible values are: \n" << endl;
      cout << "\n\tvelocity \n\tangle" << endl;
      exit(EXIT_FAILURE);
   }
   else if ( windVarType_ == "velocity" )
   {
      windVarIsVelocity_ = true;
   }

   if ( master_proc )
      cout << "\nCouette flow model selected." << endl;
 }

 void CouetteFlow::picard(double theta0, double theta, double uref,  double& uStar, double& thetaStar, double& L)
 {
    double x0[3], x[3];
    x0[0] = uref/log(zref/z0);
    x0[1] = (theta-theta0)/log(zref_T/zref_T0);
    x0[2] = x0[0]*x0[0]*theta0/(g*x[1]);

    int noit=0, MAX_NOIT = 1000;
    double res=1.0, f_SMALL_FLOAT = 1e-04;
    while ( res > f_SMALL_FLOAT && noit < MAX_NOIT)
    {
      x[0] = uref/(log(zref/z0)-stabKernel(zref/x0[2],"M",true)+stabKernel(z0/x0[2],"M",true));
      x[1] = (theta-theta0)/(log(zref_T/zref_T0)-stabKernel(zref_T/x0[2],"H",true)+stabKernel(zref_T0/x0[2],"H",true));
      x[2] = x[0]*x[0]*theta0/(g*x[1]);
      res = (x[0]-x0[0])*(x[0]-x0[0])+(x[1]-x0[1])*(x[1]-x0[1])+(x[2]-x0[2])*(x[2]-x0[2]);
      res = sqrt(res);
      x0[0] = x[0]; x0[1] = x[1]; x0[2] = x[2];
      noit += 1;
    }

    if ( noit == MAX_NOIT )
    {
	cerr << "WARNING: maximum number of iterations is reached." << endl;
    }

    uStar = x[0]*kappa; thetaStar = x[1]*kappa; L = x[2];
 }

 void CouetteFlow::plumeRise(double T, double T0, double uStar, double L, double& heff)
 {
   double dT = T - T0;
   double Tair = T0 + dT/(zref_T-zref_T0)*(stackHeight-zref_T0) + 273.15;
   double v = flowRate/(PI*0.5*0.5*d*d);
   double F = 0.25*g*v*d*d*(Ts-Tair)/Ts;
   double s = g*(dT/(zref_T-zref_T0)+0.0098)/Tair;
   double U = uStar/kappa*(log(stackHeight/z0)-stabKernel(stackHeight/L,"M",true)+stabKernel(z0/L,"M",true));

   double xf;
   if (F < 55)
    xf = 49*pow(pow(F,10.0),1.0/16.0);
   else
    xf = 119*pow(F,0.4);

   if (L < 500 && L > 0)
   {
     if (xf > 1.84*U/pow(s,0.5))
	heff = 2.4*pow(F/(U*s),1.0/3.0);
     else
	heff = 1.6*pow(F,1.0/3.0)*pow(xf,2.0/3.0)/U;
   }
   else
   {
     int sign; (F > 0) ? (sign=1) : ((F < 0) ? (sign=-1) : (sign=0));
     heff = 1.6*sign*pow(abs(F),1.0/3.0)*pow(xf,2.0/3.0)/U;
   }

   heff = heff+stackHeight;
 }

 void CouetteFlow::calcParams(int timeIndex,double* data)
 {
   double T0 = data[0];
   double T = data[1];
   double uref=data[2];
   double windDirection = data[3];
   double elevation = data[4];

   double eps=0.622;
   double ev, ev0;
   if ( !humidityUniform_ && humidityType_ == "vaporPressure" )
   {
     ev0 = data[8]; ev = data[9];
   }
   else if ( !humidityUniform_ && humidityType_ == "RH" )
   {
     double es = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T/(240.97+T)); //611.21*exp((18.678-T/234.5)*T/(257.14+T));
     double es0 = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T0/(240.97+T0)); //(18.678-T0/234.5)*T0/(257.14+T0));
     ev0 = data[8]*es0;
     ev = data[9]*es;
   }
   else if ( humidityType_ == "RH" )
   {
     double es0 = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T0/(240.97+T0));
     double es = (1.0007+3.46*1e-06*pressure)*611.21*exp(17.502*T/(240.97+T));
     ev0 = humidityValue_*es0;
     ev = humidityValue_*es;
   }
   else
   {
     ev0 = humidityValue_; ev = ev0;
   }

   double rv0 = eps*ev0/(pressure-ev0);
   double theta0 = (T0+273.15)*(1+rv0/eps)/(1+rv0)*pow(pressure/(pressure-g*rho*zref_T0),R/cp);

   double rv = eps*ev/(pressure-ev);
   double theta = (T+273.15)*(1+rv/eps)/(1+rv)*pow(pressure/(pressure-g*rho*zref_T),R/cp);

   double uStar, thetaStar, L, heff, hABL = 0.0;
   picard(theta0, theta, uref, uStar, thetaStar, L);
   plumeRise(T,T0,uStar,L,heff);

   int width=10;
   ofstream Ufile("output/PARAMS_UCouette_" + meteofile_, ios::app);
   if (Ufile.is_open())
   {
        Ufile << left << setw(width) << timeIndex
             << left << setw(width) << uStar
	     << left << setw(width) << uref // L 
	     << left << setw(width) << heff
             << left << setw(width) << windDirection
	     << left << setw(width) << elevation << "\n";
        Ufile.close();
   }

     double sigmaAz = data[5];
     double sigmaEl = data[7];
     double sigma_v = uref*sigmaAz*PI/180;

     double sigma_u = 0.0, sigma_w = 0.0, propConst = 0.0;
     if ( windVarIsVelocity_ )
     {
	sigma_u = data[5]; sigma_v = data[6]; sigma_w = data[7];
	if (abs(L) > 200)
	  hABL = 0.1*uStar/fCori;
	else if ( L > 0 && L < 200 )
	  hABL =  0.4*sqrt(uStar*L/fCori);
	else if ( L > -100 || 0.5*(sigma_u+sigma_v) <= uStar*pow(12,1.0/3.0) ) // || (pow(sigma_v/uStar,3.0)-12)*2*(-L) < 100)
        {
          double wStar = 0.5*(sigma_u+sigma_v)/0.6;
	  double buoyancy = g/theta0; // + thetaStar/kappa*(log(zref/zref_T0)-stabKernel(zref/L,"H",true)+stabKernel(zref_T0/L,"H",true)));
          hABL = wStar*wStar*wStar/(-uStar*thetaStar*buoyancy);
	  propConst = sigma_w/pow(zref_wvar,1.0/3.0);
        }
        else
	{
          hABL = (pow(0.5*(sigma_u+sigma_v)/uStar,3.0)-12)*2*(-L);
	  propConst = sigma_w/pow(zref_wvar,1.0/3.0);
	}
     }
     else if ( abs(L) > 200 )
     {
	hABL = 0.1*uStar/fCori;
	sigma_u = 1.2*sigma_v;
	if ( sigmaEl == -1 )
	  sigma_w = 1.3*uStar;
	else
	  sigma_w = uref*sigmaEl*PI/180;
     }
     else if ( L < 0 && L > -200 )
     {
	sigma_u = sigma_v;
        double buoyancy = g/theta0; // + thetaStar/kappa*(log(zref/zref_T0)-stabKernel(zref/L,"H",true)+stabKernel(zref_T0/L,"H",true)));
	if ( sigmaEl == -1 )
          propConst = 1.4*pow(-uStar*thetaStar*buoyancy,1.0/3.0);
	else
	{
	  sigma_w = uref*sigmaEl*PI/180;
	  propConst = sigma_w/pow(zref_wvar,1.0/3.0);
	}

	if (L > -100 || sigma_v <= uStar*pow(12,1.0/3.0)) // || (pow(sigma_v/uStar,3.0)-12)*2*(-L) < 100)
        {
	  double wStar = sigma_v/0.6;
          hABL = wStar*wStar*wStar/(-uStar*thetaStar*buoyancy);
	}
	else
	  hABL = (pow(sigma_v/uStar,3.0)-12)*2*(-L); 
     }
     else
     {
	hABL =  0.4*sqrt(uStar*L/fCori); // (7.2921e-5*sin(latitude*PI/180)));
	sigma_u = 1.5*sigma_v;
	if ( sigmaEl == -1 )
	  sigma_w = 1.6*uStar;
	else
	  sigma_w = uref*sigmaEl*PI/180;
     }

   if (useLangevin)
   {
     ofstream Kfile("output/PARAMS_KTaylorCF_" + meteofile_, ios::app);
     if (Kfile.is_open())
     {
        Kfile << left << setw(width) << timeIndex
              << left << setw(width) << sigma_u
	      << left << setw(width) << sigma_v 
	      << left << setw(width) << sigma_w
              << left << setw(width) << propConst
	      << left << setw(width) << uStar 
              << left << setw(width) << L
	      << left << setw(width) << heff 
              << left << setw(width) << hABL << "\n";
        Kfile.close();
     }
   }
   else
   {
     ofstream Kfile("output/PARAMS_KCouetteFlow_" + meteofile_, ios::app);
     if (Kfile.is_open())
     {
        Kfile << left << setw(width) << timeIndex
              << left << setw(width) << sigma_u
              << left << setw(width) << sigma_v
              << left << setw(width) << sigma_w
              << left << setw(width) << L
              << left << setw(width) << propConst << "\n";
/*
              << left << setw(width) << timeStep[0]
              << left << setw(width) << timeStep[2] 
              << left << setw(width) << tau[0]
              << left << setw(width) << tau[1]
              << left << setw(width) << tau[2] << "\n";
*/  
      Kfile.close();
     }
   }
 }

 double CouetteFlow::stabKernel(const double& z, const string& type, const bool& integrate)
 {
   if (integrate)
   {
     if (z > 0)
       return -5*z;
     else if ( type.compare("H") == 0 )
     {
       double x = pow(1-16*z,0.25);
       x = 0.5*(1+x*x);
       return 2*log(x);
     }
     else if ( type.compare("M") == 0 )
     {
       double sqrt_x = sqrt(1-16*z), x = sqrt(sqrt_x);
       double atan_x = atan(x);
       x = 0.5*(1+sqrt_x)*0.25*(1+2*x+sqrt_x);
       return log(x)-2*atan_x+0.5*PI;
     }
     else
     {
	cout << "\ntype " << type << " in function 'CouetteFlow::stabKernel ( cont double&, const string&, const bool& )' is not implemented.\n" << endl;
	exit(EXIT_FAILURE);
     }
   }
   else
   {
     if (z > 0)
       return 1.0+5.0*z;
     else if ( type.compare("H") == 0 )
       return pow(1.0-16.0*z,-0.5);
     else if ( type.compare("M") == 0 )
       return pow(1.0-16.0*z,-0.25);
     else
     {
        cout << "\ntype " << type << " in function 'CouetteFlow::stabKernel ( cont double&, const string&, const bool& )' is not implemented.\n" << endl;
        exit(EXIT_FAILURE);
     }
   }
 }

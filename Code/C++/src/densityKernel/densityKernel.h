#ifndef DENSITYKERNEL_H
#define DENSITYKERNEL_H

 #include <string>
 using namespace std;

 namespace densityKernel
 {
    void bandWidth
    (
	double**,
	double*,
	int,
	unsigned long,
	bool
    );

    void bandWidth3DKernel
    ( 
	double**, 
	double*, 
	int, 
	unsigned long 
    );

    void bandWidthProdKernel
    (
        double**,
        double*,
        int,
        unsigned long
    );

    void bandWidthProdEpanechnikov
    (
        double**,
        double*,
        int,
        unsigned long
    );

    void bandWidthGaussian
    (
        double**,
        double*,
        int,
        unsigned long
    );

    double weight
    ( 
	string 
    );

    unsigned fac
    (
	unsigned
    );

    double integralKernel
    ( 
	bool 
    );

    double integral3DKernel
    (
    );

   double integralProdKernel
   ( 
   );

   double integralProdEpanechnikov
   (
   );

   double integral3DEpanechnikov
   (
   );

   double integralGaussian
   (
   );

  };

#endif


#define Pi (4.*atan(1.))
#include "structures.h"
#include "functions.h"
#include "pt_kernels.h"
#include "tns_kernels.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Function_A11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

        double r=X[0];
        double x=X[1];

        double Pk=P_interpolLOG(k0,k_L,Pk_L,N);
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pk=P_interpol(k0,k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

if( fabs(1+r*r-2.*r*x)<1e-10 )
{
fval[0]=pow(k0,3)*pow(2.*Pi,-2)*0.5*Pk*a11(r,x)*Pkr;
}
else
{
   fval[0]=pow(k0,3)*pow(2.*Pi,-2)*( (A11(r,x)*Pk+A11t(r,x)*Pkr)*Pksqrt*pow(1+r*r-2.*r*x,-2)+ 0.5*Pk*a11(r,x)*Pkr  );

}

}

void Function_A12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pk=P_interpolLOG(k0,k_L,Pk_L,N);
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pk=P_interpol(k0,k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
        fval[0]=pow(k0,3)*pow(2.*Pi,-2)*0.5*Pk*a12(r,x)*Pkr;
	}
	else
	{
        fval[0]=pow(k0,3)*pow(2.*Pi,-2)*( (A12(r,x)*Pk+A12t(r,x)*Pkr)*Pksqrt*pow(1+r*r-2.*r*x,-2) + 0.5*Pkr*Pk*a12(r,x) );

        }
}

void Function_A22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pk=P_interpolLOG(k0,k_L,Pk_L,N);
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pk=P_interpol(k0,k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
	fval[0]=pow(k0,3)*pow(2.*Pi,-2)*0.5*Pk*a22(r,x)*Pkr;
	}
	else
	{
        fval[0]=pow(k0,3)*pow(2.*Pi,-2)*( (A22(r,x)*Pk+A22t(r,x)*Pkr)*Pksqrt*pow(1+r*r-2.*r*x,-2) + 0.5*Pkr*Pk*a22(r,x) );//+0.5*Pk*a22(r)*Pkr;
	}
}

void Function_A23(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pk=P_interpolLOG(k0,k_L,Pk_L,N);
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pk=P_interpol(k0,k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
	fval[0]=+pow(k0,3)*pow(2.*Pi,-2)*0.5*Pk*a23(r,x)*Pkr;
	}
	else
	{
          fval[0]=pow(k0,3)*pow(2.*Pi,-2)*( (A23(r,x)*Pk+A23t(r,x)*Pkr)*Pksqrt*pow(1+r*r-2.*r*x,-2) + 0.5*Pkr*Pk*a23(r,x) );//+0.5*Pk*a23(r)*Pkr;

	}
}

void Function_A33(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pk=P_interpolLOG(k0,k_L,Pk_L,N);
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pk=P_interpol(k0,k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
	fval[0]=pow(k0,3)*pow(2.*Pi,-2)*0.5*Pk*a33(r,x)*Pkr;
	}
	else
	{
      fval[0]=pow(k0,3)*pow(2.*Pi,-2)*( (A33(r,x)*Pk+A33t(r,x)*Pkr)*Pksqrt*pow(1+r*r-2.*r*x,-2) + 0.5*Pkr*Pk*a33(r,x) );//+0.5*Pk*a33(r)*Pkr;

	}
}

void Function_B1_11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+1)*pow(k0,3)*pow(2.*Pi,-2)*( B1_11(r,x)*pow(1.+r*r-2.*r*x,-1)*Pkr*Pksqrt );
	}
}

void Function_B1_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+2)*pow(k0,3)*pow(2.*Pi,-2)*( B1_12(r,x)*pow(1.+r*r-2.*r*x,-1)*Pkr*Pksqrt );
	}
}

void Function_B1_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,2+1)*pow(k0,3)*pow(2.*Pi,-2)*( B1_21(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );
	}
}


void Function_B1_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];
	
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,2+2)*pow(k0,3)*pow(2.*Pi,-2)*( B1_22(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );
	}
}


void Function_B2_11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];
        
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+1)*pow(k0,3)*pow(2.*Pi,-2)*( B2_11(r,x)*pow(1.+r*r-2.*r*x,-1)*Pkr*Pksqrt );
	}
}


void Function_B2_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];
	
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+2)*pow(k0,3)*pow(2.*Pi,-2)*( B2_12(r,x)*pow(1.+r*r-2.*r*x,-1)*Pkr*Pksqrt );
	}
}

void Function_B2_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];
        
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,2+1)*pow(k0,3)*pow(2.*Pi,-2)*( B2_21(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );
	}
}

void Function_B2_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,2+2)*pow(k0,3)*pow(2.*Pi,-2)*( B2_22(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );

	}
}

void Function_B3_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+2)*pow(k0,3)*pow(2.*Pi,-2)*( B3_12(r,x)*pow(1.+r*r-2.*r*x,-1)*Pkr*Pksqrt );
	}
}

void Function_B3_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];
        
        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,1+2)*pow(k0,3)*pow(2.*Pi,-2)*( B3_21(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );
	}
}

void Function_B3_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);

	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{

fval[0]=pow(-1,2+2)*pow(k0,3)*pow(2.*Pi,-2)*( B3_22(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );

	}
}

void Function_B4_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params

        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L_2pi3;
        int N= params_function.Nlin;

	double r=X[0];
	double x=X[1];

        double Pkr=P_interpolLOG(k0*r,k_L,Pk_L,N);
        double Pksqrt=P_interpolLOG(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
//        double Pkr=P_interpol(k0*r,k_L,Pk_L,N);
//        double Pksqrt=P_interpol(k0*sqrt(1.+r*r-2.*r*x),k_L,Pk_L,N);
	if( fabs(1+r*r-2.*r*x)<1e-10 )
	{
		fval[0]=0;
	}
	else
	{
 fval[0]=pow(-1,2+2)*pow(k0,3)*pow(2.*Pi,-2)*( B4_22(r,x)*pow(1.+r*r-2.*r*x,-2)*Pkr*Pksqrt );

	}
}



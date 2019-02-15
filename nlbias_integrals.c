#define Pi (4.*atan(1.))
#include "structures.h"
#include "functions.h"
#include "pt_kernels.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void intPb2_delta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpolLOG(k_q,k_lin,Pk_lin,N)*P_interpolLOG(q,k_lin,Pk_lin,N)*F2(q,k_q,(k*ctheta-q)/k_q);}
//        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpol(k_q,k_lin,Pk_lin,N)*P_interpol(q,k_lin,Pk_lin,N)*F2(q,k_q,(k*ctheta-q)/k_q);}


}


void intPb2_theta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpolLOG(k_q,k_lin,Pk_lin,N)*P_interpolLOG(q,k_lin,Pk_lin,N)*G2(q,k_q,(k*ctheta-q)/k_q);}
//        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpol(k_q,k_lin,Pk_lin,N)*P_interpol(q,k_lin,Pk_lin,N)*G2(q,k_q,(k*ctheta-q)/k_q);}


}


void intPbs2_delta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpolLOG(k_q,k_lin,Pk_lin,N)*P_interpolLOG(q,k_lin,Pk_lin,N)*F2(q,k_q,(k*ctheta-q)/k_q)*S2((k*ctheta-q)/k_q);}
//        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpol(k_q,k_lin,Pk_lin,N)*P_interpol(q,k_lin,Pk_lin,N)*F2(q,k_q,(k*ctheta-q)/k_q)*S2((k*ctheta-q)/k_q);}

}

void intPbs2_theta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpolLOG(k_q,k_lin,Pk_lin,N)*P_interpolLOG(q,k_lin,Pk_lin,N)*G2(q,k_q,(k*ctheta-q)/k_q)*S2((k*ctheta-q)/k_q);}
//        else{fval[0]=pow(2.*Pi,-2)*(q*q)*P_interpol(k_q,k_lin,Pk_lin,N)*P_interpol(q,k_lin,Pk_lin,N)*G2(q,k_q,(k*ctheta-q)/k_q)*S2((k*ctheta-q)/k_q);}

}




void intPb2s2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpolLOG(q,k_lin,Pk_lin,N)*(2./3.*P_interpolLOG(q,k_lin,Pk_lin,N)-P_interpolLOG(k_q,k_lin,Pk_lin,N)*S2((k*ctheta-q)/k_q));}
//        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpol(q,k_lin,Pk_lin,N)*(2./3.*P_interpol(q,k_lin,Pk_lin,N)-P_interpol(k_q,k_lin,Pk_lin,N)*S2((k*ctheta-q)/k_q));}

}



void intPbs22(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpolLOG(q,k_lin,Pk_lin,N)*(4./9.*P_interpolLOG(q,k_lin,Pk_lin,N)-P_interpolLOG(k_q,k_lin,Pk_lin,N)*pow(S2((k*ctheta-q)/k_q),2));}
//        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpol(q,k_lin,Pk_lin,N)*(4./9.*P_interpol(q,k_lin,Pk_lin,N)-P_interpol(k_q,k_lin,Pk_lin,N)*pow(S2((k*ctheta-q)/k_q),2));}

}


void intPb22(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
 //   double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpolLOG(q,k_lin,Pk_lin,N)*(P_interpolLOG(q,k_lin,Pk_lin,N)-P_interpolLOG(k_q,k_lin,Pk_lin,N));}
//        else{fval[0]=-0.5*pow(2.*Pi,-2)*(q*q)*P_interpol(q,k_lin,Pk_lin,N)*(P_interpol(q,k_lin,Pk_lin,N)-P_interpol(k_q,k_lin,Pk_lin,N));}

}



void intsigma3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;

    double k=params_function.k;
//    double muk=params_function.muk;
    double *Pk_lin=params_function.PK_L_2pi3;
    double *k_lin=params_function.K_L;
    int N=params_function.Nlin;



        double q=x[0];

        double ctheta=x[1];
        double k_q=sqrt(k*k+q*q-2.*q*k*ctheta);

        if(k_q<k_lin[0] || k_q>k_lin[N-1]){fval[0]=0;}
        else{fval[0]=pow(2.*Pi,-2)*q*q*P_interpolLOG(q,k_lin,Pk_lin,N)*(5./6.+15./8.*S2((k*ctheta-q)/k_q)*S2(-ctheta)-5./4.*S2((k*ctheta-q)/k_q));}
//        else{fval[0]=pow(2.*Pi,-2)*q*q*P_interpol(q,k_lin,Pk_lin,N)*(5./6.+15./8.*S2((k*ctheta-q)/k_q)*S2(-ctheta)-5./4.*S2((k*ctheta-q)/k_q));}


}




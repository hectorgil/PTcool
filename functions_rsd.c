#define Pi (4.*atan(1.))
#include "structures.h"
#include "functions.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void sigma2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

f_params params_function = *(f_params *) fdata;

double k0=params_function.k;
double *k_L=params_function.K_L;
double *Pk_L=params_function.PK_L;
int N=params_function.Nlin;



double q=x[0];

fval[0]=P_interpolLOG(q,k_L,Pk_L,N);
//fval[0]=P_interpol(q,k_L,Pk_L,N);


}

double min(int *inicio, double k1,int *a1,double k2, int *a2,double k3, int *a3,double k4, int *a4,double k5,int *a5, int n)
{

        if (a1==0 && a2==0 && a3==0 && a4==0 && a5==0){printf("Error 1 en min1");}


        double f;

        if(*inicio==1)
        {

        if(*a1==0){f=k1;*a1=n;}
        else
        {
                if (*a2==0){f=k2;*a2=n;}
                else
                {
                        if(*a3==0){f=k3;*a3=n;}
                        else
                        {
                                if (*a4==0){f=k4;*a4=n;}
                                else
                                {
                                        f=k5;*a5=n;
                                }

                        }

                }

        }

        }
     if (*inicio==2)
        {

        if(*a2==0){f=k2;*a2=n;}
        else
        {
                if (*a3==0){f=k3;*a3=n;}
                else
                {
                        if(*a4==0){f=k4;*a4=n;}
                        else
                        {
                                if (*a5==0){f=k5;*a5=n;}
                                else
                                {
                                        f=k1;*a1=n;
                                }

                        }

                }

        }

  }
        if (*inicio==3)
        {

        if(*a3==0){f=k3;*a3=n;}
        else
        {
                if (*a4==0){f=k4;*a4=n;}
                else
                {
                        if(*a5==0){f=k5;*a5=n;}
                        else
                        {
                                if (*a1==0){f=k1;*a1=n;}
                                else
                                {
                                        f=k2;*a2=n;
                                }

                        }

                }

        }

}
        if (*inicio==4)
        {

        if(*a4==0){f=k4;*a4=n;}
        else
        {
                if (*a5==0){f=k5;*a5=n;}
                else
                {
                        if(*a1==0){f=k1;*a1=n;}
                        else
                        {
                                if (*a2==0){f=k2;*a2=n;}
                                else
                                {
                                        f=k3;*a3=n;
                                }

                        }

                }

        }
}
       if (*inicio==5)
        {

        if(*a5==0){f=k5;*a5=n;}
        else
        {
                if (*a1==0){f=k1;*a1=n;}
                else
                {
                        if(*a2==0){f=k2;*a2=n;}
                        else
                        {
                                if (*a3==0){f=k3;*a3=n;}
                                else
                                {
                                        f=k4;*a4=n;
                                }

                        }

                }

        }

        }
        if (*a1==n){*inicio=2;}
        if (*a2==n){*inicio=3;}
    if (*a3==n){*inicio=4;}
    if (*a4==n){*inicio=5;}
        if (*a5==n){*inicio=1;}

        return f;


}

double N2(double pl,double p13,double p15)
{
       double f;
    if (p15>=0)
    {
                f=cosh(sqrt(p15/pl))+0.5*p13/pl*sqrt(pl/p15)*sinh(sqrt(p15/pl));
    }
        else
        {
        f=cos(sqrt(-p15/pl))+0.5*p13/pl*sqrt(-pl/p15)*sin(sqrt(-p15/pl));
        }

        f=f*f;
        return f;

}

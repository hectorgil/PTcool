#define Pi (4.*atan(1.))
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double A11 (double r, double x)
{
double f;
f=-pow(r,3)/7.*(x+6.*x*x*x+r*r*x*(-3.+10*x*x)+r*(-3.+x*x-12.*pow(x,4)));
return f;
}

double A12 (double r, double x)
{
double f;
f=pow(r,4)/14.*(x*x-1)*(-1.+7.*r*x-6.*x*x);
return f;
}

double A22 (double r, double x)
{
double f;
f=pow(r,3)/14.*(r*r*x*(13.-41.*x*x)-4.*(x+6.*x*x*x)+r*(5.+9.*x*x+42.*x*x*x*x)    );
return f;
}

double A23 (double r, double x)
{
double f;
f=pow(r,4)/14.*(x*x-1)*(-1.+7.*r*x-6.*x*x);
return f;
}

double A33 (double r, double x)
{
double f;
f=pow(r,3)/14.*(1.-7.*r*x+6.*x*x)*(-2.*x+r*(-1.+3.*x*x));
return f;
}


double A11t ( double r,double x)
{
double f;
f=1./7.*(x+r-2.*r*x*x)*(3.*r+7.*x-10.*r*x*x);
return f;
}

double A12t ( double r,double x)
{
double f;
f=r/14.*(x*x-1.)*(3.*r+7.*x-10.*r*x*x);
return f;
}


double A22t ( double r,double x)
{
double f;
f=1./14.*(28.*x*x+r*x*(25.-81.*x*x)+r*r*(1.-27.*x*x+54.*x*x*x*x ));
return f;
}


double A23t ( double r,double x)
{
double f;
f=r/14.*(1.-x*x)*(r-7.*x+6.*r*x*x);
return f;
}


double A33t ( double r,double x)
{
double f;
f=1./14.*(r-7.*x+6.*r*x*x)*(-2.*x-r+3.*r*x*x);
return f;
}


double a11 ( double r, double x)
{
double f;
if( fabs(r-1.0)<1e-10 )
{
f=-1./(84.*r)*(2.*r*(19.-24.*r*r+9.*r*r*r*r));
}
else
{
f=-1./(84.*r)*(2.*r*(19.-24.*r*r+9.*r*r*r*r)-9.*pow(r*r-1,3)*log(fabs( (r+1)/(r-1))));
}
//f=1./7.*(-7.*x*x+r*r*r*x*(-3.+10*x*x)+3.*r*(x+6.*x*x*x)+r*r*(6.-19*x*x-8.*x*x*x*x));
return f;
}

double a12 ( double r, double x)
{
double f;

if( fabs(r-1.0)<1e-10 )
{
f=1./(112.*r*r*r)*(2.*r*(r*r+1.)*(3.-14.*r*r+3.*r*r*r*r));
}
else
{
f=1./(112.*r*r*r)*(2.*r*(r*r+1.)*(3.-14.*r*r+3.*r*r*r*r)-3*pow(r*r-1,4)*log(fabs( (r+1)/(r-1))));
}

//f=1./14.*r*(-1.+x*x)*(6.*r-7.*(1.+r*r)*x+8*r*x*x);
return f;
}

double a22 ( double r, double x)
{
double f;

if( fabs(r-1.0)<1e-10 )
{
f=1./(336.*r*r*r)*(2.*r*(9.-185*r*r+159*r*r*r*r-63.*r*r*r*r*r*r));
}
else
{
f=1./(336.*r*r*r)*(2.*r*(9.-185*r*r+159*r*r*r*r-63.*r*r*r*r*r*r)+9.*pow(r*r-1,3)*(7.*r*r+1)*log(fabs((r+1)/(r-1))));
}
//f=1./14.*(-28.*x*x+r*r*r*x*(-13.+41.*x*x)+r*x*(11.+73.*x*x)-2.*r*r*(-9.+31.*x*x+20.*x*x*x*x));
return f;
}

double a23 ( double r, double x)
{
double f;
if( fabs(r-1.0)<1e-10 )
{
f=1./(112.*r*r*r)*(2.*r*(r*r+1.)*(3.-14.*r*r+3.*r*r*r*r));
}
else
{
f=1./(112.*r*r*r)*(2.*r*(r*r+1.)*(3.-14.*r*r+3.*r*r*r*r)-3*pow(r*r-1,4)*log(fabs( (r+1)/(r-1))));
}
//f=1./14.*r*(-1.+x*x)*(6.*r-7.*(1.+r*r)*x+8*r*x*x);
return f;
}

double a33 ( double r, double x)
{
double f;
if( fabs(r-1.0)<1e-10 )
{
f=f=1./(336.*r*r*r)*(2.*r*(9.-109*r*r+63*r*r*r*r-27*r*r*r*r*r*r));
}
else
{
f=1./(336.*r*r*r)*(2.*r*(9.-109*r*r+63*r*r*r*r-27*r*r*r*r*r*r)+9*pow(r*r-1,3)*(3.*r*r+1)*log(fabs( (r+1)/(r-1))));
}
//f=1./14.*(7.*x+r*(-6.+7.*r*x-8.*x*x))*(-2.*x+r*(-1+3.*x*x));
return f;
}


double B1_11 ( double r, double x)
{
double f;
f=r*r/2.*(x*x-1.);
return f;
}

double B1_12 ( double r, double x)
{
double f;
f=3.*r*r/8.*pow(x*x-1,2);
return f;
}

double B1_21 ( double r, double x)
{
double f;
f=3.*r*r*r*r/8.*pow(x*x-1,2);
return f;
}

double B1_22 ( double r, double x)
{
double f;
f=5.*r*r*r*r/16.*pow(x*x-1,3);
return f;
}

double B2_11 ( double r, double x)
{
double f;
f=r/2.*(r+2.*x-3.*r*x*x);
return f;
}

double B2_12 ( double r, double x)
{
double f;
f=-3.*r/4.*(x*x-1)*(-r-2.*x+5.*r*x*x);
return f;
}


double B2_21 ( double r, double x)
{
double f;
f=3.*r*r/4.*(x*x-1)*(-2.+r*r+6.*r*x-5.*r*r*x*x);
return f;
}

double B2_22 ( double r, double x)
{
double f;
f=-3.*r*r/16.*pow(x*x-1,2)*(6.-30.*r*x-5.*r*r+35*r*r*x*x);
return f;
}

double B3_12 ( double r, double x)
{
double f;
f=r/8.*( 4.*x*(3.-5.*x*x)+r*(3.-30.*x*x+35*x*x*x*x));
return f;
}

double B3_21 ( double r, double x)
{
double f;
f=r/8.*(-8.*x+r*(-12.+36.*x*x+12*r*x*(3-5.*x*x)+r*r*(3.-30.*x*x+35*x*x*x*x)));
return f;
}


double B3_22 ( double r, double x)
{
double f;
f=3.*r/16.*(x*x-1)*(-8.*x+r*(-12.+60.*x*x+20.*r*x*(3-7.*x*x)+5.*r*r*(1.-14.*x*x+21*x*x*x*x)));
return f;
}

double B4_22 ( double r, double x)
{
double f;
f=r/16.*(8.*x*(-3.+5.*x*x)-6.*r*(3.-30.*x*x+35*x*x*x*x)+6.*r*r*x*(15.-70.*x*x+63*x*x*x*x)+r*r*r*(5.-21.*x*x*(5.-15.*x*x+11.*x*x*x*x)));
return f;
}


//What do these terms mean?
/*

        A11=BIAS*BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_A11);
        A12=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_A12);
        A22=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_A22);
        A23=interpol_taruya(Kmodes[i],k_taruya,Taruya_A23);
        A33=interpol_taruya(Kmodes[i],k_taruya,Taruya_A33);

        B1_11=BIAS*BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B1_11);
        B1_12=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B1_12);
        B1_21=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B1_21);
        B1_22=interpol_taruya(Kmodes[i],k_taruya,Taruya_B1_22);
        B2_11=BIAS*BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B2_11);
        B2_12=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B2_12);
        B2_21=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B2_21);
        B2_22=interpol_taruya(Kmodes[i],k_taruya,Taruya_B2_22);
        B3_12=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B3_12);
        B3_21=BIAS*interpol_taruya(Kmodes[i],k_taruya,Taruya_B3_21);
        B3_22=interpol_taruya(Kmodes[i],k_taruya,Taruya_B3_22);
        B4_22=interpol_taruya(Kmodes[i],k_taruya,Taruya_B4_22);
        A_taruya=pow(MUmodes[i],2)*pow(a[1],1)*A11+pow(MUmodes[i],2)*pow(a[1],2)*A12+pow(MUmodes[i],4)*pow(a[1],2)*A22+pow(MUmodes[i],4)*pow(a[1],3)*A23+pow(MUmodes[i],6)*pow(a[1],3)*A33;
        B_taruya=pow(MUmodes[i],2)*pow(a[1],2)*B1_11+pow(MUmodes[i],2)*pow(a[1],3)*B1_12+pow(MUmodes[i],2)*pow(a[1],3)*B1_21+pow(MUmodes[i],2)*pow(a[1],4)*B1_22+pow(MUmodes[i],4)*pow(a[1],2)*B2_11+pow(MUmodes[i],4)*pow(a[1],3)*B2_12+pow(MUmodes[i],4)*pow(a[1],3)*B2_21+pow(MUmodes[i],4)*pow(a[1],4)*B2_22+pow(MUmodes[i],6)*pow(a[1],3)*B3_12+pow(MUmodes[i],6)*pow(a[1],3)*B3_21+pow(MUmodes[i],6)*pow(a[1],4)*B3_22+pow(MUmodes[i],8)*pow(a[1],4)*B4_22;


*/










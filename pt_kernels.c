#define Pi (4.*atan(1.))
#include "functions_rsd.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Important: Definitions of alpha, beta and kernels according to Scoccimarro 97

double alpha(double q1, double q2, double calpha12) //alpha function from SC97
{
	double f;
        double epsilon=1e-5;
if(q1<epsilon || q2<epsilon)
{
f=0;  //infrared divergence
}
else
{
        f=(q1/q2)*calpha12;      
}
	return f;
}

double beta(double q1,double q2, double q3, double calpha23) //beta function from SC97
{
	double f;
        double epsilon=1e-10;
if(q1<epsilon || q2<epsilon || q3<epsilon)
{
f=0; //infrared divergence
}
else
{
        f=pow(q1,2)*pow(2.*q2*q3,-1)*calpha23;      
}
	return f;
}

double S2(double ctheta)
{
double f;
f=ctheta*ctheta-1./3.;
return f;
}

double F2(double k1, double k2, double calpha12) //2-point kernel (for matter)
{
	double f;
        double epsilon=1e-10;

double mod12=pow(fabs(k1*k1+k2*k2+2.*k1*k2*calpha12),0.5);
double calpha_12_1=(k1+k2*calpha12)/mod12;

if(mod12==0){calpha_12_1=0;}

	f=1./7.*(5.*alpha(mod12,k1,calpha_12_1)+2.*beta(mod12,k1,k2,calpha12));

	return f;
}

double G2(double k1, double k2, double calpha12) //2-point kernel (for velocity)
{
	double f;
	double epsilon=1e-10;

double mod12=pow(fabs(k1*k1+k2*k2+2.*k1*k2*calpha12),0.5);
double calpha_12_1=(k1+k2*calpha12)/mod12;
if(mod12==0){calpha_12_1=0;}

	f=1./7.*(3.*alpha(mod12,k1,calpha_12_1)+4.*beta(mod12,k1,k2,calpha12));
	return f;
}

double F3(double k1, double k2, double k3, double calpha12, double calpha13, double calpha23 ) //3-point kernel (for matter)
{
if( fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha23)>1 ){ printf("Error F3 0.0 cos(alpha12)=%lf cos(alpha13)=%lf cos(alpha23)=%lf\n", calpha12, calpha13, calpha23); /* return 0; */}

double f1,f2,f;

double mod123, mod12, mod23; //mod123 es el modulo de q1+q2+q3, mod12 es el moduloa de q1+q2, mod23 de q2+q3

double calpha123_1, calpha23_1, calpha123_12, calpha12_3; //alpha123_1 es el angulo entre q1+q2+q3 y q1, alpha23_1 es el angulo entre q2+q3 y q1 etc...

double epsilon=5e-2;
	
mod123=pow(fabs(k1*k1+k2*k2+k3*k3+2.*(k1*k2*calpha12+k1*k3*calpha13+k2*k3*calpha23)),0.5);
mod12=pow(fabs(k1*k1+k2*k2+2.0*k1*k2*calpha12),0.5);
mod23=pow(fabs(k2*k2+k3*k3+2.0*k2*k3*calpha23),0.5);

calpha123_1=1./(mod123)*(k1+k2*calpha12+k3*calpha13);
calpha23_1=(1./mod23)*(k2*calpha12+k3*calpha13);
calpha123_12=(1./(mod12*mod123))*(k1*k1+2.*k1*k2*calpha12+k1*k3*calpha13+k2*k2+k2*k3*calpha23);
calpha12_3=(1./mod12)*(k1*calpha13+k2*calpha23);

if(mod123==0){calpha123_1=0; calpha123_12=0;}
if(mod23==0){calpha23_1=0;}
if(mod12==0){ calpha123_12=0; calpha12_3=0; }

	
//precission correction
if(calpha123_1>1.0 && calpha123_1<1.0+epsilon){ calpha123_1=1.0;}
if(calpha123_1<-1.0 && calpha123_1>-1.0-epsilon){ calpha123_1=-1.0;}

if(calpha23_1>1.0 && calpha23_1<1.0+epsilon){ calpha23_1=1.0;}
if(calpha23_1<-1.0 && calpha23_1>-1.0-epsilon){ calpha23_1=-1.0;}

if(calpha123_12>1.0 && calpha123_12<1.0+epsilon){ calpha123_12=1.0;}
if(calpha123_12<-1.0 && calpha123_12>-1.0-epsilon){ calpha123_12=-1.0;}

if(calpha12_3>1.0 && calpha12_3<1.0+epsilon){ calpha12_3=1.0;}
if(calpha12_3<-1.0 && calpha12_3>-1.0-epsilon){ calpha12_3=-1.0;}


if(fabs(calpha123_1)>1.0){printf("\nError F3 1.0 cos(alpha123_1)=%.10lf \n",calpha123_1); /* return 0; */}
if(fabs(calpha23_1)>1.0){printf("\nError F3 2.0 cos(alpha23_1)=%lf \n",calpha23_1); /* return 0; */}
if(fabs(calpha123_12)>1.0){printf("\nError F3 3.0 cos(alpha123_12)=%lf \n", calpha123_12 ); /* return 0; */}
if(fabs(calpha12_3)>1.0){printf("\nError F3 4.0 cos(alpha12_3)=%lf \n",calpha12_3); /* return 0; */}

         //m=1
        f1=(2.*3.+1.)*alpha(mod123,k1,calpha123_1)*F2(k2,k3,calpha23)+2.*beta(mod123,k1,mod23,calpha23_1)*G2(k2,k3,calpha23);
        
        //m=2
        f2=G2(k1,k2,calpha12)*( (2.*3.+1.)*alpha(mod123,mod12,calpha123_12)+2.*beta(mod123,mod12,k3,calpha12_3)  );

         f=pow((2.*3.+3.)*(3.-1.),-1)*(f1+f2);



	return f;
}

//non-symetryzed-4-points kernel (for velocity fields)
double G3(double k1,double k2, double k3, double calpha12, double calpha13, double calpha23) 
{

	if( fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha23)>1 ){ printf("Error G3 0.0 cos(alpha12)=%lf cos(alpha13)=%lf cos(alpha23)=%lf\n", calpha12, calpha13, calpha23);}
	
	double f1,f2,f;
	
	double mod123, mod12, mod23; //mod123 es el modulo de q1+q2+q3, mod12 es el moduloa de q1+q2, mod23 de q2+q3
	
	double calpha123_1, calpha23_1, calpha123_12, calpha12_3; //alpha123_1 es el angulo entre q1+q2+q3 y q1, alpha23_1 es el angulo entre q2+q3 y q1 etc...
	
	double epsilon=5e-2;
	
	mod123=pow(fabs(k1*k1+k2*k2+k3*k3+2.*(k1*k2*calpha12+k1*k3*calpha13+k2*k3*calpha23)),0.5);
	mod12=pow(fabs(k1*k1+k2*k2+2.0*k1*k2*calpha12),0.5);
	mod23=pow(fabs(k2*k2+k3*k3+2.0*k2*k3*calpha23),0.5);
	
	calpha123_1=1./(mod123)*(k1+k2*calpha12+k3*calpha13);
	calpha23_1=(1./mod23)*(k2*calpha12+k3*calpha13);
	calpha123_12=(1./(mod12*mod123))*(k1*k1+2.*k1*k2*calpha12+k1*k3*calpha13+k2*k2+k2*k3*calpha23);
	calpha12_3=(1./mod12)*(k1*calpha13+k2*calpha23);
	
	if(mod123==0){calpha123_1=0; calpha123_12=0;}
	if(mod23==0){calpha23_1=0;}
	if(mod12==0){ calpha123_12=0; calpha12_3=0; }

	
	//precission correction
	if(calpha123_1>1.0 && calpha123_1<1.0+epsilon){ calpha123_1=1.0;}
	if(calpha123_1<-1.0 && calpha123_1>-1.0-epsilon){ calpha123_1=-1.0;}
	
	if(calpha23_1>1.0 && calpha23_1<1.0+epsilon){ calpha23_1=1.0;}
	if(calpha23_1<-1.0 && calpha23_1>-1.0-epsilon){ calpha23_1=-1.0;}
	
	if(calpha123_12>1.0 && calpha123_12<1.0+epsilon){ calpha123_12=1.0;}
	if(calpha123_12<-1.0 && calpha123_12>-1.0-epsilon){ calpha123_12=-1.0;}
	
	if(calpha12_3>1.0 && calpha12_3<1.0+epsilon){ calpha12_3=1.0;}
	if(calpha12_3<-1.0 && calpha12_3>-1.0-epsilon){ calpha12_3=-1.0;}
	
	
	if(fabs(calpha123_1)>1.0){printf("\nError G3 1.0 cos(alpha123_1)=%.10lf (%lf,%lf,%lf) (%lf,%lf,%lf) \n",calpha123_1,k1,k2,k3,calpha12,calpha13,calpha23); }
	if(fabs(calpha23_1)>1.0){printf("\nError G3 2.0 cos(alpha23_1)=%lf (%lf,%lf,%lf) (%lf,%lf,%lf) \n",calpha23_1,k1,k2,k3,calpha12,calpha13,calpha23); }
	if(fabs(calpha123_12)>1.0){printf("\nError G3 3.0 cos(alpha123_12)=%lf (%lf,%lf,%lf) (%lf,%lf,%lf) \n", calpha123_12,k1,k2,k3,calpha12,calpha13,calpha23);}
	if(fabs(calpha12_3)>1.0){printf("\nError G3 4.0 cos(alpha12_3)=%lf (%lf,%lf,%lf) (%lf,%lf,%lf)\n",calpha12_3,k1,k2,k3,calpha12,calpha13,calpha23);}
	
	
 
/* m=1 */     f1=3.*alpha(mod123,k1,calpha123_1)*F2(k2,k3,calpha23)+6.*beta(mod123,k1,mod23,calpha23_1)*G2(k2,k3,calpha23);

/* m=2 */     f2=G2(k1,k2,calpha12)*( 3.*alpha(mod123,mod12,calpha123_12)+6.*beta(mod123,mod12,k3,calpha12_3) );

              f=pow((2.*3.+3.)*(3.-1.),-1)*(f1+f2);

	      return f;
}


double F4(double k1, double k2, double k3, double k4, double calpha12, double calpha13, double calpha14, double calpha23, double calpha24, double calpha34)
{
double f;
double f1,f2,f3;

double epsilon=5e-2;

double mod1234, mod234, mod12, mod34, mod123;
	
	if (fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha14)>1 || fabs(calpha23)>1 || fabs(calpha24)>1 || fabs(calpha34)>1){printf("Error en la entrada de F4\n");}

mod1234=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
mod234=pow(fabs(k2*k2 + k3*k3 + k4*k4 + 2.*(k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
mod12=pow(fabs(k1*k1 + k2*k2 + 2.*k1*k2*calpha12),0.5);
mod34=pow(fabs(k3*k3 + k4*k4 + 2.*k3*k4*calpha34),0.5);
mod123=pow(fabs(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23)),0.5);


double calpha234_1, calpha12_34, calpha123_4, calpha1234_123, calpha1234_12, calpha1234_1;

calpha234_1=1./mod234*(k2*calpha12 + k3*calpha13 + k4*calpha14);

calpha12_34=1./(mod12*mod34)*(k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24);

calpha123_4=1./mod123*(k1*calpha14 + k2*calpha24 + k3*calpha34);

calpha1234_1=1./mod1234*(k1 + k2*calpha12 + k3*calpha13 + k4*calpha14);

calpha1234_12=1./(mod1234*mod12)*(k1*k1 + k2*k2 + 2.*k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23 + k1*k4*calpha14 + k2*k4*calpha24);

calpha1234_123=1./(mod1234*mod123)*(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23) + k1*k4*calpha14 + k2*k4*calpha24 + k3*k4*calpha34);

if(mod1234==0){calpha1234_1=0; calpha1234_12=0; calpha1234_123=0;}
if(mod234==0){calpha234_1=0;}
if(mod12==0){ calpha12_34=0;calpha1234_12=0; }
if(mod34==0){calpha12_34=0;}
if(mod123==0){calpha1234_123=0; calpha123_4=0;}

//precission correction
if(calpha234_1>1.0 && calpha234_1<1.0+epsilon){ calpha234_1=1.0;}
if(calpha234_1<-1.0 && calpha234_1>-1.0-epsilon){ calpha234_1=-1.0;}

if(calpha12_34>1.0 && calpha12_34<1.0+epsilon){ calpha12_34=1.0;}
if(calpha12_34<-1.0 && calpha12_34>-1.0-epsilon){ calpha12_34=-1.0;}

if(calpha123_4>1.0 && calpha123_4<1.0+epsilon){ calpha123_4=1.0;}
if(calpha123_4<-1.0 && calpha123_4>-1.0-epsilon){ calpha123_4=-1.0;}

if(calpha1234_1>1.0 && calpha1234_1<1.0+epsilon){ calpha1234_1=1.0;}
if(calpha1234_1<-1.0 && calpha1234_1>-1.0-epsilon){ calpha1234_1=-1.0;}

if(calpha1234_12>1.0 && calpha1234_12<1.0+epsilon){ calpha1234_12=1.0;}
if(calpha1234_12<-1.0 && calpha1234_12>-1.0-epsilon){ calpha1234_12=-1.0;}

if(calpha1234_123>1.0 && calpha1234_123<1.0+epsilon){ calpha1234_123=1.0;}
if(calpha1234_123<-1.0 && calpha1234_123>-1.0-epsilon){ calpha1234_123=-1.0;}

if(fabs(calpha234_1)>1.0){printf("\nError F4 1.0 cos(alpha234_1)=%.20lf \n",calpha234_1); /* return 0; */}
if(fabs(calpha12_34)>1.0){printf("\nError F4 2.0 cos(alpha12_34)=%.20lf \n",calpha12_34); /* return 0; */}
if(fabs(calpha123_4)>1.0){printf("\nError F4 3.0 cos(alpha123_4)=%.20lf \n",calpha123_4); /* return 0; */}
if(fabs(calpha1234_1)>1.0){printf("\nError F4 4.0 cos(alpha1234_1)=%.20lf \n",calpha1234_1); /* return 0; */}
if(fabs(calpha1234_12)>1.0){printf("\nError F4 5.0 cos(alpha1234_12)=%.20lf \n",calpha1234_12); /* return 0; */}
if(fabs(calpha1234_123)>1.0){printf("\nError F4 6.0 cos(alpha1234_123)=%.20lf \n",calpha1234_123); /* return 0; */}

/*m=1*/ f1=(2.*4.+1.)*alpha(mod1234,k1,calpha1234_1)*F3(k2,k3,k4,calpha23,calpha24,calpha34)+2.*beta(mod1234,k1,mod234,calpha234_1)*G3(k2,k3,k4,calpha23,calpha24,calpha34);

/*m=2*/ f2=G2(k1,k2,calpha12)*((2.*4.+1)*alpha(mod1234,mod12,calpha1234_12)*F2(k3,k4,calpha34)+2.*beta(mod1234,mod12,mod34,calpha12_34)*G2(k3,k4,calpha34));

/*m=3*/ f3=G3(k1,k2,k3,calpha12,calpha13,calpha23)*((2.*4.+1)*alpha(mod1234,mod123,calpha1234_123)*1.+2.*beta(mod1234,mod123,k4,calpha123_4)*1.);

f=pow((2.*4.+3.)*(4.-1.),-1)*(f1+f2+f3);

return f;
}


double G4(double k1, double k2, double k3, double k4, double calpha12, double calpha13, double calpha14, double calpha23, double calpha24, double calpha34)
{
	double f;
	double f1,f2,f3;
	
	double epsilon=5e-2;
	
	if (fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha14)>1 || fabs(calpha23)>1 || fabs(calpha24)>1 || fabs(calpha34)>1){printf("Error en la entrada de G4\n");}
	
	double mod1234, mod234, mod12, mod34, mod123;
	
	mod1234=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	mod234=pow(fabs(k2*k2 + k3*k3 + k4*k4 + 2.*(k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	mod12=pow(fabs(k1*k1 + k2*k2 + 2.*k1*k2*calpha12),0.5);
	mod34=pow(fabs(k3*k3 + k4*k4 + 2.*k3*k4*calpha34),0.5);
	mod123=pow(fabs(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23)),0.5);
	
	
	double calpha234_1, calpha12_34, calpha123_4, calpha1234_123, calpha1234_12, calpha1234_1;
	
	calpha234_1=1./mod234*(k2*calpha12 + k3*calpha13 + k4*calpha14);
	
	calpha12_34=1./(mod12*mod34)*(k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24);
	
	calpha123_4=1./mod123*(k1*calpha14 + k2*calpha24 + k3*calpha34);
	
	calpha1234_1=1./mod1234*(k1 + k2*calpha12 + k3*calpha13 + k4*calpha14);
	
	calpha1234_12=1./(mod1234*mod12)*(k1*k1 + k2*k2 + 2.*k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23 + k1*k4*calpha14 + k2*k4*calpha24);
	
	calpha1234_123=1./(mod1234*mod123)*(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23) + k1*k4*calpha14 + k2*k4*calpha24 + k3*k4*calpha34);
	
	if(mod1234==0){calpha1234_1=0; calpha1234_12=0; calpha1234_123=0;}
	if(mod234==0){calpha234_1=0;}
	if(mod12==0){ calpha12_34=0;calpha1234_12=0; }
	if(mod34==0){calpha12_34=0;}
	if(mod123==0){calpha1234_123=0; calpha123_4=0;}
	
	//precission correction
	if(calpha234_1>1.0 && calpha234_1<1.0+epsilon){ calpha234_1=1.0;}
	if(calpha234_1<-1.0 && calpha234_1>-1.0-epsilon){ calpha234_1=-1.0;}
	
	if(calpha12_34>1.0 && calpha12_34<1.0+epsilon){ calpha12_34=1.0;}
	if(calpha12_34<-1.0 && calpha12_34>-1.0-epsilon){ calpha12_34=-1.0;}
	
	if(calpha123_4>1.0 && calpha123_4<1.0+epsilon){ calpha123_4=1.0;}
	if(calpha123_4<-1.0 && calpha123_4>-1.0-epsilon){ calpha123_4=-1.0;}
	
	if(calpha1234_1>1.0 && calpha1234_1<1.0+epsilon){ calpha1234_1=1.0;}
	if(calpha1234_1<-1.0 && calpha1234_1>-1.0-epsilon){ calpha1234_1=-1.0;}
	
	if(calpha1234_12>1.0 && calpha1234_12<1.0+epsilon){ calpha1234_12=1.0;}
	if(calpha1234_12<-1.0 && calpha1234_12>-1.0-epsilon){ calpha1234_12=-1.0;}
	
	if(calpha1234_123>1.0 && calpha1234_123<1.0+epsilon){ calpha1234_123=1.0;}
	if(calpha1234_123<-1.0 && calpha1234_123>-1.0-epsilon){ calpha1234_123=-1.0;}
	
	if(fabs(calpha234_1)>1.0){printf("\nError G4 1.0 cos(alpha234_1)=%.20lf \n",calpha234_1); /* return 0; */}
	if(fabs(calpha12_34)>1.0){printf("\nError G4 2.0 cos(alpha12_34)=%.20lf \n",calpha12_34); /* return 0; */}
	if(fabs(calpha123_4)>1.0){printf("\nError G4 3.0 cos(alpha123_4)=%.20lf \n",calpha123_4); /* return 0; */}
	if(fabs(calpha1234_1)>1.0){printf("\nError G4 4.0 cos(alpha1234_1)=%.20lf \n",calpha1234_1); /* return 0; */}
	if(fabs(calpha1234_12)>1.0){printf("\nError G4 5.0 cos(alpha1234_12)=%.20lf \n",calpha1234_12); /* return 0; */}
	if(fabs(calpha1234_123)>1.0){printf("\nError G4 6.0 cos(alpha1234_123)=%.20lf \n",calpha1234_123); /* return 0; */}
	
	/*m=1*/ f1=3.*alpha(mod1234,k1,calpha1234_1)*F3(k2,k3,k4,calpha23,calpha24,calpha34)+8.*beta(mod1234,k1,mod234,calpha234_1)*G3(k2,k3,k4,calpha23,calpha24,calpha34);
	
	/*m=2*/ f2=G2(k1,k2,calpha12)*(3.*alpha(mod1234,mod12,calpha1234_12)*F2(k3,k4,calpha34)+8.*beta(mod1234,mod12,mod34,calpha12_34)*G2(k3,k4,calpha34));
		
	/*m=3*/ f3=G3(k1,k2,k3,calpha12,calpha13,calpha23)*(3.*alpha(mod1234,mod123,calpha1234_123)*1.+8.*beta(mod1234,mod123,k4,calpha123_4)*1.);
	
	f=pow((2.*4.+3.)*(4.-1.),-1)*(f1+f2+f3);
	
	return f;
}
double F5(double k1, double k2, double k3, double k4, double k5, double calpha12, double calpha13, double calpha14, double calpha15, double calpha23, double calpha24, double calpha25, double calpha34, double calpha35, double calpha45)
{
	double f,f1,f2,f3,f4;
	f1=0;
	f2=0;
	f3=0;
	f4=0;
	
	if (fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha14)>1 || fabs(calpha15)>1 || fabs(calpha23)>1 || fabs(calpha24)>1 || fabs(calpha25)>1 || fabs(calpha34)>1 || fabs(calpha35)>1 || fabs(calpha45)>1 ){printf("Error en la entrada de F5\n");}
	
	double epsilon=5e-2;
	double epsilon2=1e-4;
	double mod12345,mod12,mod123,mod1234, mod2345, mod345, mod45;
	double calpha12345_1,calpha12345_12,calpha12345_123,calpha12345_1234,calpha1_2345,calpha12_345,calpha123_45,calpha1234_5;
	
	mod12345=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + k5*k5 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k3*k4*calpha34 + k3*k5*calpha35 + k4*k5*calpha45)),0.5);
	
	mod12=pow(fabs(k1*k1 + k2*k2 + 2.*k1*k2*calpha12),0.5);
	
	mod123=pow(fabs(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23)),0.5);
	
	mod1234=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	
	mod2345=pow(fabs(k5*k5 + k2*k2 + k3*k3 + k4*k4 + 2.*(k5*k2*calpha25 + k5*k3*calpha35 + k5*k4*calpha45 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	
	mod345=pow(fabs(k3*k3 + k4*k4 + k5*k5 + 2.*(k3*k4*calpha34 + k3*k5*calpha35 + k4*k5*calpha45)),0.5);
	
	mod45=pow(fabs(k4*k4 + k5*k5 + 2.*k4*k5*calpha45),0.5);
	
	calpha12345_1=( k1 + k2*calpha12 + k3*calpha13 + k4*calpha14 + k5*calpha15 )/(mod12345);
	
	calpha12345_12=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25)/(mod12345*mod12);
	
	calpha12345_123=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k1*k3*calpha13 + k2*k3*calpha23 + k3*k3 + k3*k4*calpha34 + k3*k5*calpha35 )/(mod12345*mod123);
	
	calpha12345_1234=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k1*k3*calpha13 + k2*k3*calpha23 + k3*k3 + k3*k4*calpha34 + k3*k5*calpha35 + k4*k1*calpha14 + k4*k2*calpha24 + k4*k3*calpha34 + k4*k4 + k4*k5*calpha45)/(mod12345*mod1234);

	calpha1_2345=(k2*calpha12 + k3*calpha13 + k4*calpha14 + k5*calpha15)/mod2345;
	
	calpha12_345=(k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25)/(mod12*mod345);
	
	calpha123_45=(k1*k4*calpha14 + k2*k4*calpha24 + k4*k3*calpha34 + k1*k5*calpha15 + k2*k5*calpha25 + k3*k5*calpha35)/(mod45*mod123);
	
	calpha1234_5=(k1*calpha15 + k2*calpha25 + k3*calpha35 + k4*calpha45)/(mod1234);
	
	if(mod12345==0){calpha12345_1=0; calpha12345_12=0; calpha12345_123=0;calpha12345_1234=0;}
	if(mod12==0){calpha12345_12=0;calpha12_345=0;}
	if(mod123==0){calpha12345_123=0;calpha123_45=0;}
	if(mod1234==0){calpha12345_1234=0;calpha1234_5=0;}
	if(mod2345==0){calpha1_2345=0;}
	if(mod345==0){calpha12_345=0;}
	if(mod45==0){calpha123_45=0;}
	if(k1==0){calpha12345_1=0;}
	if(k5==0){calpha1234_5=0;}
	
	
	//precision correction
	if(calpha12345_1>1.0 && calpha12345_1<1.0+epsilon){ calpha12345_1=1.0;}
	if(calpha12345_1<-1.0 && calpha12345_1>-1.0-epsilon){ calpha12345_1=-1.0;}

	if(calpha12345_12>1.0 && calpha12345_12<1.0+epsilon){ calpha12345_12=1.0;}
	if(calpha12345_12<-1.0 && calpha12345_12>-1.0-epsilon){ calpha12345_12=-1.0;}
	
	if(calpha12345_123>1.0 && calpha12345_123<1.0+epsilon){ calpha12345_123=1.0;}
	if(calpha12345_123<-1.0 && calpha12345_123>-1.0-epsilon){ calpha12345_123=-1.0;}
	
	if(calpha12345_1234>1.0 && calpha12345_1234<1.0+epsilon){ calpha12345_1234=1.0;}
	if(calpha12345_1234<-1.0 && calpha12345_1234>-1.0-epsilon){ calpha12345_1234=-1.0;}
	
	if(calpha1_2345>1.0 && calpha1_2345<1.0+epsilon){ calpha1_2345=1.0;}
	if(calpha1_2345<-1.0 && calpha1_2345>-1.0-epsilon){ calpha1_2345=-1.0;}
	
	if(calpha12_345>1.0 && calpha12_345<1.0+epsilon){ calpha12_345=1.0;}
	if(calpha12_345<-1.0 && calpha12_345>-1.0-epsilon){ calpha12_345=-1.0;}
	
	if(calpha123_45>1.0 && calpha123_45<1.0+epsilon){ calpha123_45=1.0;}
	if(calpha123_45<-1.0 && calpha123_45>-1.0-epsilon){ calpha123_45=-1.0;}
	
	if(calpha1234_5>1.0 && calpha1234_5<1.0+epsilon){ calpha1234_5=1.0;}
	if(calpha1234_5<-1.0 && calpha1234_5>-1.0-epsilon){ calpha1234_5=-1.0;}
	
	if(fabs(calpha12345_1)>1.0){printf("\nError F5 1.0 cos(calpha12345_1)=%.20lf (%lf,%lf) \n",calpha12345_1,mod12345,k1); }
	if(fabs(calpha12345_12)>1.0){printf("\nError F5 1.0 cos(calpha12345_12)=%.20lf (%lf,%lf) \n",calpha12345_12,mod12345,mod12); }
	if(fabs(calpha12345_123)>1.0){printf("\nError F5 1.0 cos(calpha12345_123)=%.20lf (%lf,%lf) \n",calpha12345_123,mod12345,mod123); }
	if(fabs(calpha12345_1234)>1.0){printf("\nError F5 1.0 cos(calpha12345_1234)=%.20lf (%lf,%lf) \n",calpha12345_1234,mod12345,mod1234); }
	if(fabs(calpha1_2345)>1.0){printf("\nError F5 1.0 cos(calpha1_2345)=%.20lf (%lf,%lf) \n",calpha1_2345,k1,mod2345); }
	if(fabs(calpha12_345)>1.0){printf("\nError F5 1.0 cos(calpha12_345)=%.20lf (%lf,%lf) \n",calpha12_345,mod12,mod345); }
	if(fabs(calpha123_45)>1.0){printf("\nError F5 1.0 cos(calpha123_45)=%.20lf (%lf,%lf)\n",calpha123_45,mod123,mod45); }
	if(fabs(calpha1234_5)>1.0){printf("\nError F5 1.0 cos(calpha1234_5)=%.20lf (%lf,%lf) \n",calpha1234_5,mod1234,k5); }

	
	f1=(2.*5.+1.)*alpha(mod12345,k1,calpha12345_1)*F4(k2,k3,k4,k5,calpha23,calpha24,calpha25,calpha34,calpha35,calpha45)+2.*beta(mod12345,k1,mod2345,calpha1_2345)*G4(k2,k3,k4,k5,calpha23,calpha24,calpha25,calpha34,calpha35,calpha45);
		
	f2=G2(k1,k2,calpha12)*((2.*5.+1.)*alpha(mod12345,mod12,calpha12345_12)*F3(k3,k4,k5,calpha34,calpha35,calpha45)+2.*beta(mod12345,mod12,mod345,calpha12_345)*G3(k3,k4,k5,calpha34,calpha35,calpha45));
	
	f3=G3(k1,k2,k3,calpha12,calpha13,calpha23)*((2.*5.+1.)*alpha(mod12345,mod123,calpha12345_123)*F2(k4,k5,calpha45)+2.*beta(mod12345,mod123,mod45,calpha123_45)*G2(k4,k5,calpha45));
	
	f4=G4(k1,k2,k3,k4,calpha12,calpha13,calpha14,calpha23,calpha24,calpha34)*((2.*5.+1.)*alpha(mod12345,mod1234,calpha12345_1234)*1.+2.*beta(mod12345,mod1234,k5,calpha1234_5));
	
	
	f=pow((2.*5.+3)*(5.-1.),-1)*(f4+f3+f2+f1);
	

	return f;
	
}

double G5(double k1, double k2, double k3, double k4, double k5, double calpha12, double calpha13, double calpha14, double calpha15, double calpha23, double calpha24, double calpha25, double calpha34, double calpha35, double calpha45)
{
	double f,f1,f2,f3,f4;
	f1=0;
	f2=0;
	f3=0;
	f4=0;
	
	if (fabs(calpha12)>1 || fabs(calpha13)>1 || fabs(calpha14)>1 || fabs(calpha15)>1 || fabs(calpha23)>1 || fabs(calpha24)>1 || fabs(calpha25)>1 || fabs(calpha34)>1 || fabs(calpha35)>1 || fabs(calpha45)>1 ){printf("Error en la entrada de F5\n");}
	
	double epsilon=5e-2;
	double epsilon2=1e-4;
	double mod12345,mod12,mod123,mod1234, mod2345, mod345, mod45;
	double calpha12345_1,calpha12345_12,calpha12345_123,calpha12345_1234,calpha1_2345,calpha12_345,calpha123_45,calpha1234_5;
	
	mod12345=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + k5*k5 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k3*k4*calpha34 + k3*k5*calpha35 + k4*k5*calpha45)),0.5);
	
	mod12=pow(fabs(k1*k1 + k2*k2 + 2.*k1*k2*calpha12),0.5);
	
	mod123=pow(fabs(k1*k1 + k2*k2 + k3*k3 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k2*k3*calpha23)),0.5);
	
	mod1234=pow(fabs(k1*k1 + k2*k2 + k3*k3 + k4*k4 + 2.*(k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	
	mod2345=pow(fabs(k5*k5 + k2*k2 + k3*k3 + k4*k4 + 2.*(k5*k2*calpha25 + k5*k3*calpha35 + k5*k4*calpha45 + k2*k3*calpha23 + k2*k4*calpha24 + k3*k4*calpha34)),0.5);
	
	mod345=pow(fabs(k3*k3 + k4*k4 + k5*k5 + 2.*(k3*k4*calpha34 + k3*k5*calpha35 + k4*k5*calpha45)),0.5);
	
	mod45=pow(fabs(k4*k4 + k5*k5 + 2.*k4*k5*calpha45),0.5);
	
	calpha12345_1=( k1 + k2*calpha12 + k3*calpha13 + k4*calpha14 + k5*calpha15 )/(mod12345);
	
	calpha12345_12=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25)/(mod12345*mod12);
	
	calpha12345_123=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k1*k3*calpha13 + k2*k3*calpha23 + k3*k3 + k3*k4*calpha34 + k3*k5*calpha35 )/(mod12345*mod123);
	
	calpha12345_1234=( k1*k1 + k1*k2*calpha12 + k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k1*k2*calpha12 + k2*k2 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25 + k1*k3*calpha13 + k2*k3*calpha23 + k3*k3 + k3*k4*calpha34 + k3*k5*calpha35 + k4*k1*calpha14 + k4*k2*calpha24 + k4*k3*calpha34 + k4*k4 + k4*k5*calpha45)/(mod12345*mod1234);
	
	calpha1_2345=(k2*calpha12 + k3*calpha13 + k4*calpha14 + k5*calpha15)/mod2345;
	
	calpha12_345=(k1*k3*calpha13 + k1*k4*calpha14 + k1*k5*calpha15 + k2*k3*calpha23 + k2*k4*calpha24 + k2*k5*calpha25)/(mod12*mod345);
	
	calpha123_45=(k1*k4*calpha14 + k2*k4*calpha24 + k4*k3*calpha34 + k1*k5*calpha15 + k2*k5*calpha25 + k3*k5*calpha35)/(mod45*mod123);
	
	calpha1234_5=(k1*calpha15 + k2*calpha25 + k3*calpha35 + k4*calpha45)/(mod1234);
	
	if(mod12345==0){calpha12345_1=0; calpha12345_12=0; calpha12345_123=0;calpha12345_1234=0;}
	if(mod12==0){calpha12345_12=0;calpha12_345=0;}
	if(mod123==0){calpha12345_123=0;calpha123_45=0;}
	if(mod1234==0){calpha12345_1234=0;calpha1234_5=0;}
	if(mod2345==0){calpha1_2345=0;}
	if(mod345==0){calpha12_345=0;}
	if(mod45==0){calpha123_45=0;}
	if(k1==0){calpha12345_1=0;}
	if(k5==0){calpha1234_5=0;}
	
	
	//precision correction
	if(calpha12345_1>1.0 && calpha12345_1<1.0+epsilon){ calpha12345_1=1.0;}
	if(calpha12345_1<-1.0 && calpha12345_1>-1.0-epsilon){ calpha12345_1=-1.0;}
	
	if(calpha12345_12>1.0 && calpha12345_12<1.0+epsilon){ calpha12345_12=1.0;}
	if(calpha12345_12<-1.0 && calpha12345_12>-1.0-epsilon){ calpha12345_12=-1.0;}
	
	if(calpha12345_123>1.0 && calpha12345_123<1.0+epsilon){ calpha12345_123=1.0;}
	if(calpha12345_123<-1.0 && calpha12345_123>-1.0-epsilon){ calpha12345_123=-1.0;}
	
	if(calpha12345_1234>1.0 && calpha12345_1234<1.0+epsilon){ calpha12345_1234=1.0;}
	if(calpha12345_1234<-1.0 && calpha12345_1234>-1.0-epsilon){ calpha12345_1234=-1.0;}
	
	if(calpha1_2345>1.0 && calpha1_2345<1.0+epsilon){ calpha1_2345=1.0;}
	if(calpha1_2345<-1.0 && calpha1_2345>-1.0-epsilon){ calpha1_2345=-1.0;}
	
	if(calpha12_345>1.0 && calpha12_345<1.0+epsilon){ calpha12_345=1.0;}
	if(calpha12_345<-1.0 && calpha12_345>-1.0-epsilon){ calpha12_345=-1.0;}
	
	if(calpha123_45>1.0 && calpha123_45<1.0+epsilon){ calpha123_45=1.0;}
	if(calpha123_45<-1.0 && calpha123_45>-1.0-epsilon){ calpha123_45=-1.0;}
	
	if(calpha1234_5>1.0 && calpha1234_5<1.0+epsilon){ calpha1234_5=1.0;}
	if(calpha1234_5<-1.0 && calpha1234_5>-1.0-epsilon){ calpha1234_5=-1.0;}
	
	if(fabs(calpha12345_1)>1.0){printf("\nError G5 1.0 cos(calpha12345_1)=%.20lf (%lf,%lf) \n",calpha12345_1,mod12345,k1); }
	if(fabs(calpha12345_12)>1.0){printf("\nError G5 1.0 cos(calpha12345_12)=%.20lf (%lf,%lf) \n",calpha12345_12,mod12345,mod12); }
	if(fabs(calpha12345_123)>1.0){printf("\nError G5 1.0 cos(calpha12345_123)=%.20lf (%lf,%lf) \n",calpha12345_123,mod12345,mod123); }
	if(fabs(calpha12345_1234)>1.0){printf("\nError G5 1.0 cos(calpha12345_1234)=%.20lf (%lf,%lf) \n",calpha12345_1234,mod12345,mod1234); }
	if(fabs(calpha1_2345)>1.0){printf("\nError G5 1.0 cos(calpha1_2345)=%.20lf (%lf,%lf) \n",calpha1_2345,k1,mod2345); }
	if(fabs(calpha12_345)>1.0){printf("\nError G5 1.0 cos(calpha12_345)=%.20lf (%lf,%lf) \n",calpha12_345,mod12,mod345); }
	if(fabs(calpha123_45)>1.0){printf("\nError G5 1.0 cos(calpha123_45)=%.20lf (%lf,%lf)\n",calpha123_45,mod123,mod45); }
	if(fabs(calpha1234_5)>1.0){printf("\nError G5 1.0 cos(calpha1234_5)=%.20lf (%lf,%lf) \n",calpha1234_5,mod1234,k5); }
	
	
	f1=(3.)*alpha(mod12345,k1,calpha12345_1)*F4(k2,k3,k4,k5,calpha23,calpha24,calpha25,calpha34,calpha35,calpha45)+10.*beta(mod12345,k1,mod2345,calpha1_2345)*G4(k2,k3,k4,k5,calpha23,calpha24,calpha25,calpha34,calpha35,calpha45);
	
	f2=G2(k1,k2,calpha12)*((3.)*alpha(mod12345,mod12,calpha12345_12)*F3(k3,k4,k5,calpha34,calpha35,calpha45)+10.*beta(mod12345,mod12,mod345,calpha12_345)*G3(k3,k4,k5,calpha34,calpha35,calpha45));
	
	f3=G3(k1,k2,k3,calpha12,calpha13,calpha23)*((3.)*alpha(mod12345,mod123,calpha12345_123)*F2(k4,k5,calpha45)+10.*beta(mod12345,mod123,mod45,calpha123_45)*G2(k4,k5,calpha45));
	
	f4=G4(k1,k2,k3,k4,calpha12,calpha13,calpha14,calpha23,calpha24,calpha34)*((3.)*alpha(mod12345,mod1234,calpha12345_1234)*1.+10.*beta(mod12345,mod1234,k5,calpha1234_5));
	
	
	f=pow((2.*5.+3)*(5.-1.),-1)*(f4+f3+f2+f1);
	
	
	return f;
	
}



//symetryzed-2-points kernel (for matter fields)
double Fs2(double k1, double k2, double calpha12) 
{
	
	double f;
	f=(1./2.)*(F2(k1,k2,calpha12)+F2(k2,k1,calpha12));
	return f;
	 
	
}
double Gs2(double k1, double k2, double calpha12)
{
	
double f;
f=(1./2.)*(G2(k1,k2,calpha12)+G2(k2,k1,calpha12));
return f;
	 
}



//symetryzed-3-points kernel (for matter fields)
double Fs3(double k1, double k2, double k3, double calpha12, double calpha13, double calpha23 )
{
	
	double f;
        f=(1./(1.*2.*3.))*( F3(k1, k2, k3, calpha12, calpha13, calpha23) + F3(k3, k1, k2, calpha13, calpha23, calpha12) + F3(k2,k3,k1, calpha23, calpha12, calpha13) + F3(k1,k3,k2, calpha13, calpha12, calpha23) + F3(k3,k2,k1,calpha23,calpha13,calpha12) + F3(k2,k1,k3,calpha12,calpha23,calpha13) ); 

//3!=6 terms
//123, 312, 231, 132, 321, 213

	return f;
	
}

double Gs3(double k1, double k2, double k3, double calpha12, double calpha13, double calpha23 )
{
	double f;
        f=(1./(1.*2.*3.))*( G3(k1, k2, k3, calpha12, calpha13, calpha23) + G3(k3, k1, k2, calpha13, calpha23, calpha12) + G3(k2,k3,k1, calpha23, calpha12, calpha13) + G3(k1,k3,k2, calpha13, calpha12, calpha23) + G3(k3,k2,k1,calpha23,calpha13,calpha12) + G3(k2,k1,k3,calpha12,calpha23,calpha13) ); 

//3!=6 terms
//123, 312, 231, 132, 321, 213

	return f;
}

double Fs4(double k1, double k2, double k3, double k4, double calpha12, double calpha13, double calpha14, double calpha23, double calpha24, double calpha34)
{
	double f;
	
	f=(1./(1.*2.*3.*4.))*(F4(k1,k4,k2,k3,calpha14,calpha12,calpha13,calpha24,calpha34,calpha23)+F4(k1,k3,k4,k2,calpha13,calpha14,calpha12,calpha34,calpha23,calpha24)+F4(k1,k2,k3,k4,calpha12,calpha13,calpha14,calpha23,calpha24,calpha34)+F4(k1,k4,k3,k2,calpha14,calpha13,calpha12,calpha34,calpha24,calpha23)+F4(k1,k3,k2,k4,calpha13,calpha12,calpha14,calpha23,calpha34,calpha24)+F4(k1,k2,k4,k3,calpha12,calpha14,calpha13,calpha24,calpha23,calpha34)+F4(k2,k1,k4,k3,calpha12,calpha24,calpha23,calpha14,calpha13,calpha34)+F4(k2,k3,k1,k4,calpha23,calpha12,calpha24,calpha13,calpha34,calpha14)+F4(k2,k4,k3,k1,calpha24,calpha23,calpha12,calpha34,calpha14,calpha13)+F4(k2,k1,k3,k4,calpha12,calpha23,calpha24,calpha13,calpha14,calpha34)+F4(k2,k3,k4,k1,calpha23,calpha24,calpha12,calpha34,calpha13,calpha14)+F4(k2,k4,k1,k3,calpha24,calpha12,calpha23,calpha14,calpha34,calpha13)+F4(k3,k1,k2,k4,calpha13,calpha23,calpha34,calpha12,calpha14,calpha24)+F4(k3,k4,k1,k2,calpha34,calpha13,calpha23,calpha14,calpha24,calpha12)+F4(k3,k2,k4,k1,calpha23,calpha34,calpha13,calpha24,calpha12,calpha14)+F4(k3,k1,k4,k2,calpha13,calpha34,calpha23,calpha14,calpha12,calpha24)+F4(k3,k4,k2,k1,calpha34,calpha23,calpha13,calpha24,calpha14,calpha12)+F4(k3,k2,k1,k4,calpha23,calpha13,calpha34,calpha12,calpha24,calpha14)+F4(k4,k1,k2,k3,calpha14,calpha24,calpha34,calpha12,calpha13,calpha23)+F4(k4,k3,k1,k2,calpha34,calpha14,calpha24,calpha13,calpha23,calpha12)+F4(k4,k2,k3,k1,calpha24,calpha34,calpha14,calpha23,calpha12,calpha13)+F4(k4,k1,k3,k2,calpha14,calpha34,calpha24,calpha13,calpha12,calpha23)+F4(k4,k3,k2,k1,calpha34,calpha24,calpha14,calpha23,calpha13,calpha12)+F4(k4,k2,k1,k3,calpha24,calpha14,calpha34,calpha12,calpha23,calpha13));
	
	//4!=24 terms
	
	//1423, 1342, 1234, 1432, 1324, 1243, 2143, 2314, 2431, 2134, 2341, 2413, 3124, 3412, 3241, 3142, 3421, 3214, 4123, 4312, 4231, 4132, 4321, 4213

	return f;
}

double Gs4(double k1, double k2, double k3, double k4, double calpha12, double calpha13, double calpha14, double calpha23, double calpha24, double calpha34)
{
	double f;
	
	f=(1./(1.*2.*3.*4.))*(G4(k1,k4,k2,k3,calpha14,calpha12,calpha13,calpha24,calpha34,calpha23)+G4(k1,k3,k4,k2,calpha13,calpha14,calpha12,calpha34,calpha23,calpha24)+G4(k1,k2,k3,k4,calpha12,calpha13,calpha14,calpha23,calpha24,calpha34)+G4(k1,k4,k3,k2,calpha14,calpha13,calpha12,calpha34,calpha24,calpha23)+G4(k1,k3,k2,k4,calpha13,calpha12,calpha14,calpha23,calpha34,calpha24)+G4(k1,k2,k4,k3,calpha12,calpha14,calpha13,calpha24,calpha23,calpha34)+G4(k2,k1,k4,k3,calpha12,calpha24,calpha23,calpha14,calpha13,calpha34)+G4(k2,k3,k1,k4,calpha23,calpha12,calpha24,calpha13,calpha34,calpha14)+G4(k2,k4,k3,k1,calpha24,calpha23,calpha12,calpha34,calpha14,calpha13)+G4(k2,k1,k3,k4,calpha12,calpha23,calpha24,calpha13,calpha14,calpha34)+G4(k2,k3,k4,k1,calpha23,calpha24,calpha12,calpha34,calpha13,calpha14)+G4(k2,k4,k1,k3,calpha24,calpha12,calpha23,calpha14,calpha34,calpha13)+G4(k3,k1,k2,k4,calpha13,calpha23,calpha34,calpha12,calpha14,calpha24)+G4(k3,k4,k1,k2,calpha34,calpha13,calpha23,calpha14,calpha24,calpha12)+G4(k3,k2,k4,k1,calpha23,calpha34,calpha13,calpha24,calpha12,calpha14)+G4(k3,k1,k4,k2,calpha13,calpha34,calpha23,calpha14,calpha12,calpha24)+G4(k3,k4,k2,k1,calpha34,calpha23,calpha13,calpha24,calpha14,calpha12)+G4(k3,k2,k1,k4,calpha23,calpha13,calpha34,calpha12,calpha24,calpha14)+G4(k4,k1,k2,k3,calpha14,calpha24,calpha34,calpha12,calpha13,calpha23)+G4(k4,k3,k1,k2,calpha34,calpha14,calpha24,calpha13,calpha23,calpha12)+G4(k4,k2,k3,k1,calpha24,calpha34,calpha14,calpha23,calpha12,calpha13)+G4(k4,k1,k3,k2,calpha14,calpha34,calpha24,calpha13,calpha12,calpha23)+G4(k4,k3,k2,k1,calpha34,calpha24,calpha14,calpha23,calpha13,calpha12)+G4(k4,k2,k1,k3,calpha24,calpha14,calpha34,calpha12,calpha23,calpha13));
	
	//4!=24 terms
	
	//1423, 1342, 1234, 1432, 1324, 1243, 2143, 2314, 2431, 2134, 2341, 2413, 3124, 3412, 3241, 3142, 3421, 3214, 4123, 4312, 4231, 4132, 4321, 4213
	
	return f;
}

double Fs5(double k1, double k2, double k3, double k4, double k5, double calpha12, double calpha13, double calpha14, double calpha15, double calpha23, double calpha24, double calpha25, double calpha34, double calpha35, double calpha45)
{
	double f=0;
	double k0[6];
	
	double theta[6][6];
	
	int a1,a2,a3,a4,a5;//asignadores. son 0 si no han sido asignados, son 1 si han sido asignados
	int inicio_j,inicio_k,inicio_l,inicio_m;
	
	int i,j,k,l,m;
	for(i=0;i<5;i++)
	{
		a1=0;
		a2=0;
		a3=0;
		a4=0;
		a5=0;
		if(i==0){k0[1]=k1;a1=1;}
		if(i==1){k0[1]=k2;a2=1;}
		if(i==2){k0[1]=k3;a3=1;}
		if(i==3){k0[1]=k4;a4=1;}
		if(i==4){k0[1]=k5;a5=1;}
		
		for (j=0; j<4; j++)
		{
			if (j==0){inicio_j=1;}
			
			//asignar a 0 todos excepto los del loop exterior
			if(a1>=2){a1=0;}
			if(a2>=2){a2=0;}
			if(a3>=2){a3=0;}
			if(a4>=2){a4=0;}
			if(a5>=2){a5=0;}
			
			if(j==0){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min1 da preferencia a k1 si no a sido a signado, luego a k2, etc...
			if(j==1){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min2 da preferencia a k2 si no a sido a signado, luego a k3, etc...
			if(j==2){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min3 da preferencia a k3 si no a sido a signado, luego a k4, etc...
			if(j==3){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min4 da preferencia a k4 si no a sido a signado, luego a k1, etc...
			
			if (a1==1 && a2==2){theta[1][2]=calpha12;}
			if (a1==1 && a3==2){theta[1][2]=calpha13;}
			if (a1==1 && a4==2){theta[1][2]=calpha14;}
			if (a1==1 && a5==2){theta[1][2]=calpha15;}
			if (a2==1 && a3==2){theta[1][2]=calpha23;}
			if (a2==1 && a4==2){theta[1][2]=calpha24;}
			if (a2==1 && a5==2){theta[1][2]=calpha25;}
			if (a3==1 && a4==2){theta[1][2]=calpha34;}
			if (a3==1 && a5==2){theta[1][2]=calpha35;}
			if (a4==1 && a5==2){theta[1][2]=calpha45;}
			
			if (a1==2 && a2==1){theta[1][2]=calpha12;}
			if (a1==2 && a3==1){theta[1][2]=calpha13;}
			if (a1==2 && a4==1){theta[1][2]=calpha14;}
			if (a1==2 && a5==1){theta[1][2]=calpha15;}
			if (a2==2 && a3==1){theta[1][2]=calpha23;}
			if (a2==2 && a4==1){theta[1][2]=calpha24;}
			if (a2==2 && a5==1){theta[1][2]=calpha25;}
			if (a3==2 && a4==1){theta[1][2]=calpha34;}
			if (a3==2 && a5==1){theta[1][2]=calpha35;}
			if (a4==2 && a5==1){theta[1][2]=calpha45;}
			
            for (k=0; k<3; k++)
			{
				if(k==0){inicio_k=1;}
				
				if(a1>=3){a1=0;}
				if(a2>=3){a2=0;}
				if(a3>=3){a3=0;}
				if(a4>=3){a4=0;}
				if(a5>=3){a5=0;}
				
				if(k==0){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				if(k==1){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				if(k==2){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				
				if (a1==1 && a2==3){theta[1][3]=calpha12;}
				if (a1==1 && a3==3){theta[1][3]=calpha13;}
				if (a1==1 && a4==3){theta[1][3]=calpha14;}
				if (a1==1 && a5==3){theta[1][3]=calpha15;}
				if (a2==1 && a3==3){theta[1][3]=calpha23;}
				if (a2==1 && a4==3){theta[1][3]=calpha24;}
				if (a2==1 && a5==3){theta[1][3]=calpha25;}
				if (a3==1 && a4==3){theta[1][3]=calpha34;}
				if (a3==1 && a5==3){theta[1][3]=calpha35;}
				if (a4==1 && a5==3){theta[1][3]=calpha45;}
				
				if (a1==3 && a2==1){theta[1][3]=calpha12;}
				if (a1==3 && a3==1){theta[1][3]=calpha13;}
				if (a1==3 && a4==1){theta[1][3]=calpha14;}
				if (a1==3 && a5==1){theta[1][3]=calpha15;}
				if (a2==3 && a3==1){theta[1][3]=calpha23;}
				if (a2==3 && a4==1){theta[1][3]=calpha24;}
				if (a2==3 && a5==1){theta[1][3]=calpha25;}
				if (a3==3 && a4==1){theta[1][3]=calpha34;}
				if (a3==3 && a5==1){theta[1][3]=calpha35;}
				if (a4==3 && a5==1){theta[1][3]=calpha45;}
				
				if (a1==2 && a2==3){theta[2][3]=calpha12;}
				if (a1==2 && a3==3){theta[2][3]=calpha13;}
				if (a1==2 && a4==3){theta[2][3]=calpha14;}
				if (a1==2 && a5==3){theta[2][3]=calpha15;}
				if (a2==2 && a3==3){theta[2][3]=calpha23;}
				if (a2==2 && a4==3){theta[2][3]=calpha24;}
				if (a2==2 && a5==3){theta[2][3]=calpha25;}
				if (a3==2 && a4==3){theta[2][3]=calpha34;}
				if (a3==2 && a5==3){theta[2][3]=calpha35;}
				if (a4==2 && a5==3){theta[2][3]=calpha45;}
				
				if (a1==3 && a2==2){theta[2][3]=calpha12;}
				if (a1==3 && a3==2){theta[2][3]=calpha13;}
				if (a1==3 && a4==2){theta[2][3]=calpha14;}
				if (a1==3 && a5==2){theta[2][3]=calpha15;}
				if (a2==3 && a3==2){theta[2][3]=calpha23;}
				if (a2==3 && a4==2){theta[2][3]=calpha24;}
				if (a2==3 && a5==2){theta[2][3]=calpha25;}
				if (a3==3 && a4==2){theta[2][3]=calpha34;}
				if (a3==3 && a5==2){theta[2][3]=calpha35;}
				if (a4==3 && a5==2){theta[2][3]=calpha45;}
				
				for (l=0; l<2; l++)
				{
					if (l==0) {inicio_l=1;}
					
					if(a1>=4){a1=0;}
					if(a2>=4){a2=0;}
					if(a3>=4){a3=0;}
					if(a4>=4){a4=0;}
					if(a5>=4){a5=0;}
					
					if(l==0){k0[4]=min(&inicio_l,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,4);}
					if(l==1){k0[4]=min(&inicio_l,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,4);}
					
					if (a1==1 && a2==4){theta[1][4]=calpha12;}
					if (a1==1 && a3==4){theta[1][4]=calpha13;}
					if (a1==1 && a4==4){theta[1][4]=calpha14;}
					if (a1==1 && a5==4){theta[1][4]=calpha15;}
					if (a2==1 && a3==4){theta[1][4]=calpha23;}
					if (a2==1 && a4==4){theta[1][4]=calpha24;}
					if (a2==1 && a5==4){theta[1][4]=calpha25;}
					if (a3==1 && a4==4){theta[1][4]=calpha34;}
					if (a3==1 && a5==4){theta[1][4]=calpha35;}
					if (a4==1 && a5==4){theta[1][4]=calpha45;}
					
					if (a1==4 && a2==1){theta[1][4]=calpha12;}
					if (a1==4 && a3==1){theta[1][4]=calpha13;}
					if (a1==4 && a4==1){theta[1][4]=calpha14;}
					if (a1==4 && a5==1){theta[1][4]=calpha15;}
					if (a2==4 && a3==1){theta[1][4]=calpha23;}
					if (a2==4 && a4==1){theta[1][4]=calpha24;}
					if (a2==4 && a5==1){theta[1][4]=calpha25;}
					if (a3==4 && a4==1){theta[1][4]=calpha34;}
					if (a3==4 && a5==1){theta[1][4]=calpha35;}
					if (a4==4 && a5==1){theta[1][4]=calpha45;}
				   
					if (a1==2 && a2==4){theta[2][4]=calpha12;}
					if (a1==2 && a3==4){theta[2][4]=calpha13;}
					if (a1==2 && a4==4){theta[2][4]=calpha14;}
					if (a1==2 && a5==4){theta[2][4]=calpha15;}
					if (a2==2 && a3==4){theta[2][4]=calpha23;}
					if (a2==2 && a4==4){theta[2][4]=calpha24;}
					if (a2==2 && a5==4){theta[2][4]=calpha25;}
					if (a3==2 && a4==4){theta[2][4]=calpha34;}
					if (a3==2 && a5==4){theta[2][4]=calpha35;}
					if (a4==2 && a5==4){theta[2][4]=calpha45;}
					
					if (a1==4 && a2==2){theta[2][4]=calpha12;}
					if (a1==4 && a3==2){theta[2][4]=calpha13;}
					if (a1==4 && a4==2){theta[2][4]=calpha14;}
					if (a1==4 && a5==2){theta[2][4]=calpha15;}
					if (a2==4 && a3==2){theta[2][4]=calpha23;}
					if (a2==4 && a4==2){theta[2][4]=calpha24;}
					if (a2==4 && a5==2){theta[2][4]=calpha25;}
					if (a3==4 && a4==2){theta[2][4]=calpha34;}
					if (a3==4 && a5==2){theta[2][4]=calpha35;}
					if (a4==4 && a5==2){theta[2][4]=calpha45;}
					
					if (a1==3 && a2==4){theta[3][4]=calpha12;}
					if (a1==3 && a3==4){theta[3][4]=calpha13;}
					if (a1==3 && a4==4){theta[3][4]=calpha14;}
					if (a1==3 && a5==4){theta[3][4]=calpha15;}
					if (a2==3 && a3==4){theta[3][4]=calpha23;}
					if (a2==3 && a4==4){theta[3][4]=calpha24;}
					if (a2==3 && a5==4){theta[3][4]=calpha25;}
					if (a3==3 && a4==4){theta[3][4]=calpha34;}
					if (a3==3 && a5==4){theta[3][4]=calpha35;}
					if (a4==3 && a5==4){theta[3][4]=calpha45;}
					
					if (a1==4 && a2==3){theta[3][4]=calpha12;}
					if (a1==4 && a3==3){theta[3][4]=calpha13;}
					if (a1==4 && a4==3){theta[3][4]=calpha14;}
					if (a1==4 && a5==3){theta[3][4]=calpha15;}
					if (a2==4 && a3==3){theta[3][4]=calpha23;}
					if (a2==4 && a4==3){theta[3][4]=calpha24;}
					if (a2==4 && a5==3){theta[3][4]=calpha25;}
					if (a3==4 && a4==3){theta[3][4]=calpha34;}
					if (a3==4 && a5==3){theta[3][4]=calpha35;}
					if (a4==4 && a5==3){theta[3][4]=calpha45;}
					
					
					for (m=0; m<1; m++)
					{
						if(m==0){inicio_m=1;}
						if(a1>=5){a1=0;}
						if(a2>=5){a2=0;}
						if(a3>=5){a3=0;}
						if(a4>=5){a4=0;}
						if(a5>=5){a5=0;}
						
				      	if(m==0){k0[5]=min(&inicio_m,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,5);}
						
						if (a1==1 && a2==5){theta[1][5]=calpha12;}
						if (a1==1 && a3==5){theta[1][5]=calpha13;}
						if (a1==1 && a4==5){theta[1][5]=calpha14;}
						if (a1==1 && a5==5){theta[1][5]=calpha15;}
						if (a2==1 && a3==5){theta[1][5]=calpha23;}
						if (a2==1 && a4==5){theta[1][5]=calpha24;}
						if (a2==1 && a5==5){theta[1][5]=calpha25;}
						if (a3==1 && a4==5){theta[1][5]=calpha34;}
						if (a3==1 && a5==5){theta[1][5]=calpha35;}
						if (a4==1 && a5==5){theta[1][5]=calpha45;}
						
						if (a1==5 && a2==1){theta[1][5]=calpha12;}
						if (a1==5 && a3==1){theta[1][5]=calpha13;}
						if (a1==5 && a4==1){theta[1][5]=calpha14;}
						if (a1==5 && a5==1){theta[1][5]=calpha15;}
						if (a2==5 && a3==1){theta[1][5]=calpha23;}
						if (a2==5 && a4==1){theta[1][5]=calpha24;}
						if (a2==5 && a5==1){theta[1][5]=calpha25;}
						if (a3==5 && a4==1){theta[1][5]=calpha34;}
						if (a3==5 && a5==1){theta[1][5]=calpha35;}
						if (a4==5 && a5==1){theta[1][5]=calpha45;}
						
						if (a1==2 && a2==5){theta[2][5]=calpha12;}
						if (a1==2 && a3==5){theta[2][5]=calpha13;}
						if (a1==2 && a4==5){theta[2][5]=calpha14;}
						if (a1==2 && a5==5){theta[2][5]=calpha15;}
						if (a2==2 && a3==5){theta[2][5]=calpha23;}
						if (a2==2 && a4==5){theta[2][5]=calpha24;}
						if (a2==2 && a5==5){theta[2][5]=calpha25;}
						if (a3==2 && a4==5){theta[2][5]=calpha34;}
						if (a3==2 && a5==5){theta[2][5]=calpha35;}
						if (a4==2 && a5==5){theta[2][5]=calpha45;}
						
						if (a1==5 && a2==2){theta[2][5]=calpha12;}
						if (a1==5 && a3==2){theta[2][5]=calpha13;}
						if (a1==5 && a4==2){theta[2][5]=calpha14;}
						if (a1==5 && a5==2){theta[2][5]=calpha15;}
						if (a2==5 && a3==2){theta[2][5]=calpha23;}
						if (a2==5 && a4==2){theta[2][5]=calpha24;}
						if (a2==5 && a5==2){theta[2][5]=calpha25;}
						if (a3==5 && a4==2){theta[2][5]=calpha34;}
						if (a3==5 && a5==2){theta[2][5]=calpha35;}
						if (a4==5 && a5==2){theta[2][5]=calpha45;}
						
						if (a1==3 && a2==5){theta[3][5]=calpha12;}
						if (a1==3 && a3==5){theta[3][5]=calpha13;}
						if (a1==3 && a4==5){theta[3][5]=calpha14;}
						if (a1==3 && a5==5){theta[3][5]=calpha15;}
						if (a2==3 && a3==5){theta[3][5]=calpha23;}
						if (a2==3 && a4==5){theta[3][5]=calpha24;}
						if (a2==3 && a5==5){theta[3][5]=calpha25;}
						if (a3==3 && a4==5){theta[3][5]=calpha34;}
						if (a3==3 && a5==5){theta[3][5]=calpha35;}
						if (a4==3 && a5==5){theta[3][5]=calpha45;}
						
						if (a1==5 && a2==3){theta[3][5]=calpha12;}
						if (a1==5 && a3==3){theta[3][5]=calpha13;}
						if (a1==5 && a4==3){theta[3][5]=calpha14;}
						if (a1==5 && a5==3){theta[3][5]=calpha15;}
						if (a2==5 && a3==3){theta[3][5]=calpha23;}
						if (a2==5 && a4==3){theta[3][5]=calpha24;}
						if (a2==5 && a5==3){theta[3][5]=calpha25;}
						if (a3==5 && a4==3){theta[3][5]=calpha34;}
						if (a3==5 && a5==3){theta[3][5]=calpha35;}
						if (a4==5 && a5==3){theta[3][5]=calpha45;}
						
						if (a1==4 && a2==5){theta[4][5]=calpha12;}
						if (a1==4 && a3==5){theta[4][5]=calpha13;}
						if (a1==4 && a4==5){theta[4][5]=calpha14;}
						if (a1==4 && a5==5){theta[4][5]=calpha15;}
						if (a2==4 && a3==5){theta[4][5]=calpha23;}
						if (a2==4 && a4==5){theta[4][5]=calpha24;}
						if (a2==4 && a5==5){theta[4][5]=calpha25;}
						if (a3==4 && a4==5){theta[4][5]=calpha34;}
						if (a3==4 && a5==5){theta[4][5]=calpha35;}
						if (a4==4 && a5==5){theta[4][5]=calpha45;}
						
						if (a1==5 && a2==4){theta[4][5]=calpha12;}
						if (a1==5 && a3==4){theta[4][5]=calpha13;}
						if (a1==5 && a4==4){theta[4][5]=calpha14;}
						if (a1==5 && a5==4){theta[4][5]=calpha15;}
						if (a2==5 && a3==4){theta[4][5]=calpha23;}
						if (a2==5 && a4==4){theta[4][5]=calpha24;}
						if (a2==5 && a5==4){theta[4][5]=calpha25;}
						if (a3==5 && a4==4){theta[4][5]=calpha34;}
						if (a3==5 && a5==4){theta[4][5]=calpha35;}
						if (a4==5 && a5==4){theta[4][5]=calpha45;}
						
			            f=f+F5(k0[1],k0[2],k0[3],k0[4],k0[5],theta[1][2],theta[1][3],theta[1][4],theta[1][5],theta[2][3],theta[2][4],theta[2][5],theta[3][4],theta[3][5],theta[4][5]);
				//	printf("(%lf %lf %lf %lf %lf, %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf (%d,%d,%d,%d,%d) \n",k0[1],k0[2],k0[3],k0[4],k0[5],theta[1][2],theta[1][3],theta[1][4],theta[1][5],theta[2][3],theta[2][4],theta[2][5],theta[3][4],theta[3][5],theta[4][5],a1,a2,a3,a4,a5);
					}
				}
			}
		}
	}
	
	f=f/120.; 
return f;
}

double Gs5(double k1, double k2, double k3, double k4, double k5, double calpha12, double calpha13, double calpha14, double calpha15, double calpha23, double calpha24, double calpha25, double calpha34, double calpha35, double calpha45)
{
	double f=0;
	double k0[6];
	
	double theta[6][6];
	
	int a1,a2,a3,a4,a5;//asignadores. son 0 si no han sido asignados, son 1 si han sido asignados
	int inicio_j,inicio_k,inicio_l,inicio_m;
	
	int i,j,k,l,m;
	for(i=0;i<5;i++)
	{
		a1=0;
		a2=0;
		a3=0;
		a4=0;
		a5=0;
		if(i==0){k0[1]=k1;a1=1;}
		if(i==1){k0[1]=k2;a2=1;}
		if(i==2){k0[1]=k3;a3=1;}
		if(i==3){k0[1]=k4;a4=1;}
		if(i==4){k0[1]=k5;a5=1;}
		
		for (j=0; j<4; j++)
		{
			if (j==0){inicio_j=1;}
			
			//asignar a 0 todos excepto los del loop exterior
			if(a1>=2){a1=0;}
			if(a2>=2){a2=0;}
			if(a3>=2){a3=0;}
			if(a4>=2){a4=0;}
			if(a5>=2){a5=0;}
			
			if(j==0){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min1 da preferencia a k1 si no a sido a signado, luego a k2, etc...
			if(j==1){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min2 da preferencia a k2 si no a sido a signado, luego a k3, etc...
			if(j==2){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min3 da preferencia a k3 si no a sido a signado, luego a k4, etc...
			if(j==3){k0[2]=min(&inicio_j,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,2);}//min4 da preferencia a k4 si no a sido a signado, luego a k1, etc...
			
			if (a1==1 && a2==2){theta[1][2]=calpha12;}
			if (a1==1 && a3==2){theta[1][2]=calpha13;}
			if (a1==1 && a4==2){theta[1][2]=calpha14;}
			if (a1==1 && a5==2){theta[1][2]=calpha15;}
			if (a2==1 && a3==2){theta[1][2]=calpha23;}
			if (a2==1 && a4==2){theta[1][2]=calpha24;}
			if (a2==1 && a5==2){theta[1][2]=calpha25;}
			if (a3==1 && a4==2){theta[1][2]=calpha34;}
			if (a3==1 && a5==2){theta[1][2]=calpha35;}
			if (a4==1 && a5==2){theta[1][2]=calpha45;}
			
			if (a1==2 && a2==1){theta[1][2]=calpha12;}
			if (a1==2 && a3==1){theta[1][2]=calpha13;}
			if (a1==2 && a4==1){theta[1][2]=calpha14;}
			if (a1==2 && a5==1){theta[1][2]=calpha15;}
			if (a2==2 && a3==1){theta[1][2]=calpha23;}
			if (a2==2 && a4==1){theta[1][2]=calpha24;}
			if (a2==2 && a5==1){theta[1][2]=calpha25;}
			if (a3==2 && a4==1){theta[1][2]=calpha34;}
			if (a3==2 && a5==1){theta[1][2]=calpha35;}
			if (a4==2 && a5==1){theta[1][2]=calpha45;}
			
            for (k=0; k<3; k++)
			{
				if(k==0){inicio_k=1;}
				
				if(a1>=3){a1=0;}
				if(a2>=3){a2=0;}
				if(a3>=3){a3=0;}
				if(a4>=3){a4=0;}
				if(a5>=3){a5=0;}
				
				if(k==0){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				if(k==1){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				if(k==2){k0[3]=min(&inicio_k,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,3);}
				
				if (a1==1 && a2==3){theta[1][3]=calpha12;}
				if (a1==1 && a3==3){theta[1][3]=calpha13;}
				if (a1==1 && a4==3){theta[1][3]=calpha14;}
				if (a1==1 && a5==3){theta[1][3]=calpha15;}
				if (a2==1 && a3==3){theta[1][3]=calpha23;}
				if (a2==1 && a4==3){theta[1][3]=calpha24;}
				if (a2==1 && a5==3){theta[1][3]=calpha25;}
				if (a3==1 && a4==3){theta[1][3]=calpha34;}
				if (a3==1 && a5==3){theta[1][3]=calpha35;}
				if (a4==1 && a5==3){theta[1][3]=calpha45;}
				
				if (a1==3 && a2==1){theta[1][3]=calpha12;}
				if (a1==3 && a3==1){theta[1][3]=calpha13;}
				if (a1==3 && a4==1){theta[1][3]=calpha14;}
				if (a1==3 && a5==1){theta[1][3]=calpha15;}
				if (a2==3 && a3==1){theta[1][3]=calpha23;}
				if (a2==3 && a4==1){theta[1][3]=calpha24;}
				if (a2==3 && a5==1){theta[1][3]=calpha25;}
				if (a3==3 && a4==1){theta[1][3]=calpha34;}
				if (a3==3 && a5==1){theta[1][3]=calpha35;}
				if (a4==3 && a5==1){theta[1][3]=calpha45;}
				
				if (a1==2 && a2==3){theta[2][3]=calpha12;}
				if (a1==2 && a3==3){theta[2][3]=calpha13;}
				if (a1==2 && a4==3){theta[2][3]=calpha14;}
				if (a1==2 && a5==3){theta[2][3]=calpha15;}
				if (a2==2 && a3==3){theta[2][3]=calpha23;}
				if (a2==2 && a4==3){theta[2][3]=calpha24;}
				if (a2==2 && a5==3){theta[2][3]=calpha25;}
				if (a3==2 && a4==3){theta[2][3]=calpha34;}
				if (a3==2 && a5==3){theta[2][3]=calpha35;}
				if (a4==2 && a5==3){theta[2][3]=calpha45;}
				
				if (a1==3 && a2==2){theta[2][3]=calpha12;}
				if (a1==3 && a3==2){theta[2][3]=calpha13;}
				if (a1==3 && a4==2){theta[2][3]=calpha14;}
				if (a1==3 && a5==2){theta[2][3]=calpha15;}
				if (a2==3 && a3==2){theta[2][3]=calpha23;}
				if (a2==3 && a4==2){theta[2][3]=calpha24;}
				if (a2==3 && a5==2){theta[2][3]=calpha25;}
				if (a3==3 && a4==2){theta[2][3]=calpha34;}
				if (a3==3 && a5==2){theta[2][3]=calpha35;}
				if (a4==3 && a5==2){theta[2][3]=calpha45;}
				
				for (l=0; l<2; l++)
				{
					if (l==0) {inicio_l=1;}
					
					if(a1>=4){a1=0;}
					if(a2>=4){a2=0;}
					if(a3>=4){a3=0;}
					if(a4>=4){a4=0;}
					if(a5>=4){a5=0;}
					
					if(l==0){k0[4]=min(&inicio_l,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,4);}
					if(l==1){k0[4]=min(&inicio_l,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,4);}
					
					if (a1==1 && a2==4){theta[1][4]=calpha12;}
					if (a1==1 && a3==4){theta[1][4]=calpha13;}
					if (a1==1 && a4==4){theta[1][4]=calpha14;}
					if (a1==1 && a5==4){theta[1][4]=calpha15;}
					if (a2==1 && a3==4){theta[1][4]=calpha23;}
					if (a2==1 && a4==4){theta[1][4]=calpha24;}
					if (a2==1 && a5==4){theta[1][4]=calpha25;}
					if (a3==1 && a4==4){theta[1][4]=calpha34;}
					if (a3==1 && a5==4){theta[1][4]=calpha35;}
					if (a4==1 && a5==4){theta[1][4]=calpha45;}
					
					if (a1==4 && a2==1){theta[1][4]=calpha12;}
					if (a1==4 && a3==1){theta[1][4]=calpha13;}
					if (a1==4 && a4==1){theta[1][4]=calpha14;}
					if (a1==4 && a5==1){theta[1][4]=calpha15;}
					if (a2==4 && a3==1){theta[1][4]=calpha23;}
					if (a2==4 && a4==1){theta[1][4]=calpha24;}
					if (a2==4 && a5==1){theta[1][4]=calpha25;}
					if (a3==4 && a4==1){theta[1][4]=calpha34;}
					if (a3==4 && a5==1){theta[1][4]=calpha35;}
					if (a4==4 && a5==1){theta[1][4]=calpha45;}
					
					if (a1==2 && a2==4){theta[2][4]=calpha12;}
					if (a1==2 && a3==4){theta[2][4]=calpha13;}
					if (a1==2 && a4==4){theta[2][4]=calpha14;}
					if (a1==2 && a5==4){theta[2][4]=calpha15;}
					if (a2==2 && a3==4){theta[2][4]=calpha23;}
					if (a2==2 && a4==4){theta[2][4]=calpha24;}
					if (a2==2 && a5==4){theta[2][4]=calpha25;}
					if (a3==2 && a4==4){theta[2][4]=calpha34;}
					if (a3==2 && a5==4){theta[2][4]=calpha35;}
					if (a4==2 && a5==4){theta[2][4]=calpha45;}
					
					if (a1==4 && a2==2){theta[2][4]=calpha12;}
					if (a1==4 && a3==2){theta[2][4]=calpha13;}
					if (a1==4 && a4==2){theta[2][4]=calpha14;}
					if (a1==4 && a5==2){theta[2][4]=calpha15;}
					if (a2==4 && a3==2){theta[2][4]=calpha23;}
					if (a2==4 && a4==2){theta[2][4]=calpha24;}
					if (a2==4 && a5==2){theta[2][4]=calpha25;}
					if (a3==4 && a4==2){theta[2][4]=calpha34;}
					if (a3==4 && a5==2){theta[2][4]=calpha35;}
					if (a4==4 && a5==2){theta[2][4]=calpha45;}
					
					if (a1==3 && a2==4){theta[3][4]=calpha12;}
					if (a1==3 && a3==4){theta[3][4]=calpha13;}
					if (a1==3 && a4==4){theta[3][4]=calpha14;}
					if (a1==3 && a5==4){theta[3][4]=calpha15;}
					if (a2==3 && a3==4){theta[3][4]=calpha23;}
					if (a2==3 && a4==4){theta[3][4]=calpha24;}
					if (a2==3 && a5==4){theta[3][4]=calpha25;}
					if (a3==3 && a4==4){theta[3][4]=calpha34;}
					if (a3==3 && a5==4){theta[3][4]=calpha35;}
					if (a4==3 && a5==4){theta[3][4]=calpha45;}
					
					if (a1==4 && a2==3){theta[3][4]=calpha12;}
					if (a1==4 && a3==3){theta[3][4]=calpha13;}
					if (a1==4 && a4==3){theta[3][4]=calpha14;}
					if (a1==4 && a5==3){theta[3][4]=calpha15;}
					if (a2==4 && a3==3){theta[3][4]=calpha23;}
					if (a2==4 && a4==3){theta[3][4]=calpha24;}
					if (a2==4 && a5==3){theta[3][4]=calpha25;}
					if (a3==4 && a4==3){theta[3][4]=calpha34;}
					if (a3==4 && a5==3){theta[3][4]=calpha35;}
					if (a4==4 && a5==3){theta[3][4]=calpha45;}
					
					
					for (m=0; m<1; m++)
					{
						if(m==0){inicio_m=1;}
						if(a1>=5){a1=0;}
						if(a2>=5){a2=0;}
						if(a3>=5){a3=0;}
						if(a4>=5){a4=0;}
						if(a5>=5){a5=0;}
						
				      	if(m==0){k0[5]=min(&inicio_m,k1,&a1,k2,&a2,k3,&a3,k4,&a4,k5,&a5,5);}
						
						if (a1==1 && a2==5){theta[1][5]=calpha12;}
						if (a1==1 && a3==5){theta[1][5]=calpha13;}
						if (a1==1 && a4==5){theta[1][5]=calpha14;}
						if (a1==1 && a5==5){theta[1][5]=calpha15;}
						if (a2==1 && a3==5){theta[1][5]=calpha23;}
						if (a2==1 && a4==5){theta[1][5]=calpha24;}
						if (a2==1 && a5==5){theta[1][5]=calpha25;}
						if (a3==1 && a4==5){theta[1][5]=calpha34;}
						if (a3==1 && a5==5){theta[1][5]=calpha35;}
						if (a4==1 && a5==5){theta[1][5]=calpha45;}
						
						if (a1==5 && a2==1){theta[1][5]=calpha12;}
						if (a1==5 && a3==1){theta[1][5]=calpha13;}
						if (a1==5 && a4==1){theta[1][5]=calpha14;}
						if (a1==5 && a5==1){theta[1][5]=calpha15;}
						if (a2==5 && a3==1){theta[1][5]=calpha23;}
						if (a2==5 && a4==1){theta[1][5]=calpha24;}
						if (a2==5 && a5==1){theta[1][5]=calpha25;}
						if (a3==5 && a4==1){theta[1][5]=calpha34;}
						if (a3==5 && a5==1){theta[1][5]=calpha35;}
						if (a4==5 && a5==1){theta[1][5]=calpha45;}
						
						if (a1==2 && a2==5){theta[2][5]=calpha12;}
						if (a1==2 && a3==5){theta[2][5]=calpha13;}
						if (a1==2 && a4==5){theta[2][5]=calpha14;}
						if (a1==2 && a5==5){theta[2][5]=calpha15;}
						if (a2==2 && a3==5){theta[2][5]=calpha23;}
						if (a2==2 && a4==5){theta[2][5]=calpha24;}
						if (a2==2 && a5==5){theta[2][5]=calpha25;}
						if (a3==2 && a4==5){theta[2][5]=calpha34;}
						if (a3==2 && a5==5){theta[2][5]=calpha35;}
						if (a4==2 && a5==5){theta[2][5]=calpha45;}
						
						if (a1==5 && a2==2){theta[2][5]=calpha12;}
						if (a1==5 && a3==2){theta[2][5]=calpha13;}
						if (a1==5 && a4==2){theta[2][5]=calpha14;}
						if (a1==5 && a5==2){theta[2][5]=calpha15;}
						if (a2==5 && a3==2){theta[2][5]=calpha23;}
						if (a2==5 && a4==2){theta[2][5]=calpha24;}
						if (a2==5 && a5==2){theta[2][5]=calpha25;}
						if (a3==5 && a4==2){theta[2][5]=calpha34;}
						if (a3==5 && a5==2){theta[2][5]=calpha35;}
						if (a4==5 && a5==2){theta[2][5]=calpha45;}
						
						if (a1==3 && a2==5){theta[3][5]=calpha12;}
						if (a1==3 && a3==5){theta[3][5]=calpha13;}
						if (a1==3 && a4==5){theta[3][5]=calpha14;}
						if (a1==3 && a5==5){theta[3][5]=calpha15;}
						if (a2==3 && a3==5){theta[3][5]=calpha23;}
						if (a2==3 && a4==5){theta[3][5]=calpha24;}
						if (a2==3 && a5==5){theta[3][5]=calpha25;}
						if (a3==3 && a4==5){theta[3][5]=calpha34;}
						if (a3==3 && a5==5){theta[3][5]=calpha35;}
						if (a4==3 && a5==5){theta[3][5]=calpha45;}
						
						if (a1==5 && a2==3){theta[3][5]=calpha12;}
						if (a1==5 && a3==3){theta[3][5]=calpha13;}
						if (a1==5 && a4==3){theta[3][5]=calpha14;}
						if (a1==5 && a5==3){theta[3][5]=calpha15;}
						if (a2==5 && a3==3){theta[3][5]=calpha23;}
						if (a2==5 && a4==3){theta[3][5]=calpha24;}
						if (a2==5 && a5==3){theta[3][5]=calpha25;}
						if (a3==5 && a4==3){theta[3][5]=calpha34;}
						if (a3==5 && a5==3){theta[3][5]=calpha35;}
						if (a4==5 && a5==3){theta[3][5]=calpha45;}
						
						if (a1==4 && a2==5){theta[4][5]=calpha12;}
						if (a1==4 && a3==5){theta[4][5]=calpha13;}
						if (a1==4 && a4==5){theta[4][5]=calpha14;}
						if (a1==4 && a5==5){theta[4][5]=calpha15;}
						if (a2==4 && a3==5){theta[4][5]=calpha23;}
						if (a2==4 && a4==5){theta[4][5]=calpha24;}
						if (a2==4 && a5==5){theta[4][5]=calpha25;}
						if (a3==4 && a4==5){theta[4][5]=calpha34;}
						if (a3==4 && a5==5){theta[4][5]=calpha35;}
						if (a4==4 && a5==5){theta[4][5]=calpha45;}
						
						if (a1==5 && a2==4){theta[4][5]=calpha12;}
						if (a1==5 && a3==4){theta[4][5]=calpha13;}
						if (a1==5 && a4==4){theta[4][5]=calpha14;}
						if (a1==5 && a5==4){theta[4][5]=calpha15;}
						if (a2==5 && a3==4){theta[4][5]=calpha23;}
						if (a2==5 && a4==4){theta[4][5]=calpha24;}
						if (a2==5 && a5==4){theta[4][5]=calpha25;}
						if (a3==5 && a4==4){theta[4][5]=calpha34;}
						if (a3==5 && a5==4){theta[4][5]=calpha35;}
						if (a4==5 && a5==4){theta[4][5]=calpha45;}
						
			            f=f+G5(k0[1],k0[2],k0[3],k0[4],k0[5],theta[1][2],theta[1][3],theta[1][4],theta[1][5],theta[2][3],theta[2][4],theta[2][5],theta[3][4],theta[3][5],theta[4][5]);
						//	printf("(%lf %lf %lf %lf %lf, %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf (%d,%d,%d,%d,%d) \n",k0[1],k0[2],k0[3],k0[4],k0[5],theta[1][2],theta[1][3],theta[1][4],theta[1][5],theta[2][3],theta[2][4],theta[2][5],theta[3][4],theta[3][5],theta[4][5],a1,a2,a3,a4,a5);
					}
				}
			}
		}
	}
	
	f=f/120.; 
	return f;
}



#define Pi (4.*atan(1.))
#include "structures.h"
#include "functions.h"
#include "pt_kernels.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void Function_P44_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	//full P44

	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
	int N= params_function.Nlin;
	double q1=x[0];
	double q2=x[1];
	double q3=x[2];
	double ctheta1=x[3];
	double ctheta2=x[4];
	double ctheta3=x[5];
	double phi2=x[6];
	double phi3=x[7];
	
	double P3=P_interpolLOG(q3,k_L,Pk_L,N);
	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
	
//        double P3=P_interpol(q3,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);

	double epsilon=10e-5;
	
	double modk123;
	double ctheta12;
	double ctheta13;
	double ctheta23;
	double ctheta1_123;
	double ctheta2_123;
	double ctheta3_123;
	
	ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	ctheta13=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta3*ctheta3)*cos(phi3)+ctheta1*ctheta3;
	ctheta23=sqrt(1-ctheta2*ctheta2)*sqrt(1-ctheta3*ctheta3)*(cos(phi2)*cos(phi3)+sin(phi2)*sin(phi3))+ctheta2*ctheta3;
	
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}
	
	if(ctheta13>1.0 && ctheta13<1.0+epsilon){ ctheta13=1.0;}
	if(ctheta13<-1.0 && ctheta13>-1.0-epsilon){ ctheta13=-1.0;}
	
	if(ctheta23>1.0 && ctheta23<1.0+epsilon){ ctheta23=1.0;}
	if(ctheta23<-1.0 && ctheta23>-1.0-epsilon){ ctheta23=-1.0;}
	
	modk123=sqrt(k0*k0+q1*q1+q2*q2+q3*q3-2.*k0*(q1*ctheta1+q2*ctheta2+q3*ctheta3)+2.*(q1*q2*ctheta12+q1*q3*ctheta13+q2*q3*ctheta23));
	
	ctheta1_123=(k0*ctheta1-q1-q2*ctheta12-q3*ctheta13)/modk123;
	
	if(ctheta1_123>1.0 && ctheta1_123<1.0+epsilon){ ctheta1_123=1.0;}
	if(ctheta1_123<-1.0 && ctheta1_123>-1.0-epsilon){ ctheta1_123=-1.0;}
	
	ctheta2_123=(k0*ctheta2-q1*ctheta12-q2-q3*ctheta23)/modk123;
	
	if(ctheta2_123>1.0 && ctheta2_123<1.0+epsilon){ ctheta2_123=1.0;}
	if(ctheta2_123<-1.0 && ctheta2_123>-1.0-epsilon){ ctheta2_123=-1.0;}
	
	ctheta3_123=(k0*ctheta3-q1*ctheta13-q2*ctheta23-q3)/modk123;
	
	if(ctheta3_123>1.0 && ctheta3_123<1.0+epsilon){ ctheta3_123=1.0;}
	if(ctheta3_123<-1.0 && ctheta3_123>-1.0-epsilon){ ctheta3_123=-1.0;}
	
	if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p44 ctheta12=%lf\n",ctheta12); }
	if( fabs(ctheta13)>1.0 ){ printf("\nError ang. a p44 ctheta13=%lf\n",ctheta13); }
	if( fabs(ctheta23)>1.0 ){ printf("\nError ang. a p44 ctheta23=%lf\n",ctheta23); }
	
	if( fabs(ctheta1_123)>1.0 ){ printf("\nError ang. a p44 ctheta1_123=%lf\n",ctheta1_123); }
	if( fabs(ctheta2_123)>1.0 ){ printf("\nError ang. a p44 ctheta2_123=%lf\n",ctheta2_123); }
	if( fabs(ctheta3_123)>1.0 ){ printf("\nError ang. a p44 ctheta3_123=%lf\n",ctheta3_123); }
	
	double f=pow(Fs4(q1,q2,q3,modk123,ctheta12,ctheta13,ctheta1_123,ctheta23,ctheta2_123,ctheta3_123),2)*P_interpolLOG(modk123,k_L,Pk_L,N);
//        double f=pow(Fs4(q1,q2,q3,modk123,ctheta12,ctheta13,ctheta1_123,ctheta23,ctheta2_123,ctheta3_123),2)*P_interpol(modk123,k_L,Pk_L,N);	
	fval[0]=24.*(2.*Pi)*q1*q1*P1*q2*q2*P2*q3*q3*P3*f;
        if(modk123==0){fval[0]=0;}

}


void Function_P33_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;	
        int N= params_function.Nlin;
	double q1=x[0];
	double q2=x[1];
	double ctheta1=x[2];
	double ctheta2=x[3];
	double phi2=x[4];
	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);

	double epsilon=10e-10;
	
	double modk12;
	double cthetak12_1;
	double cthetak12_2;
	double ctheta12;
	
	ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}
	
	modk12=pow(fabs(k0*k0+q1*q1+q2*q2-2.*k0*q1*ctheta1-2.*k0*q2*ctheta2+2.*q1*q2*ctheta12),0.5);
	
	cthetak12_1=(k0*ctheta1-q1-q2*ctheta12)/modk12;
	
	if(cthetak12_1>1.0 && cthetak12_1<1.0+epsilon){ cthetak12_1=1.0;}
	if(cthetak12_1<-1.0 && cthetak12_1>-1.0-epsilon){ cthetak12_1=-1.0;}
	
	cthetak12_2=(k0*ctheta2-q2-q1*ctheta12)/modk12;
	
	if(cthetak12_2>1.0 && cthetak12_2<1.0+epsilon){ cthetak12_2=1.0;}
	if(cthetak12_2<-1.0 && cthetak12_2>-1.0-epsilon){ cthetak12_2=-1.0;}
	
	if( fabs(cthetak12_1)>1.0 ){ printf("\nError ang. a p33 cthetak12_1=%lf\n",cthetak12_1); /*return 0;*/}
	if( fabs(cthetak12_2)>1.0 ){ printf("\nError ang. a p33 cthetak12_2=%lf\n",cthetak12_2); /*return 0;*/}
	if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p33 ctheta12=%lf\n",ctheta12); /*return 0;*/}
	
	double f;
        double Pkqq=P_interpolLOG(modk12,k_L,Pk_L,N);
//        double Pkqq=P_interpol(modk12,k_L,Pk_L,N);

        f=Fs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Fs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Pkqq;

	fval[0]=6.*4.*Pi*q1*q1*P1*q2*q2*P2*f;
	
	
}

void Function_P33_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q1=x[0];
	double q2=x[1];
	double ctheta1=x[2];
	double ctheta2=x[3];
	double phi2=x[4];
	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
	
	double epsilon=10e-10;
	
	double modk12;
	double cthetak12_1;
	double cthetak12_2;
	double ctheta12;
	
	ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}
	
	modk12=pow(fabs(k0*k0+q1*q1+q2*q2-2.*k0*q1*ctheta1-2.*k0*q2*ctheta2+2.*q1*q2*ctheta12),0.5);
	
	cthetak12_1=(k0*ctheta1-q1-q2*ctheta12)/modk12;
	
	if(cthetak12_1>1.0 && cthetak12_1<1.0+epsilon){ cthetak12_1=1.0;}
	if(cthetak12_1<-1.0 && cthetak12_1>-1.0-epsilon){ cthetak12_1=-1.0;}
	
	cthetak12_2=(k0*ctheta2-q2-q1*ctheta12)/modk12;
	
	if(cthetak12_2>1.0 && cthetak12_2<1.0+epsilon){ cthetak12_2=1.0;}
	if(cthetak12_2<-1.0 && cthetak12_2>-1.0-epsilon){ cthetak12_2=-1.0;}
	
	if( fabs(cthetak12_1)>1.0 ){ printf("\nError ang. a p33 cthetak12_1=%lf\n",cthetak12_1); /*return 0;*/}
	if( fabs(cthetak12_2)>1.0 ){ printf("\nError ang. a p33 cthetak12_2=%lf\n",cthetak12_2); /*return 0;*/}
	if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p33 ctheta12=%lf\n",ctheta12); /*return 0;*/}
	
	double f;
        double Pkqq=P_interpolLOG(modk12,k_L,Pk_L,N);
//        double Pkqq=P_interpol(modk12,k_L,Pk_L,N);

                     f=Gs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Gs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Pkqq;

	fval[0]=6.*4.*Pi*q1*q1*P1*q2*q2*P2*f;
	
	
}

void Function_P33_do(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q1=x[0];
	double q2=x[1];
	double ctheta1=x[2];
	double ctheta2=x[3];
	double phi2=x[4];
	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
	
	double epsilon=10e-10;
	
	double modk12;
	double cthetak12_1;
	double cthetak12_2;
	double ctheta12;
	
	ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}
	
	modk12=pow(fabs(k0*k0+q1*q1+q2*q2-2.*k0*q1*ctheta1-2.*k0*q2*ctheta2+2.*q1*q2*ctheta12),0.5);
	
	cthetak12_1=(k0*ctheta1-q1-q2*ctheta12)/modk12;
	
	if(cthetak12_1>1.0 && cthetak12_1<1.0+epsilon){ cthetak12_1=1.0;}
	if(cthetak12_1<-1.0 && cthetak12_1>-1.0-epsilon){ cthetak12_1=-1.0;}
	
	cthetak12_2=(k0*ctheta2-q2-q1*ctheta12)/modk12;
	
	if(cthetak12_2>1.0 && cthetak12_2<1.0+epsilon){ cthetak12_2=1.0;}
	if(cthetak12_2<-1.0 && cthetak12_2>-1.0-epsilon){ cthetak12_2=-1.0;}
	
	if( fabs(cthetak12_1)>1.0 ){ printf("\nError ang. a p33 cthetak12_1=%lf\n",cthetak12_1); /*return 0;*/}
	if( fabs(cthetak12_2)>1.0 ){ printf("\nError ang. a p33 cthetak12_2=%lf\n",cthetak12_2); /*return 0;*/}
	if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p33 ctheta12=%lf\n",ctheta12); /*return 0;*/}
	
	double f;
        double Pkqq=P_interpolLOG(modk12,k_L,Pk_L,N);
//        double Pkqq=P_interpol(modk12,k_L,Pk_L,N);

                     f=Fs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Gs3(q1,q2,modk12,ctheta12,cthetak12_1,cthetak12_2)*Pkqq;

	fval[0]=6.*4.*Pi*q1*q1*P1*q2*q2*P2*f;
	
	
}

void Function_P22_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q=x[0];
	double ctheta=x[1];
	double P=P_interpolLOG(q,k_L,Pk_L,N);
//        double P=P_interpol(q,k_L,Pk_L,N);

	
	double epsilon=1e-10;
	double mod_kq=pow(q*q+k0*k0-2.*q*k0*ctheta,0.5);//modulo de k-q
	//if(mod_kq<epsilon){printf("\n Error mod_kq->0 (%lf)\n",mod_kq);}
	
	double calpha_q_kq=(k0*ctheta-q)/mod_kq; //anglre entre q i k-q
	
	if(calpha_q_kq>1.0 && calpha_q_kq<1.0+epsilon){ calpha_q_kq=1.0;}
	if(calpha_q_kq<-1.0 && calpha_q_kq>-1.0-epsilon){ calpha_q_kq=-1.0;}
	
	
	if( fabs(calpha_q_kq)>1.0 ){ printf("\nError ang. a p22 calpha_k_q=%lf\n",calpha_q_kq);}
	
	
	double P_kq=P_interpolLOG(mod_kq,k_L,Pk_L,N);
//        double P_kq=P_interpol(mod_kq,k_L,Pk_L,N);

	
	fval[0]=2.*Pi*q*q*2.*P_interpolLOG(q,k_L,Pk_L,N)*pow(Fs2(q,mod_kq,calpha_q_kq),2)*P_kq;
//        fval[0]=2.*Pi*q*q*2.*P_interpol(q,k_L,Pk_L,N)*pow(Fs2(q,mod_kq,calpha_q_kq),2)*P_kq;
	
	
}

void Function_P22_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q=x[0];
	double ctheta=x[1];
	double P=P_interpolLOG(q,k_L,Pk_L,N);
//        double P=P_interpol(q,k_L,Pk_L,N);
	
	double epsilon=1e-10;
	double mod_kq=pow(q*q+k0*k0-2.*q*k0*ctheta,0.5);//modulo de k-q
	//if(mod_kq<epsilon){printf("\n Error mod_kq->0 (%lf)\n",mod_kq);}
	
	double calpha_q_kq=(k0*ctheta-q)/mod_kq; //anglre entre q i k-q
	
	if(calpha_q_kq>1.0 && calpha_q_kq<1.0+epsilon){ calpha_q_kq=1.0;}
	if(calpha_q_kq<-1.0 && calpha_q_kq>-1.0-epsilon){ calpha_q_kq=-1.0;}
	
	
	if( fabs(calpha_q_kq)>1.0 ){ printf("\nError ang. a p22 calpha_k_q=%lf\n",calpha_q_kq);}
	
	
	double P_kq=P_interpolLOG(mod_kq,k_L,Pk_L,N);
//        double P_kq=P_interpol(mod_kq,k_L,Pk_L,N);

	
	fval[0]=2.*Pi*q*q*2.*P_interpolLOG(q,k_L,Pk_L,N)*pow(Gs2(q,mod_kq,calpha_q_kq),2)*P_kq;
//        fval[0]=2.*Pi*q*q*2.*P_interpol(q,k_L,Pk_L,N)*pow(Gs2(q,mod_kq,calpha_q_kq),2)*P_kq;

	
	
}

void Function_P22_do(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q=x[0];
	double ctheta=x[1];
	double P=P_interpolLOG(q,k_L,Pk_L,N);
//        double P=P_interpol(q,k_L,Pk_L,N);
        
	
	double epsilon=1e-10;
	double mod_kq=pow(q*q+k0*k0-2.*q*k0*ctheta,0.5);//modulo de k-q
	//if(mod_kq<epsilon){printf("\n Error mod_kq->0 (%lf)\n",mod_kq);}
	
	double calpha_q_kq=(k0*ctheta-q)/mod_kq; //anglre entre q i k-q
	
	if(calpha_q_kq>1.0 && calpha_q_kq<1.0+epsilon){ calpha_q_kq=1.0;}
	if(calpha_q_kq<-1.0 && calpha_q_kq>-1.0-epsilon){ calpha_q_kq=-1.0;}
	
	
	if( fabs(calpha_q_kq)>1.0 ){ printf("\nError ang. a p22 calpha_k_q=%lf\n",calpha_q_kq);}
	
	
	double P_kq=P_interpolLOG(mod_kq,k_L,Pk_L,N);
//        double P_kq=P_interpol(mod_kq,k_L,Pk_L,N);

	
	fval[0]=2.*Pi*q*q*2.*P_interpolLOG(q,k_L,Pk_L,N)*Fs2(q,mod_kq,calpha_q_kq)*Gs2(q,mod_kq,calpha_q_kq)*P_kq;
//        fval[0]=2.*Pi*q*q*2.*P_interpol(q,k_L,Pk_L,N)*Fs2(q,mod_kq,calpha_q_kq)*Gs2(q,mod_kq,calpha_q_kq)*P_kq;
	
	
}

void Function_P13_d(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q=x[0];
	double ctheta=x[1];
	double P=P_interpolLOG(q,k_L,Pk_L,N);
//        double P=P_interpol(q,k_L,Pk_L,N);
	
	fval[0]=2.*Pi*q*q*6.*P_interpolLOG(q,k_L,Pk_L,N)*P_interpolLOG(k0,k_L,Pk_L,N)*Fs3(k0,q,q,ctheta,-ctheta,-1.0);
//        fval[0]=2.*Pi*q*q*6.*P_interpolLOG(q,k_L,Pk_L,N)*P_interpolLOG(k0,k_L,Pk_L,N)*Fs3(k0,q,q,ctheta,-ctheta,-1.0);
	
	
}

void Function_P13_o(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q=x[0];
	double ctheta=x[1];
	double P=P_interpolLOG(q,k_L,Pk_L,N);
//      double P=P_interpol(q,k_L,Pk_L,N);

	
	fval[0]=2.*Pi*q*q*6.*P_interpolLOG(q,k_L,Pk_L,N)*P_interpolLOG(k0,k_L,Pk_L,N)*Gs3(k0,q,q,ctheta,-ctheta,-1.0);
//        fval[0]=2.*Pi*q*q*6.*P_interpolLOG(q,k_L,Pk_L,N)*P_interpolLOG(k0,k_L,Pk_L,N)*Gs3(k0,q,q,ctheta,-ctheta,-1.0);
	
	
}

void Function_P15_d(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	double epsilon=1e-10;
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q1=x[0];
	double q2=x[1];
	double ctheta1=x[2];
	double ctheta2=x[3];
	double phi2=x[4];

	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
	double P0=P_interpolLOG(k0,k_L,Pk_L,N);

//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P0=P_interpol(k0,k_L,Pk_L,N);

	double ctheta12;
	
	ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}
	
	if (fabs(ctheta1)>1 || fabs(ctheta2)>1 || fabs(ctheta12)>1) {printf("error en P15, (%lf,%lf,%lf)\n",ctheta1,ctheta2, ctheta12);	}


	fval[0]=30.*(2.*2.*Pi)*q1*q1*q2*q2*P1*P2*P0*Fs5(k0,q1,q1,q2,q2,ctheta1,-ctheta1, ctheta2,-ctheta2,-1.,ctheta12,-ctheta12,-ctheta12,ctheta12,-1.);

}

void Function_P15_o(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double epsilon=1e-10;
        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;

        double q1=x[0];
        double q2=x[1];
        double ctheta1=x[2];
        double ctheta2=x[3];
        double phi2=x[4];

        double P2=P_interpolLOG(q2,k_L,Pk_L,N);
        double P1=P_interpolLOG(q1,k_L,Pk_L,N);
        double P0=P_interpolLOG(k0,k_L,Pk_L,N);

//        double P2=P_interpol(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P0=P_interpol(k0,k_L,Pk_L,N);


        double ctheta12;

        ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
        if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
        if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}

        if (fabs(ctheta1)>1 || fabs(ctheta2)>1 || fabs(ctheta12)>1) {printf("error en P15, (%lf,%lf,%lf)\n",ctheta1,ctheta2, ctheta12); }


        fval[0]=30.*(2.*2.*Pi)*q1*q1*q2*q2*P1*P2*P0*Gs5(k0,q1,q1,q2,q2,ctheta1,-ctheta1, ctheta2,-ctheta2,-1.,ctheta12,-ctheta12,-ctheta12,ctheta12,-1.);

}


void Function_P24_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
	f_params params_function = *(f_params *) fdata;//cast from void to f_params
	double epsilon=1e-5;
	double k0= params_function.k;
	double *k_L= params_function.K_L;
	double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;
	
	double q1=x[0];
	double q2=x[1];
	double ctheta1=x[2];
	double ctheta2=x[3];
	double phi2=x[4];
	double P1=P_interpolLOG(q1,k_L,Pk_L,N);
	double P2=P_interpolLOG(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);

	
	double ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;
	
	if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
	if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}

	
	double mod_kq1=pow(fabs(q1*q1+k0*k0-2.*q1*k0*ctheta1),0.5);//modulo de k-q1
	double Pkq1=P_interpolLOG(mod_kq1,k_L,Pk_L,N);
//        double Pkq1=P_interpol(mod_kq1,k_L,Pk_L.N);

	
	double calpha_q1_kq1=(k0*ctheta1-q1)/mod_kq1; //anglre entre q1 i k-q1
	double calpha_q2_kq1=(k0*ctheta2-q1*ctheta12)/mod_kq1; //anglre entre q2 i k-q1

	if(calpha_q1_kq1>1.0 && calpha_q1_kq1<1.0+epsilon){ calpha_q1_kq1=1.0;}
	if(calpha_q1_kq1<-1.0 && calpha_q1_kq1>-1.0-epsilon){ calpha_q1_kq1=-1.0;}

	if(calpha_q2_kq1>1.0 && calpha_q2_kq1<1.0+epsilon){ calpha_q2_kq1=1.0;}
	if(calpha_q2_kq1<-1.0 && calpha_q2_kq1>-1.0-epsilon){ calpha_q2_kq1=-1.0;}

	
	if( fabs(calpha_q1_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q1=%lf\n",calpha_q1_kq1);}
	if( fabs(calpha_q2_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q2=%lf\n",calpha_q2_kq1);}
	if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p24 calpha_k_q=%lf\n",ctheta12);}

	
	fval[0]=24.*q1*q1*q2*q2*4.*Pi*P1*P2*Pkq1*Fs2(q1,mod_kq1,calpha_q1_kq1)*Fs4(q1,mod_kq1,q2,q2,calpha_q1_kq1,ctheta12,-ctheta12,calpha_q2_kq1,-calpha_q2_kq1,-1.);
	
	
}

void Function_P24_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double epsilon=1e-5;
        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;

        double q1=x[0];
        double q2=x[1];
        double ctheta1=x[2];
        double ctheta2=x[3];
        double phi2=x[4];
        double P1=P_interpolLOG(q1,k_L,Pk_L,N);
        double P2=P_interpolLOG(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);

        double ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;

        if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
        if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}


        double mod_kq1=pow(fabs(q1*q1+k0*k0-2.*q1*k0*ctheta1),0.5);//modulo de k-q1
        double Pkq1=P_interpolLOG(mod_kq1,k_L,Pk_L,N);
//        double Pkq1=P_interpol(mod_kq1,k_L,Pk_L,N);

        double calpha_q1_kq1=(k0*ctheta1-q1)/mod_kq1; //anglre entre q1 i k-q1
        double calpha_q2_kq1=(k0*ctheta2-q1*ctheta12)/mod_kq1; //anglre entre q2 i k-q1

        if(calpha_q1_kq1>1.0 && calpha_q1_kq1<1.0+epsilon){ calpha_q1_kq1=1.0;}
        if(calpha_q1_kq1<-1.0 && calpha_q1_kq1>-1.0-epsilon){ calpha_q1_kq1=-1.0;}

        if(calpha_q2_kq1>1.0 && calpha_q2_kq1<1.0+epsilon){ calpha_q2_kq1=1.0;}
        if(calpha_q2_kq1<-1.0 && calpha_q2_kq1>-1.0-epsilon){ calpha_q2_kq1=-1.0;}


        if( fabs(calpha_q1_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q1=%lf\n",calpha_q1_kq1);}
        if( fabs(calpha_q2_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q2=%lf\n",calpha_q2_kq1);}
        if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p24 calpha_k_q=%lf\n",ctheta12);}


        fval[0]=24.*q1*q1*q2*q2*4.*Pi*P1*P2*Pkq1*Gs2(q1,mod_kq1,calpha_q1_kq1)*Gs4(q1,mod_kq1,q2,q2,calpha_q1_kq1,ctheta12,-ctheta12,calpha_q2_kq1,-calpha_q2_kq1,-1.);


}

void Function_P24_do1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double epsilon=1e-5;
        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;

        double q1=x[0];
        double q2=x[1];
        double ctheta1=x[2];
        double ctheta2=x[3];
        double phi2=x[4];
        double P1=P_interpolLOG(q1,k_L,Pk_L,N);
        double P2=P_interpolLOG(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);

        double ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;

        if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
        if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}


        double mod_kq1=pow(fabs(q1*q1+k0*k0-2.*q1*k0*ctheta1),0.5);//modulo de k-q1
        double Pkq1=P_interpolLOG(mod_kq1,k_L,Pk_L,N);
//        double Pkq1=P_interpol(mod_kq1,k_L,Pk_L,N);


        double calpha_q1_kq1=(k0*ctheta1-q1)/mod_kq1; //anglre entre q1 i k-q1
        double calpha_q2_kq1=(k0*ctheta2-q1*ctheta12)/mod_kq1; //anglre entre q2 i k-q1

        if(calpha_q1_kq1>1.0 && calpha_q1_kq1<1.0+epsilon){ calpha_q1_kq1=1.0;}
        if(calpha_q1_kq1<-1.0 && calpha_q1_kq1>-1.0-epsilon){ calpha_q1_kq1=-1.0;}

        if(calpha_q2_kq1>1.0 && calpha_q2_kq1<1.0+epsilon){ calpha_q2_kq1=1.0;}
        if(calpha_q2_kq1<-1.0 && calpha_q2_kq1>-1.0-epsilon){ calpha_q2_kq1=-1.0;}


        if( fabs(calpha_q1_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q1=%lf\n",calpha_q1_kq1);}
        if( fabs(calpha_q2_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q2=%lf\n",calpha_q2_kq1);}
        if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p24 calpha_k_q=%lf\n",ctheta12);}


        fval[0]=24.*q1*q1*q2*q2*4.*Pi*P1*P2*Pkq1*Gs2(q1,mod_kq1,calpha_q1_kq1)*Fs4(q1,mod_kq1,q2,q2,calpha_q1_kq1,ctheta12,-ctheta12,calpha_q2_kq1,-calpha_q2_kq1,-1.);


}

void Function_P24_do2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
        f_params params_function = *(f_params *) fdata;//cast from void to f_params
        double epsilon=1e-5;
        double k0= params_function.k;
        double *k_L= params_function.K_L;
        double *Pk_L= params_function.PK_L;
        int N= params_function.Nlin;

        double q1=x[0];
        double q2=x[1];
        double ctheta1=x[2];
        double ctheta2=x[3];
        double phi2=x[4];
        double P1=P_interpolLOG(q1,k_L,Pk_L,N);
        double P2=P_interpolLOG(q2,k_L,Pk_L,N);
//        double P1=P_interpol(q1,k_L,Pk_L,N);
//        double P2=P_interpol(q2,k_L,Pk_L,N);

        double ctheta12=sqrt(1.-ctheta1*ctheta1)*sqrt(1.-ctheta2*ctheta2)*cos(phi2)+ctheta1*ctheta2;

        if(ctheta12>1.0 && ctheta12<1.0+epsilon){ ctheta12=1.0;}
        if(ctheta12<-1.0 && ctheta12>-1.0-epsilon){ ctheta12=-1.0;}


        double mod_kq1=pow(fabs(q1*q1+k0*k0-2.*q1*k0*ctheta1),0.5);//modulo de k-q1
        double Pkq1=P_interpolLOG(mod_kq1,k_L,Pk_L,N);
//        double Pkq1=P_interpol(mod_kq1,k_L,Pk_L,N);

        double calpha_q1_kq1=(k0*ctheta1-q1)/mod_kq1; //anglre entre q1 i k-q1
        double calpha_q2_kq1=(k0*ctheta2-q1*ctheta12)/mod_kq1; //anglre entre q2 i k-q1

        if(calpha_q1_kq1>1.0 && calpha_q1_kq1<1.0+epsilon){ calpha_q1_kq1=1.0;}
        if(calpha_q1_kq1<-1.0 && calpha_q1_kq1>-1.0-epsilon){ calpha_q1_kq1=-1.0;}

        if(calpha_q2_kq1>1.0 && calpha_q2_kq1<1.0+epsilon){ calpha_q2_kq1=1.0;}
        if(calpha_q2_kq1<-1.0 && calpha_q2_kq1>-1.0-epsilon){ calpha_q2_kq1=-1.0;}


        if( fabs(calpha_q1_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q1=%lf\n",calpha_q1_kq1);}
        if( fabs(calpha_q2_kq1)>1.0 ){ printf("\nError ang. a p24 calpha_kq1_q2=%lf\n",calpha_q2_kq1);}
        if( fabs(ctheta12)>1.0 ){ printf("\nError ang. a p24 calpha_k_q=%lf\n",ctheta12);}


        fval[0]=24.*q1*q1*q2*q2*4.*Pi*P1*P2*Pkq1*Fs2(q1,mod_kq1,calpha_q1_kq1)*Gs4(q1,mod_kq1,q2,q2,calpha_q1_kq1,ctheta12,-ctheta12,calpha_q2_kq1,-calpha_q2_kq1,-1.);


}


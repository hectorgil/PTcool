/*
 * PT-cool 
 * All rights reserved
 * Author: Hector Gil Marin
 * Date: 8th Feb 2019
 * email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#include "structures.h"
#include "functions.h"
#include "functions_rsd.h"

#include "nlbias_integrals.h"

#include "pt_kernels.h"
#include "pt_integrals.h"

#include "tns_integrals.h"
#include "tns_kernels.h"

#define Pi (4.*atan(1.))


int main(int argc, char *argv[])
{
FILE *f;
char name_file_ini[200];
char name_input[200];
char name_id[200];
char path_output[200];
char name_output[200];
char interval_type[200];//lin or log
int i,j,imax;
int Nlines,N;
double kmin,k,kmax;
double interval;
double precision_rpt, precision_max,precision_tns,precision_nlb;
double *klin,*Plin,*Plin2pi3;
double kin,pin;
double kmin1,kmax1,kmax11;
double Sigma_v;
double trash;
double time_ini,time_fin;
double p0,p22dd,p22do,p22oo,p33dd,p33do,p33oo,p13d,p13o,p15d,p15o,p24dd,p24oo,p24do1,p24do2;

double *Pk1,*Pk2,*Pk3,*Pk4,*Pk5,*Pk6,*Pk7,*Pk8,*Pk9,*Pk10,*Pk11,*Pk12;

double Pb2_delta,Pb2_theta,Pbs2_delta,Pbs2_theta,Pb2s2,Pbs22,Pb22,sigma3;
double *Pk13,*Pk14,*Pk15,*Pk16,*Pk17,*Pk18,*Pk19,*Pk20;

double Taruya_A11,Taruya_A12,Taruya_A22,Taruya_A23,Taruya_A33,Taruya_B1_11,Taruya_B1_12,Taruya_B1_21,Taruya_B1_22,Taruya_B2_11,Taruya_B2_12,Taruya_B2_21,Taruya_B2_22,Taruya_B3_12,Taruya_B3_21,Taruya_B3_22,Taruya_B4_22;
double *Pk21,*Pk22,*Pk23,*Pk24,*Pk25,*Pk26,*Pk27,*Pk28,*Pk29,*Pk30,*Pk31,*Pk32,*Pk33,*Pk34,*Pk35,*Pk36,*Pk37;

double sigma8;

int tid;

printf("\t ==== PT-cool ====\n\n");
printf("\t *** All rights reserved (c) ***\n\n");
printf("\t *** Hector Gil Marin ***\n\n");
printf("\t *** Last Release: 8th Feb 2019 ***\n\n");
printf("\t *** email: hector.gil.marin@gmail.com or hectorgil@icc.ub.edu ***\n\n");
//return 0;
//read params file
sprintf(name_file_ini,argv[1]);
f=fopen(name_file_ini,"r");
if(f==NULL){printf("File %s not found...\t exiting now\n",name_file_ini);return 0;}
else{printf("Reading Inizialization file: %s\n\n",name_file_ini);}

fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %*s %*s %s\n",name_input);
printf("Intput linear power spectrum file: %s\n",name_input);
fscanf(f,"%*s %*s %*s %*s %s\n",path_output);
fscanf(f,"%*s %s\n\n",name_id);

fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %s\n",interval_type);
if( strcmp(interval_type, "lin") != 0 && strcmp(interval_type, "log10") != 0 && strcmp(interval_type, "linear") != 0){printf("Interval type must be either lin or log10. Read entry is %s. Exiting now...\n",interval_type);exit(0);}
printf("%s space among k-bins.\n",interval_type);
fscanf(f,"%*s %*s %lf\n",&kmin);
fscanf(f,"%*s %*s %lf\n",&kmax);
fscanf(f,"%*s %lf\n\n",&interval);
if(kmin>kmax){printf("Error, kmin>kmax (%lf>%lf). Exiting now...\n",kmin,kmax);exit(0);}
if(interval<=0){printf("Error, interval=<0. Exiting now...\n");exit(0);}
printf("k-range %lf<=k<=%lf. Interval Dk=%lf\n",kmin,kmax,interval);
fscanf(f,"%*s %*s\n");
fscanf(f,"%*s %*s %lf\n",&precision_max);
fscanf(f,"%*s %*s %lf\n",&precision_rpt);
fscanf(f,"%*s %*s %lf\n",&precision_nlb);
fscanf(f,"%*s %*s %lf\n",&precision_tns);

if(precision_max<=0 || precision_rpt<=0 || precision_nlb<=0 || precision_tns<=0){printf("Precision parameters must be >0. Exiting now...\n");exit(0);}

printf("Maximum precision: %e\n",precision_max);
printf("PT precision: %e\n",precision_rpt);
printf("NLB precision: %e\n",precision_nlb);
printf("TNS precision: %e\n",precision_tns);


fclose(f);

//////
sprintf(name_output,"%s/Perturbation_theory_%s.dat",path_output,name_id);
printf("Results will be written in %s ...\n",name_output);

f=fopen(name_output,"w");
if(f==NULL){printf("Error, %s could not be created. Exiting now...\n",name_output); exit(0);}
fclose(f);

f=fopen(name_input,"r");
if(f==NULL){printf("Error, %s not be found. Exiting now...\n",name_input); exit(0);}
fclose(f);
printf("Reading linear power spectrum file %s...\n",name_input);

Nlines=countlines(name_input);
klin =  (double*) calloc( Nlines, sizeof(double));
Plin =  (double*) calloc( Nlines, sizeof(double));
Plin2pi3 =  (double*) calloc( Nlines, sizeof(double));
f=fopen(name_input,"r");
for(i=0;i<Nlines;i++)
{
fscanf(f,"%lf %lf\n",&kin,&pin);//change this by a function
klin[i]=kin;
Plin[i]=pin*pow(2.*Pi,-3);//CAMB notation
Plin2pi3[i]=pin;//CAMB notation

}
fclose(f);
kmin1=klin[0];
kmax1=1.0;
kmax11=klin[Nlines-1];


double        XMIN44[8]={kmin1,kmin1,kmin1,-1,-1,-1,0,0};
double        XMAX44[8]={kmax1,kmax1,kmax1,1,1,1,2.*Pi,2.*Pi};
double        XMIN33[5]={kmin1,kmin1,-1,-1,0};
double        XMAX33[5]={kmax1,kmax1,1,1,Pi};
double        XMIN13[2]={kmin1,-1};
double        XMAX13[2]={kmax1,1};
double        XMIN22[2]={kmin1,-1};
double        XMAX22[2]={kmax1,1};
double        XMAX22alt[2]={kmax11,1};
double        XMIN[1]={kmin1};
double        XMAX[1]={kmax1};
double        XMAXalt[1]={kmax11};

if(strcmp(interval_type, "lin") == 0 || (strcmp(interval_type, "linear") == 0) )
{
N=(int)( (kmax-kmin)/interval  );
i=0;
imax=(int)((kmax-kmin)/interval)+1;



trash=fabs(imax-((kmax-kmin)/interval)+1);
if( fabs(trash-1)<0.0001){imax=imax+1;}
}
if(strcmp(interval_type,"log10") == 0)
{
N=(int)( (log10(kmax)-log10(kmin))/interval  );
i=0;
imax=(int)((log10(kmax)-log10(kmin))/interval)+1;

trash=fabs(imax-((log10(kmax)-log10(kmin))/interval)+1);
if( fabs(trash-1)<0.0001){imax=imax+1;}


}

if(imax<=0){printf("Error, imax=%d\n",imax);return 0;}
Pk1 = (double*) calloc( imax, sizeof(double));
Pk2 = (double*) calloc( imax, sizeof(double));
Pk3 = (double*) calloc( imax, sizeof(double));
Pk4 = (double*) calloc( imax, sizeof(double));
Pk5 = (double*) calloc( imax, sizeof(double));
Pk6 = (double*) calloc( imax, sizeof(double));
Pk7 = (double*) calloc( imax, sizeof(double));
Pk8 = (double*) calloc( imax, sizeof(double));
Pk9 = (double*) calloc( imax, sizeof(double));
Pk10 = (double*) calloc( imax, sizeof(double));
Pk11 = (double*) calloc( imax, sizeof(double));
Pk12 = (double*) calloc( imax, sizeof(double));
Pk13 = (double*) calloc( imax, sizeof(double));
Pk14 = (double*) calloc( imax, sizeof(double));
Pk15 = (double*) calloc( imax, sizeof(double));
Pk16 = (double*) calloc( imax, sizeof(double));
Pk17 = (double*) calloc( imax, sizeof(double));
Pk18 = (double*) calloc( imax, sizeof(double));
Pk19 = (double*) calloc( imax, sizeof(double));
Pk20 = (double*) calloc( imax, sizeof(double));
Pk21 = (double*) calloc( imax, sizeof(double));
Pk22 = (double*) calloc( imax, sizeof(double));
Pk23 = (double*) calloc( imax, sizeof(double));
Pk24 = (double*) calloc( imax, sizeof(double));
Pk25 = (double*) calloc( imax, sizeof(double));
Pk26 = (double*) calloc( imax, sizeof(double));
Pk27 = (double*) calloc( imax, sizeof(double));
Pk28 = (double*) calloc( imax, sizeof(double));
Pk29 = (double*) calloc( imax, sizeof(double));
Pk30 = (double*) calloc( imax, sizeof(double));
Pk31 = (double*) calloc( imax, sizeof(double));
Pk32 = (double*) calloc( imax, sizeof(double));
Pk33 = (double*) calloc( imax, sizeof(double));
Pk34 = (double*) calloc( imax, sizeof(double));
Pk35 = (double*) calloc( imax, sizeof(double));
Pk36 = (double*) calloc( imax, sizeof(double));
Pk37 = (double*) calloc( imax, sizeof(double));


        f_params *function_parameters;
time_ini=time(NULL);
printf("Performing integrals. This may take a while...\n");

#pragma omp parallel for private(sigma8,function_parameters,tid,i,k,p0,p22dd,p22do,p22oo,p33dd,p33do,p33oo,p13d,p13o,p15d,p15o,p24dd,p24oo,p24do1,p24do2,Pb2_delta,Pb2_theta,Pbs2_delta,Pbs2_theta,Pb2s2,Pbs22,Pb22,sigma3,Taruya_A11,Taruya_A12,Taruya_A22,Taruya_A23,Taruya_A33,Taruya_B1_11,Taruya_B1_12,Taruya_B1_21,Taruya_B1_22,Taruya_B2_11,Taruya_B2_12,Taruya_B2_21,Taruya_B2_22,Taruya_B3_12,Taruya_B3_21,Taruya_B3_22,Taruya_B4_22,trash) shared(interval,Nlines,klin,Plin,Plin2pi3,imax,kmin,precision_rpt,precision_max,precision_nlb,precision_tns,XMIN33,XMAX33,XMIN22,XMAX22,XMIN13,XMAX13,XMAX22alt,XMIN44,XMAX44,XMIN,XMAX,XMAXalt,Pk1,Pk2,Pk3,Pk4,Pk5,Pk6,Pk7,Pk8,Pk9,Pk10,Pk11,Pk12,Pk13,Pk14,Pk15,Pk16,Pk17,Pk18,Pk19,Pk20,Pk21,Pk22,Pk23,Pk24,Pk25,Pk26,Pk27,Pk28,Pk29,Pk30,Pk31,Pk32,Pk33,Pk34,Pk35,Pk36,Pk37,interval_type,kmin1,kmax1,kmax11)
for(i=0;i<imax;i++)
{
tid=omp_get_thread_num();

if(strcmp(interval_type, "lin") == 0 || (strcmp(interval_type, "linear") == 0) )
{
k=kmin+i*interval;
}

if(strcmp(interval_type,"log10") == 0)

{
k=log10(kmin)+i*interval;
k=pow(10,k);
}
        function_parameters = (f_params *) malloc(sizeof(f_params));        
        (*function_parameters).k=k;
        (*function_parameters).K_L=klin;
        (*function_parameters).PK_L=Plin;
        (*function_parameters).PK_L_2pi3=Plin2pi3;
        (*function_parameters).Nlin=Nlines;

if(i==0){
adapt_integrate(1, integrals8 , function_parameters, 1, XMIN, XMAXalt ,0, precision_max, precision_max, &sigma8, &trash);
Pk37[0]=sqrt(sigma8);
}

printf("[%d/%d] k->%lf (threadID=%d)\n",i,imax-1,k,tid);
//Perturbation Theory
       //adapt_integrate(1, Function_P44 , function_parameters, 8, XMIN44, XMAX44,0, precision_rpt, precision_rpt, &P44[i], &trash);//9 sin simetria

        adapt_integrate(1, Function_P33_dd , function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p33dd, &trash);//6 sin simetria
        adapt_integrate(1, Function_P22_dd , function_parameters, 2, XMIN22, XMAX22,0, precision_max, precision_max, &p22dd, &trash);//3 sin simetria
        adapt_integrate(1, Function_P13_d , function_parameters, 2, XMIN13, XMAX13,0, precision_max, precision_max, &p13d, &trash);//3 sin simetria

        adapt_integrate(1, Function_P33_oo , function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p33oo, &trash);//6 sin simetria
        adapt_integrate(1, Function_P22_oo , function_parameters, 2, XMIN22, XMAX22,0, precision_max, precision_max, &p22oo, &trash);//3 sin simetria
        adapt_integrate(1, Function_P13_o , function_parameters, 2, XMIN13, XMAX13,0, precision_max, precision_max, &p13o, &trash);//3 sin simetria

        adapt_integrate(1, Function_P33_do , function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p33do, &trash);//6 sin simetria
        adapt_integrate(1, Function_P22_do , function_parameters, 2, XMIN22, XMAX22,0, precision_max, precision_max, &p22do, &trash);//3 sin simetria

        adapt_integrate(1, Function_P24_dd, function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p24dd, &trash);//6 sin simetria
        adapt_integrate(1, Function_P24_oo, function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p24oo, &trash);//6 sin simetria
        adapt_integrate(1, Function_P24_do1, function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p24do1, &trash);//6 sin simetria
        adapt_integrate(1, Function_P24_do2, function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p24do2, &trash);//6 sin simetria
        adapt_integrate(1, Function_P15_d , function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p15d, &trash);//6 sin simetria
        adapt_integrate(1, Function_P15_o , function_parameters, 5, XMIN33, XMAX33,0, precision_rpt, precision_rpt, &p15o, &trash);//6 sin simetria

        p0=P_interpolLOG(k,klin,Plin,Nlines)*pow(2.*Pi,3);
        p13d=p13d*pow(2.*Pi,3);
        p13o=p13o*pow(2.*Pi,3);

        p22dd=p22dd*pow(2.*Pi,3);
        p33dd=p33dd*pow(2.*Pi,3);

        p22do=p22do*pow(2.*Pi,3);
        p33do=p33do*pow(2.*Pi,3);

        p22oo=p22oo*pow(2.*Pi,3);
        p33oo=p33oo*pow(2.*Pi,3);
        p24dd=p24dd*pow(2.*Pi,3);
        p24oo=p24oo*pow(2.*Pi,3);
        p24do1=p24do1*pow(2.*Pi,3);
        p24do2=p24do2*pow(2.*Pi,3);
        p15d=p15d*pow(2.*Pi,3);
        p15o=p15o*pow(2.*Pi,3);

Pk1[i]=p0;Pk2[i]=p22dd;Pk3[i]=p22do;Pk4[i]=p22oo;Pk5[i]=p33dd;Pk6[i]=p33do;Pk7[i]=p33oo;Pk8[i]=p13d;Pk9[i]=p13o;Pk10[i]=p15d;Pk11[i]=p15o;

//bias
         adapt_integrate(1,intPb2_delta,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pb2_delta,&trash);
         adapt_integrate(1,intPb2_theta,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pb2_theta,&trash);
         adapt_integrate(1,intPbs2_delta,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pbs2_delta,&trash);
         adapt_integrate(1,intPbs2_theta,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pbs2_theta,&trash);
         adapt_integrate(1,intPb2s2,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pb2s2,&trash);
         adapt_integrate(1,intPbs22,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pbs22,&trash);
         adapt_integrate(1,intPb22,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&Pb22,&trash);
         adapt_integrate(1,intsigma3,function_parameters,2,XMIN22,XMAX22alt,0,precision_nlb,precision_nlb,&sigma3,&trash);


sigma3=sigma3*P_interpolLOG(k,klin,Plin2pi3,Nlines);

Pk12[i]=Pb2_delta;Pk13[i]=Pb2_theta;Pk14[i]=Pbs2_delta;Pk15[i]=Pbs2_theta;Pk16[i]=Pb2s2;Pk17[i]=Pbs22;Pk18[i]=Pb22;Pk19[i]=sigma3;

//Taruya terms
               adapt_integrate(1, Function_A11 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_A11, &trash);
                adapt_integrate(1, Function_A12 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_A12, &trash);
                adapt_integrate(1, Function_A22 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_A22, &trash);
                adapt_integrate(1, Function_A23 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_A23, &trash);
                adapt_integrate(1, Function_A33 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_A33, &trash);
                adapt_integrate(1, Function_B1_11 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B1_11, &trash);
                adapt_integrate(1, Function_B1_12 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B1_12, &trash);
                adapt_integrate(1, Function_B1_21 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B1_21, &trash);
                adapt_integrate(1, Function_B1_22 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B1_22, &trash);
                adapt_integrate(1, Function_B2_11 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B2_11, &trash);
                adapt_integrate(1, Function_B2_12 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B2_12, &trash);
                adapt_integrate(1, Function_B2_21 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B2_21, &trash);
                adapt_integrate(1, Function_B2_22 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B2_22, &trash);
                adapt_integrate(1, Function_B3_12 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B3_12, &trash);
                adapt_integrate(1, Function_B3_21 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B3_21, &trash);
                adapt_integrate(1, Function_B3_22 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B3_22, &trash);
                adapt_integrate(1, Function_B4_22 , function_parameters, 2, XMIN22, XMAX22alt,0, precision_tns, precision_tns, &Taruya_B4_22, &trash);
Pk20[i]=Taruya_A11;Pk21[i]=Taruya_A12;Pk22[i]=Taruya_A22;Pk23[i]=Taruya_A23;Pk24[i]=Taruya_A33;Pk25[i]=Taruya_B1_11;Pk26[i]=Taruya_B1_12;Pk27[i]=Taruya_B1_21;Pk28[i]=Taruya_B1_22;Pk29[i]=Taruya_B2_11;Pk30[i]=Taruya_B2_12;Pk31[i]=Taruya_B2_21;Pk32[i]=Taruya_B2_22;Pk33[i]=Taruya_B3_12;Pk34[i]=Taruya_B3_21;Pk35[i]=Taruya_B3_22;Pk36[i]=Taruya_B4_22;

}

f=fopen(name_output,"w");
for(i=0;i<imax;i++)
{

if(strcmp(interval_type, "lin") == 0 || (strcmp(interval_type, "linear") == 0) )
{
k=kmin+i*interval;
}

if(strcmp(interval_type,"log10") == 0)
{
k=log10(kmin)+i*interval;
k=pow(10,k);
}


fprintf(f,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",k,Pk1[i],Pk2[i],Pk3[i],Pk4[i],Pk5[i],Pk6[i],Pk7[i],Pk8[i],Pk9[i],Pk10[i],Pk11[i],Pk12[i],Pk13[i],Pk14[i],Pk15[i],Pk16[i],Pk17[i],Pk18[i],Pk19[i],Pk20[i],Pk21[i],Pk22[i],Pk23[i],Pk24[i],Pk25[i],Pk26[i],Pk27[i],Pk28[i],Pk29[i],Pk30[i],Pk31[i],Pk32[i],Pk33[i],Pk34[i],Pk35[i],Pk36[i],Pk37[0]);
}
fclose(f);
time_fin=time(NULL);
time_ini=(time_fin-time_ini)/60.;//in minutes

if(time_ini<90)
{
printf("Computation finised in %lf minutes.\n",time_ini);
}
else
{
printf("Computation finised in %lf hours.\n",time_ini/60.);
}



free(Pk1);free(Pk2);free(Pk3);free(Pk4);free(Pk5);free(Pk6);free(Pk7);free(Pk8);free(Pk9);free(Pk10);free(Pk11);free(Pk12);free(Pk13);free(Pk14);free(Pk15);
free(Pk16);free(Pk17);free(Pk18);free(Pk19);free(Pk20);free(Pk21);free(Pk22);free(Pk23);free(Pk24);free(Pk25);free(Pk26);free(Pk27);free(Pk28);free(Pk29);free(Pk30);free(Pk31);
free(Pk32);free(Pk33);free(Pk34);free(Pk35);free(Pk36);free(Pk37);
free(klin);free(Plin);free(Plin2pi3);

printf("Success! Exiting...\n");

return 0;
}

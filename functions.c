#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define Pi (4.*atan(1.))

void string_copy(char *from, char *to) {

    while ((*to++ = *from++) != '\0')
        ;
}

int countlines(char *filename)
{
FILE* myfile = fopen(filename, "r");
if(myfile==NULL){printf("Error reading %s. Exiting now...\n",filename);exit(0);}
int ch, number_of_lines = -1;

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

long int countlinesLI(char *filename)
{
FILE* myfile = fopen(filename, "r");
if(myfile==NULL){printf("Error reading %s. Exiting now...\n",filename);exit(0);}
long int ch, number_of_lines = -1;

do
{
            ch = fgetc(myfile);
                    if(ch == '\n')
                                        number_of_lines++;
} while (ch != EOF);

if(ch != '\n' && number_of_lines != 0)
    number_of_lines++;

        fclose(myfile);
        return number_of_lines;
}

void freeTokens(double **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{    
     free(tokens[i]);
}
free(tokens);

}

void freeTokens2(double ***tokens, int N1, int N2)
{
int i,j;
for(i=0;i<N1;++i)
{
     for(j=0;j<N2;++j)
     {
        free(tokens[i][j]);
     }
     free(tokens[i]);
}
free(tokens);

}



void freeTokensInt(int **tokens, int N)
{
int i;
for(i=0;i<N;++i)
{
     free(tokens[i]);
}
free(tokens);
}


double P_interpolLOG(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N);
if(i==0)
{
//printf("Interpol %d %lf %lf %lf %lf %lf -> \t",i,k0,k[i],k[i+1],P[i],P[i+1]);
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
else
{
//printf("Interpol %d %lf %lf %lf %lf %lf- >\t",i,k0,k[i-1],k[i],P[i-1],P[i]);
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
P0=m*log10(k0)+n;
P0=pow(10.,P0);
}
//printf("%lf\n",P0);
return P0;
}

double P_interpol(double k0, double *k, double *P, int N)
{
double P0,m,n;
int i;

if(k0<k[0] || k0>k[N-1])
{
P0=0;//return 0 if k0 is outside k-range
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i] && i<N);
if(i==0)
{
//printf("Interpol %d %lf %lf %lf %lf %lf -> \t",i,k0,k[i],k[i+1],P[i],P[i+1]);
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}
else
{
//printf("Interpol %d %lf %lf %lf %lf %lf- >\t",i,k0,k[i-1],k[i],P[i-1],P[i]);
m=( (P[i]) - (P[i-1]) )/( (k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
P0=m*(k0)+n;
}

}
//printf("%lf\n",P0);
return P0;
}

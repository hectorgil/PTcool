#include "functions.h"
#include "pt_kernels.h"

void Function_P44_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//3loop

void Function_P33_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//2loop

void Function_P22_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//1loop

void Function_P24_dd(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P24_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P24_do2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P24_do1(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P13_d(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P15_d(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P15_o(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P33_do(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//2loop

void Function_P22_do(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//1loop

void Function_P13_o(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void Function_P33_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//2loop

void Function_P22_oo(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);//1loop


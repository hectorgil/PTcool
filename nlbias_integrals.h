#include "functions.h"
#include "pt_kernels.h"

void intPb2_delta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPb2_theta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPbs2_delta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPbs2_theta(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPb22(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPb2s2(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intPbs22(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

void intsigma3(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

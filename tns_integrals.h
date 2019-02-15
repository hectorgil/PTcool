#include "functions.h"
#include "pt_kernels.h"
#include "tns_kernels.h"

void Function_A11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_A12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_A22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_A23(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_A33(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B1_11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B1_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B1_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B1_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B2_11(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B2_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B2_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B2_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B3_12(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B3_21(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B3_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);

void Function_B4_22(unsigned ndim, const double *X, void *fdata, unsigned fdim, double *fval);


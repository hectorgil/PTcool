
void string_copy(char *from, char *to);

double P_interpol(double k0, double *k, double *P, int N);

double P_interpolLOG(double k0, double *k, double *P, int N);

int countlines(char *filename);

long int countlinesLI(char *filename);

void freeTokens(double** tokens, int N);

void freeTokensInt(int** tokens, int N);

void freeTokens2(double ***tokens, int N1, int N2);

double Wth(double x);

void integrals8(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);


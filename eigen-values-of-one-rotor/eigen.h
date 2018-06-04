#ifndef _EIGEN_H
#define _EIGEN_H 1

//Legendre polynomial "plgndr.c" is called externally
extern "C" double plgndr(int l, int m, double x);

//"gauleg.f" gives the quadrature points and weights of Legendre polynomial
extern "C" void gauleg(double x1, double x2, double *x, double *w, int n);

//vh2h2 ---> potential H2-H2: added by Tapas Sahoo
extern "C" void vinit_();
extern "C" void vh2h2_(double *rCOM, double *bondLength1, double *bondLength2, double *theta1, double *theta2, double *dihedralAngle, double *potl);

//==================================================================================================================================

/*For Matrix diagonalization (Real Symmetric)*/
//----------------------------------------------------------------------------------------------------------------------------------

extern "C" void dsyev_( char* jobz, char* uplo, int* n, double* a, int* lda, double* w, double* work, int* lwork, int* info );

void print_matrix( char* desc, int m, int n, double* a, int lda );
void matdiagsymmetric(int nSizeTotal, double *a);

//==================================================================================================================================

/*For Matrix diagonalization (Complex Hermitian)*/
//----------------------------------------------------------------------------------------------------------------------------------
/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* ZHEEV prototype */
extern "C" void zheev_( char* jobz, char* uplo, int* n, dcomplex* a, int* lda, double* w, dcomplex* work, int* lwork, double* rwork, int* info );
/* Auxiliary routines prototypes */
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );
void print_rmatrix( char* desc, int m, int n, double* a, int lda );

void matdiaghermitian(int nSizeTotal, double *c, double *d);
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );

//==================================================================================================================================


void thetagrid(int nSize, double *cosTheta, double *weightsTheta);
void phigrid(int nSize, double *phi, double &weightsPhi);

double normalizedPlm(int j,int m, double x);
double normalizedPl0(int j,double x);
void normalizationCheck(int nSize, int nSizeRot);

int sizeloop(int nSizeRot, int nSkip);
double rotEnergy(int l);

#endif

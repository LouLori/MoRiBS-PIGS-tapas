#ifndef _MC_ESTIM_H
#define _MC_ESTIM_H 1

#include <fstream>

void  InitMCEstims(void);
void  DoneMCEstims(void);
void  ResetMCEstims(void);

void SaveDensities1D(const char [], double);

void IOxyzAng(int,const char []);

void SaveDensities2D(const char [], double,int); //added by Hui Li
void SaveGraSum(const char fname [], double acount); //added by Toby Zeng
void SaveGxyzSum(const char fname [], double acount); //added by Toby Zeng
void SaveRho1D(const char fname [], double acount, int mode); //added by Toby Zeng
void SaveDensities3D(const char [], double,int); //added by Toby Zeng
void SaveRhoThetaChi(const char [], double,int); //added by Toby Zeng

double GetPotEnergy_Diff(void); //added by Hui Li

double GetPotEnergy_Densities(void);
double GetPotEnergyPIGS(void);
double GetPotEnergyCage(const double *);
double GetTotalEnergy(void);
#ifdef DDCORR
void GetDipoleCorrelationPIMC(double *, double *, double *, double *, double *);
void GetDipoleCorrelationPIGSENT(double *, double *, double *, double *, double *);
void GetDipoleCorrelationPIGS(double *, double *, double *, double *, double *);
#endif
void GetCosThetaPIGS(double &, double *);
void GetCosThetaPIMC(double &, double *);
double GetPhi();
double GetPotEnergy(void);
double GetKinEnergy(void);
double GetPotEnergy_Entanglement(int atom0, int atom1);
double GetEstimNM(void);
double GetEstimDM(void);
double *GetCosThetaEntanglement();
double *GetPhiEntanglement();
double *GetProdUvec12();

double GetConfPoten_Densities(void); // HA test

// super densities

void GetExchangeLength(void);
void GetAreaEstimators(void);
void GetAreaEstim3D(int);

extern double _areas3DMFF[6];  // block accumulated area tensor estimator in dopant-fixed frame
extern double _inert3DMFF[9]; // block accumulated classical moment of inertia in dopant-fixed frame

extern double _areas3DSFF[6];  // block accumulated area tensor estimator in dopant-fixed frame
extern double _inert3DSFF[9]; // block accumulated classical moment of inertia in dopant-fixed frame

void SaveExchangeLength (const char [], double, long int);
void SaveAreaEstimators (const char [], double, long int);
void SaveAreaEstim3D (const char [], double, long int,int);

// permutation sampling
void GetPermutation(void);

// Rotations

double GetRotEnergy(void);
double GetRotEnergyPIGS(void);
double GetRotEnergyCage(void);
double GetRotPlanarEnergy(void);
double GetRotE3D(void); // get rotational energy for nonlinear rotor, added by toby
double GetRotE3Dstep(int, int); // real step in loop of GetRotE3D, added by toby
double GetExactRotEnergy(int);

void GetRCF(void);
void SaveRCF(const char [],double,int);

void GetExactRCF(double *,int,int);

extern int  * _pflags;

extern int PrintXYZprl; // for printing instantaneous xyz coordinates and permutation table for each bloc

extern double ErotSQ; // asymmetric top rotational energy square estimator
extern double Erot_termSQ; // sum of square of each bead's rotational energy estimator

// MPI rotational variables
extern double srotchunk; // chunk summation of rotational energy
extern double srotsum; // total summation of rotational energy
//Last two lines added by Tapas Sahoo
void VectorNormalisation(double *);
double DotProduct(double *, double *);
void CrossProduct(double *, double *, double *);
void UnitVectors(const double *, double *);
double PotFunc(int , int , const double *, const double *, int );

double GetPotEnergyPIGSENT(void);
double GetTotalEnergyPIGSENT(void);
void GetCosThetaPIGSENT(double &, double *);
#ifdef HISTOGRAM
void GetDensities(void);
#endif
#ifdef NEWDENSITY
void GetDensities(void);
void GetDensitiesEndBeads(void);
#endif
#endif  // mc_estim.h

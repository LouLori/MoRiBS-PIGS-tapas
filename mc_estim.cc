// estimators for MC

#include <math.h>
#include <iomanip>

//#include <gsl/gsl_sf_legendre.h>

#include "mc_setup.h"
#include "mc_poten.h"
#include "mc_utils.h"
#include "mc_estim.h"
#include "mc_input.h"
#include "mc_qworm.h"
#include "mc_confg.h"
#include "mc_const.h"

#include <omp.h>

// -----  density distributions ----------------------------

int NUMB_DENS1D;                // number of 1D density distributions
int NUMB_DENS2D;                // number of 2D density distributions
int NUMB_DENS3D;                // number of 3D density distributions added by Toby

const int MC_BINSR = 300;       // number of bins for radius original value=500
const int MC_BINST = 50;       // number of bins for theta original value=100
const int MC_BINSC = 100;       // number of bins for chi, added by Toby original value=200

#ifndef NEWDENSITY
const double MIN_RADIUS = 0.0;  // [AA], min radius for radial distribution functions 
const double MAX_RADIUS = 15.0; // [AA], max radius for radial distribution functions original value=25
#else
const double MIN_RADIUS = -5.0;  // [AA], min radius for radial distribution functions 
const double MAX_RADIUS = 5.0; // [AA], max radius for radial distribution functions original value=25
#endif

double _delta_radius;           // \delta r for density distributions    
double _delta_theta;            // \delta\theta density distributions 
double _delta_chi;              // \delta\chi   density distributions

double _max_radius;             // max radius for density distributions
double _min_radius;             // min radius for density distributions 

#ifdef HISTOGRAM
const int MC_BINSXYZ = 100;       // number of bins for radius original value=500
double _delta_xyz;           // \delta r for density distributions    
double _max_xyz;             // max radius for density distributions
double _min_xyz;             // min radius for density distributions 
double **   _gxyz1D;              // 1D radial distribution functions
double **   _gxyz1D_sum;          // 1D radial distribution functions accumulated
#endif
double **   _gr1D;              // 1D radial distribution functions
double **   _gr1D_sum;          // 1D radial distribution functions accumulated
double ***  _gr2D;              // 2D distribution functions
double ***  _gr2D_sum;          // 2D distribution functions accumulated

double **   _gr3D;              // 3D distribution functions, added by Toby
double **   _gr3D_sum;          // 3D distribution functions accumulated, added by Toby

double *    _relthe_sum;        // Binning of relative euler angle theta between two rotational imaginary slices
double *    _relphi_sum;        // Binning of relative euler angle phi between two rotational imaginary slices
double *    _relchi_sum;        // Binning of relative euler angle chi between two rotational imaginary slices

double *    _drgid;             // radial  grid for distribution functions
double *    _trgid;             // angular grid for distribution functions

void densities_malloc(void);
void densities_init  (void);
void densities_mfree (void);
void densities_reset (int);  //revised by Hui Li

void bin_1Ddensity(double,int);
#ifdef HISTOGRAM
void binxyz_1Ddensity(double,int);
#endif
void bin_2Ddensity(double,double,int);
void bin_3Ddensity(double,double,double,int);

// -----  Rotational correlation functions ----

double ** _rcf;       
double ** _rcf_sum;

#ifdef ROTCORR
const int MAXNATOMS = 20;
const int MAXROTBEADS = 256;
double ** _rcfx;
double ** _rcfx_sum;
double ** _rcfy;
double ** _rcfy_sum;
double  _rcfijx[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
double  _rcfijx_sum[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
double  _rcfijy[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
double  _rcfijy_sum[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
double  _rcfijz[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
double  _rcfijz_sum[MAXNATOMS][MAXNATOMS][NUMB_RCF][MAXROTBEADS];
#endif

void rcf_malloc(void);
void rcf_init  (void);
void rcf_mfree (void);
void rcf_reset (int);

// SUPER (area estimator)

double _areas[2];
double _area2[2];

double _inert[2];

const int PERP = 0;
const int PARL = 1;

// SUPER (area estimator) for non-linear dopant
double _areas3DMFF[6];   // block accumulated area tensor estimator in dopant-fixed frame
double _inert3DMFF[9];  // block accumulated classical moment of inertia in dopant-fixed frame

double _areas3DSFF[6];   // block accumulated area tensor estimator in space-fixed frame
double _inert3DSFF[9];  // block accumulated classical moment of inertia in space-fixed frame

void super_reset (void);

// integer flag
int PrintXYZprl; // for printing instantaneous xyz coordinates and permutation table for each block

// rotational energy square estimator
double ErotSQ; // asymmetric top rotational energy square estimator
double Erot_termSQ; // sum of square of each bead's rotational energy estimator


//----------------------------------------------------

int     * _pflags;
double  * _ploops;

extern "C" void rotden_(double *Eulan1,double *Eulan2,double *Eulrel,double *rho,double *erot,double *esq,double *rhoprp, double *erotpr, double *erotsq, int *istop); // external fortran by toby

extern "C" void vcord_(double *Eulang, double *RCOM, double *Rpt, double *vtable, int *Rgrd,int *THgrd, int *CHgrd, double *Rvmax, double *Rvmin, double *Rvstep, double *vpot3d, double *radret, double *theret, double *chiret, double *hatx, double *haty, double *hatz, int *ivcord);

extern "C" void rsrot_(double *Eulan1,double *Eulan2,double *X_Rot,double *Y_Rot,double *Z_Rot,double *tau,int *iodevn,double *eoff,double *rho,double *erot);

extern "C" void rsline_(double *X_Rot,double *p0,double *tau,double *rho,double *erot);

extern "C" void vspher_(double *r,double *v);

extern "C" void caleng_(double *com_1, double *com_2, double *E_2H2O, double *Eulang_1, double *Eulang_2);
//vh2h2 ---> potential H2-H2: added by Tapas Sahoo

extern "C" void vh2h2_(double *rd, double *r1, double *r2, double *t1, double *t2, double *phi, double *potl);
extern "C" void cluster_(double *com_1, double *com_2, double *Eulang_1, double *Eulang_2, double *E_12);
extern "C" double plgndr(int l, int m, double x);
#ifdef MOLECULEINCAGE
//Pot of H2O-C60 one cage
extern "C" void calengy_(double *com_1, double *Eulang_1, double *E_H2OC60);
//Potential of interaction of two cages with water inside
extern "C" void cluster_(double *com_1, double *com_2, double *Eulang_1, double *Eulang_2, double *E_12);
#endif

void InitMCEstims(void)
{
   densities_init();   
   densities_malloc();
   densities_reset (MC_TOTAL);  //revised by Hui 

   if (BSTYPE >= 0)
   {
       _pflags = new int    [MCAtom[BSTYPE].numb];
       _ploops = new double [MCAtom[BSTYPE].numb];
   }

   if (ROTATION)
   {
      rcf_init();
      rcf_malloc();
      rcf_reset(MC_TOTAL);
   }

// SAVE RESULTS OF PREVIOUS SIMULATIONS

   string fname;

// energy

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_ENG; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// area estimators

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_SUP; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// 3d area estimators in dopant-fixed frame

   fname  = MCFileName.c_str();
   fname += IO_EXT_MFFSUP3D;

   if (FileExist(fname.c_str()))
   {
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// 3d area estimators in space-fixed frame

   fname  = MCFileName.c_str();
   fname += IO_EXT_SFFSUP3D;

   if (FileExist(fname.c_str()))
   {
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }

// exchange length

   fname  = MCFileName.c_str(); 
   fname += IO_EXT_PRL; 

   if (FileExist(fname.c_str()))
   { 
      IOFileBackUp(fname.c_str());
      IOFileDelete(fname.c_str());
   }
}

void DoneMCEstims(void)
{
   densities_mfree();

   if (BSTYPE >= 0)
   {
       delete [] _pflags;
       delete [] _ploops;
   }

   if (ROTATION)
   rcf_mfree();
}

void ResetMCEstims(void)  // reset BLOCK averages
{
   densities_reset(MC_BLOCK);  //revised by Hui

   if (BOSONS)
   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   {
      _pflags [atom] = 0;
      _ploops [atom] = 0.0;
   }

   if (ROTATION)
   rcf_reset(MC_BLOCK);

   super_reset();
}

// -----  BEGIN densities ------------------------------

void densities_init(void)
//
// atom-atom      (the same type) and 
// atom-molecule  distribution functions  
//
// LIMITATION:    the density type corresponds to the atom type
//                no cross densities  
{
   const char *_proc_=__func__;    // "densities_init()";

   if ((NUMB_ATOMTYPES > 1) || (NUMB_MOLCTYPES > 1))
   nrerror(_proc_,"No more then one atom/molecule type: densities and potential energy");   

// if ((NUMB_MOLCTYPES>0) && (NUMB_MOLCS!=1))
// nrerror(_proc_,"Only one molecular impurity");

#ifndef NEWDENSITY
   NUMB_DENS1D = NUMB_ATOMTYPES;  //# atom-atom densities [no cross-distributions]
#else
   NUMB_DENS1D = NDIM;  //# atom-atom densities [no cross-distributions]
#endif

   if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))                          
     NUMB_DENS2D = NUMB_ATOMTYPES;  //# molecule-atoms distributions 

   // Toby adds 3d distribution
   if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))                          
     NUMB_DENS3D = NUMB_ATOMTYPES+NUMB_MOLCTYPES;  //# non-linear molecule-atoms distributions 

  _max_radius   = MAX_RADIUS/Units.length;  // max radius for radial distributions;
  _min_radius   = MIN_RADIUS/Units.length;  // min radius for radial distributions;

  _delta_radius = (_max_radius - _min_radius)/(double)MC_BINSR;
  _delta_theta  =  M_PI/(double)(MC_BINST - 1);
// Toby adds chi bin
  _delta_chi    =  2.0*M_PI/(double)(MC_BINSC - 1);
#ifdef HISTOGRAM
  _max_xyz   = 1.0;
  _min_radius  = -1.0;
  _delta_xyz = (_max_xyz - _min_xyz)/(double)MC_BINSXYZ;
#endif
}

void densities_malloc(void)
// memory allocation for densities
{
	_gr1D  = doubleMatrix(NUMB_DENS1D,MC_BINSR);
	_gr1D_sum  = doubleMatrix(NUMB_DENS1D,MC_BINSR);
#ifdef HISTOGRAM
	_gxyz1D  = doubleMatrix(NDIM,MC_BINSXYZ);
	_gxyz1D_sum  = doubleMatrix(NDIM,MC_BINSXYZ);
#endif

	if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))
    {
    	_gr2D  = new double ** [NUMB_DENS2D];
      	_gr2D_sum  = new double ** [NUMB_DENS2D];  //added by Hui Li

      	for (int id=0;id<NUMB_DENS2D;id++) 
		{
	  		_gr2D[id] = doubleMatrix(MC_BINSR,MC_BINST);
	  		_gr2D_sum[id] = doubleMatrix(MC_BINSR,MC_BINST); //added by Hui Li
		}
    }

  	if (IMPURITY && (MCAtom[IMTYPE].molecule == 3))
    {
      	_gr2D  = new double ** [NUMB_DENS2D];
      	_gr2D_sum  = new double ** [NUMB_DENS2D];  //added by Hui Li

      	for (int id=0;id<NUMB_DENS2D;id++) 
		{
	  		_gr2D[id] = doubleMatrix(MC_BINSR,MC_BINST);
	  		_gr2D_sum[id] = doubleMatrix(MC_BINSR,MC_BINST); //added by Hui Li
		}
    }
// the following if block added by Toby
  	if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))
  	{
     	_gr3D  = doubleMatrix(NUMB_DENS3D,MC_BINSR*MC_BINST*MC_BINSC);
     	_gr3D_sum  = doubleMatrix(NUMB_DENS3D,MC_BINSR*MC_BINST*MC_BINSC);

     	_relthe_sum = new double [MC_BINST];
     	_relphi_sum = new double [MC_BINSC];
     	_relchi_sum = new double [MC_BINSC];
  	}
}

void densities_mfree(void)
{
	free_doubleMatrix(_gr1D);
   	free_doubleMatrix(_gr1D_sum);
#ifdef HISTOGRAM
	free_doubleMatrix(_gxyz1D);
   	free_doubleMatrix(_gxyz1D_sum);
#endif
   
   	if (IMPURITY && MCAtom[IMTYPE].molecule == 1)
   	{
      	for (int id=0;id<NUMB_DENS2D;id++) 
      	{
      		free_doubleMatrix(_gr2D[id]);
      		free_doubleMatrix(_gr2D_sum[id]);    //added by Hui Li
      	}
      	delete [] _gr2D;
      	delete [] _gr2D_sum;   //added by Hui Li
    }
   	if (IMPURITY && MCAtom[IMTYPE].molecule == 3)
    {
       	for (int id=0;id<NUMB_DENS2D;id++) 
	 	{
	   		free_doubleMatrix(_gr2D[id]);
	   		free_doubleMatrix(_gr2D_sum[id]);    //added by Hui Li
	 	}
       	delete [] _gr2D;
       	delete [] _gr2D_sum;   //added by Hui Li
    }
   	if (IMPURITY && MCAtom[IMTYPE].molecule == 2)
   	{
      	free_doubleMatrix(_gr3D);
      	free_doubleMatrix(_gr3D_sum);
   	}
}

void densities_reset(int mode)     //revised by Hui Li
{
#ifdef DEBUG_PIMC
   	const char *_proc_= __func__; //  rcf_reset() 

   	if ((mode != MC_BLOCK) && (mode != MC_TOTAL))
	nrerror(_proc_,"Unknow mode");
#endif

   	for (int id=0;id<NUMB_DENS1D;id++) 
    for (int ir=0;ir<MC_BINSR;ir++) 
    {
	 	if(mode == MC_BLOCK)
	   	_gr1D[id][ir] = 0.0;

	 	if(mode == MC_TOTAL)
	   	_gr1D_sum[id][ir] = 0.0;
    }
#ifdef HISTOGRAM
   	for (int id=0;id<NDIM;id++) 
    for (int ixyz=0;ixyz<MC_BINSXYZ;ixyz++) 
    {
	 	if(mode == MC_BLOCK)
	   	_gxyz1D[id][ixyz] = 0.0;

	 	if(mode == MC_TOTAL)
	   	_gxyz1D_sum[id][ixyz] = 0.0;
    }
#endif

   	if (IMPURITY && (MCAtom[IMTYPE].molecule == 1))
	for (int id=0;id<NUMB_DENS2D;id++) 
    for (int ir=0;ir<MC_BINSR;ir++) 
	for (int it=0;it<MC_BINST;it++) 
	{
		if(mode == MC_BLOCK)
	    _gr2D[id][ir][it] = 0.0;     // block average

	    if(mode == MC_TOTAL)        
	    _gr2D_sum[id][ir][it] = 0.0;     // accumulated average
	}

   	if (IMPURITY && (MCAtom[IMTYPE].molecule == 3))
    for (int id=0;id<NUMB_DENS2D;id++) 
    for (int ir=0;ir<MC_BINSR;ir++) 
	for (int it=0;it<MC_BINST;it++) 
	{
	    if(mode == MC_BLOCK)
	    _gr2D[id][ir][it] = 0.0;     // block average

	    if(mode == MC_TOTAL)        
	    _gr2D_sum[id][ir][it] = 0.0;     // accumulated average
	}
    
   	if (IMPURITY && (MCAtom[IMTYPE].molecule == 2))
	{
       	for (int id=0;id<NUMB_DENS3D;id++)
	 	for (int ir=0;ir<MC_BINSR;ir++)
	   	for (int it=0;it<MC_BINST;it++)
	    for (int ic=0;ic<MC_BINSC;ic++)
	    {

			int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;

			if(mode == MC_BLOCK)
			_gr3D[id][ijk] = 0.0;     // block average

			if(mode == MC_TOTAL)
			_gr3D_sum[id][ijk] = 0.0;     // accumulated average

	    }

       	if(mode == MC_TOTAL)
	 	{

	   		for (int it=0;it<MC_BINST;it++)
	     	_relthe_sum[it]=0.0;

	   		for (int ic=0;ic<MC_BINSC;ic++)
	     	{
	       		_relphi_sum[ic]=0.0;
	       		_relchi_sum[ic]=0.0;
	     	}

	 	}

	}
}

//------- RCF -------------------

void rcf_init(void)
{
}

void rcf_malloc(void)
// memory allocation for rotational correlation functions
{
  _rcf     = doubleMatrix(NUMB_RCF,NumbRotTimes);
  _rcf_sum = doubleMatrix(NUMB_RCF,NumbRotTimes);
#ifdef ROTCORR
  _rcfx     = doubleMatrix(NUMB_RCF,NumbRotTimes);
  _rcfx_sum = doubleMatrix(NUMB_RCF,NumbRotTimes);
  _rcfy     = doubleMatrix(NUMB_RCF,NumbRotTimes);
  _rcfy_sum = doubleMatrix(NUMB_RCF,NumbRotTimes);
#endif
}

void rcf_mfree(void)
{
   free_doubleMatrix(_rcf);
   free_doubleMatrix(_rcf_sum);
#ifdef ROTCORR
   free_doubleMatrix(_rcfx);
   free_doubleMatrix(_rcfx_sum);
   free_doubleMatrix(_rcfy);
   free_doubleMatrix(_rcfy_sum);
#endif
}

void rcf_reset(int mode)
{
#ifdef DEBUG_PIMC
   	const char *_proc_= __func__; //  rcf_reset() 

   	if ((mode != MC_BLOCK) && (mode != MC_TOTAL))
   	nrerror(_proc_,"Unknow mode"); 
#endif 

   	for (int ip=0;ip<NUMB_RCF; ip++) 
	{
   		for (int it=0;it<NumbRotTimes;it++)
   		{ 
      		_rcf[ip][it]  = 0.0;  // block average
#ifdef ROTCORR
      		_rcfx[ip][it] = 0.0;  // block average
      		_rcfy[ip][it] = 0.0;  // block average
   			for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
			{ 
   				for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
    			{
      				_rcfijx[atom0][atom1][ip][it] = 0.0;
      				_rcfijy[atom0][atom1][ip][it] = 0.0;
      				_rcfijz[atom0][atom1][ip][it] = 0.0;
    			}
			}
#endif

       		if (mode == MC_TOTAL)   // total average
			{
      			_rcf_sum[ip][it] = 0.0;
#ifdef ROTCORR
     			_rcfx_sum[ip][it] = 0.0;
       			_rcfy_sum[ip][it] = 0.0;
   				for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
				{
   			    	for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
    				{
      					_rcfijx_sum[atom0][atom1][ip][it] = 0.0;
      					_rcfijy_sum[atom0][atom1][ip][it] = 0.0;
      					_rcfijz_sum[atom0][atom1][ip][it] = 0.0;
    				}
				}
#endif
       		}
    	}
	}
}
/*
//added by Hui Li
double GetPotEnergy_Diff(void)
{
   const char *_proc_=__func__; //  GetPotEnergy_Diff()  

#ifdef DEBUG_WORM
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif

   double dr[NDIM];
   double dspot = 0.0;

   for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
   for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
   {
      int type0   = MCType[atom0];
      int type1   = MCType[atom1];

      int offset0 = NumbTimes*atom0;
      int offset1 = NumbTimes*atom1;

      for (int it=0;it<NumbTimes;it++) 	    
      {  
         int t0 = offset0 + it;
         int t1 = offset1 + it;

         double dr2 = 0.0;  		 
         for (int id=0;id<NDIM;id++)
         {
            dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

            if (MINIMAGE)
            dr[id] -= (BoxSize*rint(dr[id]/BoxSize));

            dr2    += (dr[id]*dr[id]);
         }
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
         double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

         if (MCAtom[type0].molecule||MCAtom[type1].molecule)  // 2D interaction 
         {
         //  type 1 is a molecule 
     
             int sgn   = 1;             // set to -1 to correct the orientaion of dr

             int tm    = offset1 + it/RotRatio;
//           int tm    = offset1 + floor((double)it/(double)RotRatio);

             int typep = type1;         // define the type of the potential
             int typed = type0;         // define the type of the density

         //  type 0 is a molecule ?   
             if (MCAtom[type0].molecule)  // does not work for two molecules
             {
                sgn   = -1;   
                tm    = offset0 + it/RotRatio;
                typep = type0; 
                typed = type1; 
             }

             double cost = 0.0;
             for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             cost /= r;                     // cos(theta)
             cost *= sgn;                   // correct the orientation 

             dspot += DLPot2D(r,cost,typep);  // potential differencies  
         }
      }  // LOOP OVER TIME SLICES
   }     // LOOP OVER ATOM PAIRS

   return (dspot/(double)NumbTimes);
}
*/

double GetPotEnergyPIGS(void)
{
	const char *_proc_=__func__; 

#ifdef DEBUG_WORM
	if (Worm.exists)
	nrerror(_proc_," Only for Z-configurations");
#endif

    string stype = MCAtom[IMTYPE].type;
   	int it = ((NumbRotTimes - 1)/2);
	double spot;
	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb > 1) )
	{
		double Eulang0[NDIM], Eulang1[NDIM];

        spot = 0.0;
        for (int atom0 = 0; atom0 < (NumbAtoms-1); atom0++)
		{
           	int offset0 = NumbTimes*atom0;
           	int t0 = offset0 + it;

			if (stype == HF)
			{
   				Eulang0[PHI] = MCAngles[PHI][t0];
   				Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   				Eulang0[CHI] = 0.0;
			}

        	for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
       		{
            	int offset1 = NumbTimes*atom1;
            	int t1 = offset1 + it;

            	if (stype == HF)
            	{
   					Eulang1[PHI] = MCAngles[PHI][t1];
   					Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   					Eulang1[CHI] = 0.0;
                	spot += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            	} //stype

            	if (stype == H2)
           		{
                	double s1 = 0.0;
                	double s2 = 0.0;
                	double dr2 = 0.0;
					double dr[NDIM];
                	for (int id=0;id<NDIM;id++)
                	{
                    	dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                    	dr2    += (dr[id]*dr[id]);
                    	double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t0];
                    	double cst2 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t1];
                    	s1 += cst1;
                    	s2 += cst2;
                	}
                	double r = sqrt(dr2);
                	double th1 = acos(s1/r);
                	double th2 = acos(s2/r);

                	double b1[NDIM];
                	double b2[NDIM];
                	double b3[NDIM];
                	for (int id=0;id<NDIM;id++)
                	{
                    	b1[id] = MCCosine[id][t0];
                    	b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                    	b3[id] = MCCosine[id][t1];
                	}
                	VectorNormalisation(b1);
                	VectorNormalisation(b2);
                	VectorNormalisation(b3);

                	//Calculation of dihedral angle 
                	double n1[NDIM];
                	double n2[NDIM];
                	double mm[NDIM];

                	CrossProduct(b2, b1, n1);
                	CrossProduct(b2, b3, n2);
                	CrossProduct(b2, n2, mm);

                	double xx = DotProduct(n1, n2);
                	double yy = DotProduct(n1, mm);

                	double phi = atan2(yy, xx);
                	if (phi<0.0) phi += 2.0*M_PI;

                	//Dihedral angle calculation is completed here
                	double r1 = 0.74;// bond length in Angstrom
					r1 /= BOHRRADIUS;
                	double r2 = r1;// bond length in bohr
                	double rd = r/BOHRRADIUS;
                	double potl;
                	vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                	spot += potl*CMRECIP2KL;
            	} //stype

        	}// loop over atoms1 
       	}// loop over atoms0 
    }
    if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
    {
        int offset0 = 0;
#ifndef GAUSSIANMOVE
        int t0  = offset0 + it;
        double E12     = -2.0*DipoleMomentAU2*MCCosine[2][t0]/(RR*RR*RR);
        spot    = E12*AuToKelvin;
#else
        int t0  = offset0 + (NumbTimes-1)/2;
        double spot3d = 0.0;
        for (int id = 0; id < NDIM; id++)
        {
            spot3d += 0.5*MCCoords[id][t0]*MCCoords[id][t0];
        }
        spot   = spot3d;
#endif
    }

    double spot_cage = 0.0;
#ifdef CAGEPOT
    for (int atom0 = atomStart; atom0 < NumbAtoms; atom0++)
	{
        int offset0 = NumbTimes*atom0;
        int t0 = offset0 + it;
    	double cost = MCAngles[CTH][t0];
    	double phi = MCAngles[PHI][t0];
    	if (phi < 0.0) phi = 2.0*M_PI + phi;
    	phi = fmod(phi,2.0*M_PI);
		int type0   =  MCType[atom0];
    	spot_cage += LPot2DRotDOF(cost,phi,type0);
	}
#endif
	double spotReturn = (spot + spot_cage);
	return spotReturn;
}

double GetPotEnergyPIGSENT(void)
{
	const char *_proc_=__func__; 

#ifdef DEBUG_WORM
	if (Worm.exists)
	nrerror(_proc_," Only for Z-configurations");
#endif

    string stype = MCAtom[IMTYPE].type;
   	int it = ((NumbRotTimes - 1)/2);

	double spot_sector, spot_pair, spot_cage, spot_sector_cage;
	spot_sector      = 0.0;
	spot_sector_cage = 0.0;
	int atomStart, atomEnd;
	for (int isector = 0; isector < 2; isector++)
	{
		if (isector == 0)
		{
			atomStart = 0;
			atomEnd   = NumbAtoms/2;
		}
		if (isector == 1)
		{
			atomStart = NumbAtoms/2;
			atomEnd   = NumbAtoms;
		}

		if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb > 1) )
		{
			double Eulang0[NDIM], Eulang1[NDIM];

        	spot_pair = 0.0;
        	for (int atom0=atomStart;atom0<(atomEnd-1);atom0++)
			{
           		int offset0 = NumbTimes*atom0;
           		int t0 = offset0 + it;

				if (stype == HF)
				{
   					Eulang0[PHI] = MCAngles[PHI][t0];
   					Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   					Eulang0[CHI] = 0.0;
				}

        		for (int atom1=(atom0+1);atom1<atomEnd;atom1++)
       			{
            		int offset1 = NumbTimes*atom1;
            		int t1 = offset1 + it;

            		if (stype == HF)
            		{
   						Eulang1[PHI] = MCAngles[PHI][t1];
   						Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   						Eulang1[CHI] = 0.0;
                		spot_pair += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            		} //stype

            		if (stype == H2)
           			{
                		double s1 = 0.0;
                		double s2 = 0.0;
                		double dr2 = 0.0;
						double dr[NDIM];
                		for (int id=0;id<NDIM;id++)
                		{
                    		dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                    		dr2    += (dr[id]*dr[id]);
                    		double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t0];
                    		double cst2 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t1];
                    		s1 += cst1;
                    		s2 += cst2;
                		}
                		double r = sqrt(dr2);
                		double th1 = acos(s1/r);
                		double th2 = acos(s2/r);

                		double b1[NDIM];
                		double b2[NDIM];
                		double b3[NDIM];
                		for (int id=0;id<NDIM;id++)
                		{
                    		b1[id] = MCCosine[id][t0];
                    		b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                    		b3[id] = MCCosine[id][t1];
                		}
                		VectorNormalisation(b1);
                		VectorNormalisation(b2);
                		VectorNormalisation(b3);

                		//Calculation of dihedral angle 
                		double n1[NDIM];
                		double n2[NDIM];
                		double mm[NDIM];

                		CrossProduct(b2, b1, n1);
                		CrossProduct(b2, b3, n2);
                		CrossProduct(b2, n2, mm);

                		double xx = DotProduct(n1, n2);
                		double yy = DotProduct(n1, mm);

                		double phi = atan2(yy, xx);
                		if (phi<0.0) phi += 2.0*M_PI;

                		//Dihedral angle calculation is completed here
                		double r1 = 0.74;// bond length in Angstrom
						r1 /= BOHRRADIUS;
                		double r2 = r1;// bond length in bohr
                		double rd = r/BOHRRADIUS;
                		double potl;
                		vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                		spot_pair += potl*CMRECIP2KL;
            		} //stype

        		}// loop over atoms1 
       		}// loop over atoms0 
			spot_sector += spot_pair;

    		spot_cage = 0.0;
#ifdef CAGEPOT
    		for (int atom0 = atomStart; atom0 < atomEnd; atom0++)
			{
        		int offset0 = NumbTimes*atom0;
        		int t0 = offset0 + it;
    			double cost = MCAngles[CTH][t0];
    			double phi = MCAngles[PHI][t0];
    			if (phi < 0.0) phi = 2.0*M_PI + phi;
    			phi = fmod(phi,2.0*M_PI);
				int type0   =  MCType[atom0];
    			spot_cage += LPot2DRotDOF(cost,phi,type0);
			}
			spot_sector_cage += spot_cage;
#endif
    	}
	}
	double spotReturn = 0.5*(spot_sector + spot_sector_cage);
	return spotReturn;
}

double GetPotEnergy_Densities(void)
// should be compatible with PotEnergy() from mc_piqmc.cc
{
	const char *_proc_=__func__; //  GetPotEnergy_Densities()  

#ifdef DEBUG_WORM
	if (Worm.exists)
	nrerror(_proc_," Only for Z-configurations");
#endif

	// double dr[NDIM];
    string stype = MCAtom[IMTYPE].type;
	double spot = 0.0;
	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb > 1) )
	{
        spot = 0.0;
        for (int atom0 = 0; atom0 < (NumbAtoms-1); atom0++)
		{
           	int offset0 = NumbTimes*atom0;
        	for (int atom1 = (atom0+1); atom1 < NumbAtoms; atom1++)
        	{
            	int offset1 = NumbTimes*atom1;

		    	double spot_pair = 0.0;
#pragma omp parallel for reduction(+: spot_pair)
		    	for (int it = 0; it < NumbTimes; it++) 	    
		    	{  
                	int t0 = offset0 + it;
                	int t1 = offset1 + it;

                	if (stype == H2)
                	{
                    	double s1 = 0.0;
                    	double s2 = 0.0;
                    	double dr2 = 0.0;
				    	double dr[NDIM];
                    	for (int id=0;id<NDIM;id++)
                    	{
                        	dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                        	dr2    += (dr[id]*dr[id]);
                        	double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t0];
                        	double cst2 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][t1];
                        	s1 += cst1;
                        	s2 += cst2;
                    	}
                    	double r = sqrt(dr2);
                    	double th1 = acos(s1/r);
                    	double th2 = acos(s2/r);

                    	double b1[NDIM];
                    	double b2[NDIM];
                    	double b3[NDIM];
                    	for (int id = 0; id < NDIM; id++)
                    	{
                        	b1[id] = MCCosine[id][t0];
                        	b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                        	b3[id] = MCCosine[id][t1];
                    	}
                    	VectorNormalisation(b1);
                    	VectorNormalisation(b2);
                    	VectorNormalisation(b3);

                    	//Calculation of dihedral angle 
                    	double n1[NDIM];
                    	double n2[NDIM];
                    	double mm[NDIM];

                    	CrossProduct(b2, b1, n1);
                    	CrossProduct(b2, b3, n2);
                    	CrossProduct(b2, n2, mm);

                    	double xx = DotProduct(n1, n2);
                    	double yy = DotProduct(n1, mm);

                    	double phi = atan2(yy, xx);
                    	if (phi<0.0) phi += 2.0*M_PI;

                    	//Dihedral angle calculation is completed here
                		double r1 = 0.74;// bond length in Angstrom
						r1 /= BOHRRADIUS;
                    	double r2 = r1;// bond length in bohr
                    	double rd = r/BOHRRADIUS;
                    	double potl;
                    	vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                    	spot_pair += potl*CMRECIP2KL;
                	} //stype

                	if (stype == HF)
                	{
						double Eulang0[NDIM], Eulang1[NDIM];
   						Eulang0[PHI] = MCAngles[PHI][t0];
   						Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   						Eulang0[CHI] = 0.0;
   						Eulang1[PHI] = MCAngles[PHI][t1];
   						Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   						Eulang1[CHI] = 0.0;
                		spot_pair += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
                	} //stype
            	}
            	spot += spot_pair;
        	}// loop over atoms (molecules)
        }// loop over atoms (molecules)
    }
    if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
    {
        spot = 0.0;
        double dm   = DipoleMoment/AuToDebye;
        double dm2  = dm*dm;
        double RR = Distance/BOHRRADIUS;

        int offset0 = 0;

        double spot_pair = 0.0;
#pragma omp parallel for reduction(+: spot_pair)
	    for (int it = 0; it < NumbTimes; it++) 	  
        {  
            int t0  = offset0 + it;

            double E12 = -2.0*dm2*MCCosine[2][t0]/(RR*RR*RR);
            spot_pair  = E12*AuToKelvin;
#ifdef POTZERO
			spot_pair = 0.0;
#endif
        }
      	spot += spot_pair;
    }
#ifdef IOWRITE
	if ( MCAtom[IMTYPE].numb > 1)
    {
	for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
	for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
	{
		int type0   = MCType[atom0];
		int type1   = MCType[atom1];

		int offset0 = NumbTimes*atom0;
		int offset1 = NumbTimes*atom1;

		double spot_pair=0.0;

		#pragma omp parallel for reduction(+: spot_pair)
		for (int it=0;it<NumbTimes;it++) 	    
		{  
			double dr[NDIM];
			int t0 = offset0 + it;
			int t1 = offset1 + it;

			double dr2 = 0.0;  		 
			for (int id=0;id<NDIM;id++)
			{
				dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

				if (MINIMAGE)
				dr[id] -= (BoxSize*rint(dr[id]/BoxSize));
				dr2    += (dr[id]*dr[id]);
			}
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
			double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

			if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
			{
				//  type 1 is a molecule 
     
				int sgn   = 1;             // set to -1 to correct the orientaion of dr

				int tm    = offset1 + it/RotRatio;
				//int tm    = offset1 + floor((double)it/(double)RotRatio);

				int typep = type1;         // define the type of the potential
				int typed = type0;         // define the type of the density

				//  type 0 is a molecule ?   
				if (MCAtom[type0].molecule == 1)  // does not work for two molecules
				{
					sgn   = -1;   
                	tm    = offset0 + it/RotRatio;
                	typep = type0; 
                	typed = type1; 
             	}

             	double cost = 0.0;
             	for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             	cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             	cost /= r;                     // cos(theta)
             	cost *= sgn;                   // correct the orientation 

             	bin_2Ddensity (r,cost,typed);  // densities 
				spot_pair += LPot2D(r,cost,typep);  // potential energy 
				//cout<<"it="<<it<<" spot_pair="<<spot_pair<<" LPot2D="<<LPot2D(r,cost,typep)<<endl;
			}
//----------- [ATOM - NON-LINEAR MOLECULE] ----------------------
			else if (((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && (MCAtom[type0].molecule != MCAtom[type1].molecule) ) // 3D interaction, no density is calculated now
			{

            	int tm;
            	int typed;

            	double RCOM[3];
            	double Rpt[3];
        	    double Eulang[3];
            	double vpot3d;
            	double radret;
            	double theret;
            	double chiret;
            	double hatx[3];
            	double haty[3];
            	double hatz[3];
            	int    ivcord = 0;
            	if(MCAtom[type0].molecule == 2)
            	{
					//determine type of atoms for bin_3Ddensity
					typed = type1;
               		tm  = offset0 + it/RotRatio;
             		for (int id=0;id<NDIM;id++)
              		{
                  		RCOM[id] = MCCoords[id][t0];
                 		Rpt[id]  = MCCoords[id][t1];
              		}	
            	}
            	else
            	{
				//determine type of atoms for bin_3Ddensity
					typed = type0;
               		tm  = offset1 + it/RotRatio;
               		for (int id=0;id<NDIM;id++)
               		{
                		Rpt[id]  = MCCoords[id][t0];
                		RCOM[id] = MCCoords[id][t1];
               		}
            	}
            	Eulang[PHI]=MCAngles[PHI][tm];
            	Eulang[CTH]=acos(MCAngles[CTH][tm]);
            	Eulang[CHI]=MCAngles[CHI][tm];

            	if( ISPHER == 0)
            	{
               		vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);
            	}
            	else if( ISPHER == 1)
            	{
               		radret = r;
               		vspher_(&radret,&vpot3d);
               		theret = 0.0;
               		chiret = 0.0;
            	}

            	bin_3Ddensity (radret,theret,chiret,typed);  // accumulate density

            	spot_pair += vpot3d;
        	}
//----------[NON-LINEAR - NON-LINEAR from GG]
         	else if ( ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
         	{
          	//if ( (MCType[atom0] == IMTYPE) && (MCType[atom1] == IMTYPE) )
         	//{
         	//if ( (MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2) )
         	//{
        	//cout<<"Le if de GG: MCAtom[type0].molecule MCAtom[type1].molecule"<<MCAtom[type0].molecule<<" "<<MCAtom[type1].molecule<<endl;
             	double com_1[3];
             	double com_2[3];
             	double Eulang_1[3];
             	double Eulang_2[3];
             	double E_2H2O;
             	for (int id=0;id<NDIM;id++)
             	{
                	com_1[id] = MCCoords[id][t0];
                  	com_2[id] = MCCoords[id][t1];
             	}
             	int tm0=offset0 + it/RotRatio;
             	int tm1=offset1 + it/RotRatio;
             	Eulang_1[PHI]=MCAngles[PHI][tm0];
             	Eulang_1[CTH]=acos(MCAngles[CTH][tm0]);
             	Eulang_1[CHI]=MCAngles[CHI][tm0];
             	Eulang_2[PHI]=MCAngles[PHI][tm1];
             	Eulang_2[CTH]=acos(MCAngles[CTH][tm1]);
             	Eulang_2[CHI]=MCAngles[CHI][tm1];
             	caleng_(com_1, com_2, &E_2H2O,
                	    Eulang_1, Eulang_2);
             	spot_pair += E_2H2O;
          		//}  //
          		//}
         	}
//----------[ATOM - ATOM] ------------------------------- 
         	else                // only one atom type    
         	if ((type0 == type1) && MCAtom[type0].molecule == 0) // no "cross" densities 
         	{
            	bin_1Ddensity (r,type1);    // densities 
            	spot_pair += SPot1D(r,type1);    // potential energy
         	}
      	}  // LOOP OVER TIME SLICES
      	spot += spot_pair;
   	}     // LOOP OVER ATOM PAIRS
	}
#endif

    double spot_cage = 0.0;
#ifdef CAGEPOT
    for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
    {
        int offset0 = NumbTimes*atom0;

        double spot_beads=0.0;
        #pragma omp parallel for reduction(+: spot_beads)
        for (int it = 0; it < NumbTimes; it++)
        {
            int t0 = offset0 + it;
            double cost = MCAngles[CTH][t0];
            double phi = MCAngles[PHI][t0];
            if (phi < 0.0) phi = 2.0*M_PI + phi;
            phi = fmod(phi,2.0*M_PI);
            int type0   =  MCType[atom0];
            spot_beads += LPot2DRotDOF(cost,phi,type0);
        }
        spot_cage += spot_beads;
    }
#endif
	double spotReturn = spot + spot_cage;
	return (spotReturn/(double)NumbTimes);
}

#ifdef HISTOGRAM
void GetDensities(void)
{
	const char *_proc_=__func__; //  GetPotEnergy_Densities()  

	for (int atom0 = 0; atom0<1; atom0++)
	{
		for (int it = 0; it < 1; it++)
		{
			int t0 = it + atom0*NumbRotTimes;
			for (int id=0;id<NDIM;id++)
			{
       			binxyz_1Ddensity(MCCosine[id][t0],id);    // densities 
			}
		}
	}
}
#endif

#ifdef NEWDENSITY
void GetDensities(void)
{
	const char *_proc_=__func__; //  GetPotEnergy_Densities()  

    if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
    {
		int atom0 = 0;
		int it = (NumbTimes - 1)/2;
		int t0 = it + atom0*NumbTimes;
		for (int id=0;id<NDIM;id++)
		{
			double r = MCCoords[id][t0]/BOHRRADIUS;
           	bin_1Ddensity (r,id);    // densities 
		}
	}
}

void GetDensitiesEndBeads(void)
{
    const char *_proc_=__func__; //  GetPotEnergy_Densities()

    if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
    {
        int atom0 = 0;
       	int t0 = 0 + atom0*NumbTimes;
       	int tp = (NumbTimes-1) + atom0*NumbTimes;
       	for (int id=0;id<NDIM;id++)
       	{
           	double r = (MCCoords[id][t0] + MCCoords[id][tp])/BOHRRADIUS;
           	bin_1Ddensity (r,id);    // densities
		}
    }
}
#endif

double GetTotalEnergy(void)
{
    string stype = MCAtom[IMTYPE].type;
	double spot;
	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb > 1) )
	{
        spot = 0.0;
        for (int atom0 = 0; atom0 < (NumbAtoms-1); atom0++)
		{
           	int offset0 = NumbTimes*atom0;

        	for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
        	{
            	int offset1 = NumbTimes*atom1;

        		double spot_pair=0.0;
            	#pragma omp parallel for reduction(+: spot_pair)
            	for (int it = 0; it < NumbTimes; it += (NumbTimes - 1))
				{
                	int t0 = offset0 + it;
                	int t1 = offset1 + it;
                	int tm0=offset0 + it/RotRatio;
                	int tm1=offset1 + it/RotRatio;

                	if (stype == H2)
                	{
                    	double s1 = 0.0;
                    	double s2 = 0.0;
                    	double dr2 = 0.0;
				    	double dr[NDIM];
                    	for (int id=0;id<NDIM;id++)
                    	{
                        	dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                        	dr2    += (dr[id]*dr[id]);
                        	double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][tm0];
                        	double cst2 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][tm1];
                        	s1 += cst1;
                        	s2 += cst2;
                    	}
                    	double r = sqrt(dr2);
                    	double th1 = acos(s1/r);
                    	double th2 = acos(s2/r);

                    	double b1[NDIM];
                    	double b2[NDIM];
                    	double b3[NDIM];
                    	for (int id=0;id<NDIM;id++)
                    	{
                        	b1[id] = MCCosine[id][tm0];
                        	b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                        	b3[id] = MCCosine[id][tm1];
                    	}
                    	VectorNormalisation(b1);
                    	VectorNormalisation(b2);
                    	VectorNormalisation(b3);

                    	//Calculation of dihedral angle 
                    	double n1[NDIM];
                    	double n2[NDIM];
                    	double mm[NDIM];

                    	CrossProduct(b2, b1, n1);
                    	CrossProduct(b2, b3, n2);
                    	CrossProduct(b2, n2, mm);

                    	double xx = DotProduct(n1, n2);
                    	double yy = DotProduct(n1, mm);

                    	double phi = atan2(yy, xx);
                    	if (phi<0.0) phi += 2.0*M_PI;

                    	//Dihedral angle calculation is completed here
                		double r1 = 0.74;// bond length in Angstrom
						r1 /= BOHRRADIUS;
                    	double r2 = r1;// bond length in bohr
                    	double rd = r/BOHRRADIUS;
                    	double potl;
                    	vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                    	spot_pair += potl*CMRECIP2KL;
                	} //stype
                	if (stype == HF)
                	{
						double Eulang0[NDIM], Eulang1[NDIM];
   						Eulang0[PHI] = MCAngles[PHI][t0];
   						Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   						Eulang0[CHI] = 0.0;
   						Eulang1[PHI] = MCAngles[PHI][t1];
   						Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   						Eulang1[CHI] = 0.0;
                		spot_pair += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
                	} //stype
				}//loop over beads
				spot += spot_pair;
        	}// loop over atoms (molecules)
        }// loop over atoms (molecules)
    }
    if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
    {
        int offset0 = 0;

        spot = 0.0;
        double E12;
        for (int it = 0; it < NumbTimes; it += (NumbTimes - 1))
		{
            int t0  = offset0 + it;

#ifndef GAUSSIANMOVE
            E12     = -2.0*DipoleMomentAU2*MCCosine[2][t0]/(RR*RR*RR);
            spot   += E12*AuToKelvin;
#else
			double spot3d = 0.0;
			for (int id = 0; id < NDIM; id++)
			{
            	spot3d += 0.5*MCCoords[id][t0]*MCCoords[id][t0];
			}
            spot   += spot3d;
#endif
        }
    }

    double spot_cage = 0.0;
#ifdef CAGEPOT
    for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
    {
        int offset0 = NumbTimes*atom0;

   		double spot_beads=0.0;
       	#pragma omp parallel for reduction(+: spot_beads)
       	for (int it = 0; it < NumbTimes; it += (NumbTimes - 1))
		{
        	int t0 = offset0 + it;
        	double cost = MCAngles[CTH][t0];
        	double phi = MCAngles[PHI][t0];
        	if (phi < 0.0) phi = 2.0*M_PI + phi;
        	phi = fmod(phi,2.0*M_PI);
        	int type0   =  MCType[atom0];
        	spot_beads += LPot2DRotDOF(cost,phi,type0);
		}
		spot_cage += spot_beads;
    }
#endif
	double spotReturn = 0.5*(spot + spot_cage);
	return spotReturn;
}

double GetTotalEnergyPIGSENT(void)
{
    string stype = MCAtom[IMTYPE].type;
	double spot_sector, spot_pair, spot_beads, spot_cage, spot_sector_cage;
	spot_sector = 0.0;
	spot_sector_cage = 0.0;
	int atomStart, atomEnd;
	for (int isector = 0; isector < 2; isector++)
	{
		if (isector == 0)
		{
			atomStart = 0;
			atomEnd   = NumbAtoms/2;
		}
		if (isector == 1)
		{
			atomStart = NumbAtoms/2;
			atomEnd   = NumbAtoms;
		}

		if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb > 1) )
		{
			spot_pair = 0.0;
        	for (int atom0=atomStart;atom0<(atomEnd-1);atom0++)
			{
           		int offset0 = NumbTimes*atom0;

        		for (int atom1=(atom0+1);atom1<atomEnd;atom1++)
        		{
            		int offset1 = NumbTimes*atom1;

        			double spot_beads = 0.0;
            		#pragma omp parallel for reduction(+: spot_beads)
            		for (int it = 0; it < NumbTimes; it += (NumbTimes - 1))
					{
                		int t0 = offset0 + it;
                		int t1 = offset1 + it;
                		int tm0=offset0 + it/RotRatio;
                		int tm1=offset1 + it/RotRatio;

                		if (stype == H2)
                		{
                    		double s1 = 0.0;
                    		double s2 = 0.0;
                    		double dr2 = 0.0;
				    		double dr[NDIM];
                    		for (int id=0;id<NDIM;id++)
                    		{
                        		dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                        		dr2    += (dr[id]*dr[id]);
                        		double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][tm0];
                        		double cst2 = (MCCoords[id][t1] - MCCoords[id][t0])*MCCosine[id][tm1];
                        		s1 += cst1;
                        		s2 += cst2;
                    		}
                    		double r = sqrt(dr2);
                    		double th1 = acos(s1/r);
                    		double th2 = acos(s2/r);

                    		double b1[NDIM];
                    		double b2[NDIM];
                    		double b3[NDIM];
                    		for (int id=0;id<NDIM;id++)
                    		{
                        		b1[id] = MCCosine[id][tm0];
                        		b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                        		b3[id] = MCCosine[id][tm1];
                    		}
                    		VectorNormalisation(b1);
                    		VectorNormalisation(b2);
                    		VectorNormalisation(b3);

                    		//Calculation of dihedral angle 
                    		double n1[NDIM];
                    		double n2[NDIM];
                    		double mm[NDIM];

                    		CrossProduct(b2, b1, n1);
                    		CrossProduct(b2, b3, n2);
                    		CrossProduct(b2, n2, mm);

                    		double xx = DotProduct(n1, n2);
                    		double yy = DotProduct(n1, mm);

                    		double phi = atan2(yy, xx);
                    		if (phi<0.0) phi += 2.0*M_PI;

                    		//Dihedral angle calculation is completed here
                			double r1 = 0.74;// bond length in Angstrom
							r1 /= BOHRRADIUS;
                    		double r2 = r1;// bond length in bohr
                    		double rd = r/BOHRRADIUS;
                    		double potl;
                    		vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                    		spot_beads += potl*CMRECIP2KL;
                		} //stype
                		if (stype == HF)
                		{
							double Eulang0[NDIM], Eulang1[NDIM];
   							Eulang0[PHI] = MCAngles[PHI][t0];
   							Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   							Eulang0[CHI] = 0.0;
   							Eulang1[PHI] = MCAngles[PHI][t1];
   							Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   							Eulang1[CHI] = 0.0;
                			spot_pair += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
                		} //stype
					}//loop over beads
					spot_pair += spot_beads;
        		}// loop over atoms1 (molecules)
        	}// loop over atoms0 (molecules)
			spot_sector += 0.5*spot_pair;

    		spot_cage = 0.0;
#ifdef CAGEPOT
    		for (int atom0 = atomStart; atom0 < atomEnd; atom0++)
    		{
        		int offset0 = NumbTimes*atom0;

   				spot_beads=0.0;
       			#pragma omp parallel for reduction(+: spot_beads)
       			for (int it = 0; it < NumbTimes; it += (NumbTimes - 1))
				{
        			int t0 = offset0 + it;
        			double cost = MCAngles[CTH][t0];
        			double phi = MCAngles[PHI][t0];
        			if (phi < 0.0) phi = 2.0*M_PI + phi;
        			phi = fmod(phi,2.0*M_PI);
        			int type0   =  MCType[atom0];
        			spot_beads += LPot2DRotDOF(cost,phi,type0);
				}
				spot_cage += spot_beads;
    		}
			spot_sector_cage += 0.5*spot_cage;
#endif
		}
	}
	double spotReturn = 0.5*(spot_sector + spot_sector_cage);
	return spotReturn;
}

double GetRotEnergyPIGS(void)
{
#ifdef PIGSENTBOTH
	int atomStart = NumbAtoms/2;
#else
	int atomStart = 0;
#endif

    double srot = 0.0;
    int type  = IMTYPE;
    double nslice = (((double)NumbRotTimes) - 1.0);

    for (int atom0 = atomStart; atom0 < NumbAtoms; atom0++)
    {
        int offset0 = NumbTimes*atom0;

        int it = ((NumbRotTimes - 1)/2);
        int t0 = offset0 +  it;
        int tm = (t0 - 1);
        int tp = (t0 + 1);

        double p0 = 0.0;
        double p1 = 0.0;
        for (int id=0;id<NDIM;id++)
        {
            p0 += (MCCosine[id][t0]*MCCosine[id][tm]);
            p1 += (MCCosine[id][t0]*MCCosine[id][tp]);
	    } 


#ifdef TYPE0
        double rdens0 = SRotDens(p0,type);
        double rdens1 = SRotDens(p1,type);
        if (fabs(rdens0) > RZERO)               // need to find asymptotic for small rot dens
        {
            srot += SRotDensDeriv(p0,type)/rdens0;
        }
        if (fabs(rdens1) > RZERO)               // need to find asymptotic for small rot dens
        {
            srot += SRotDensDeriv(p1,type)/rdens1;
        }
#endif
#ifdef TYPE1
        srot += SRotDensDeriv(p0,type) + SRotDensDeriv(p1,type);
#endif
	}
    return (0.5*nslice*srot);
}

void GetCosThetaPIMC(double &cosTheta, double *compxyz)
{
    const char *_proc_=__func__;

    // if user passed in a null pointer for array, bail out early!
    if (!compxyz)
        return;

    double scosTheta;
    double scompxyz[NDIM];

    if(MCAtom[IMTYPE].numb > 1)
    {
        scosTheta      = 0.0;
        for (int atom0 = 0; atom0 < (NumbAtoms - 1); atom0++)
        {
            for (int atom1 = (atom0 + 1); atom1 < NumbAtoms; atom1++)
            {
	            double sum_beads = 0.0;
    	    	for (int it = 0; it < NumbRotTimes; it++) // Rotational Time slices, P
        	    {
                	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
                    int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

                    int t0      = offset0 + it;
                    int t1      = offset1 + it;

                    double cst  = 0.0;
                    for (int id = 0; id < NDIM; id++) // X, Y, Z
                    {
                        cst    += MCCosine[id][t0]*MCCosine[id][t1];
                    }
                    sum_beads += cst;
                }
                scosTheta   += sum_beads;
            }     // LOOP OVER ATOM PAIRS
        }
        scompxyz[0] = 0.0;
        scompxyz[1] = 0.0;
        scompxyz[2] = 0.0;

        for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
        {
       		double sum_beadsx = 0.0;
            double sum_beadsy = 0.0;
            double sum_beadsz = 0.0;
            for (int it = 0; it < NumbRotTimes; it++) // Rotational Time slices, P
            {
                int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
                int t0      = offset0 + it;

                sum_beadsx += MCCosine[0][t0];
                sum_beadsy += MCCosine[1][t0];
                sum_beadsz += MCCosine[2][t0];
            }
            scompxyz[0] += sum_beadsx;
            scompxyz[1] += sum_beadsy;
            scompxyz[2] += sum_beadsz;
        }
    }
    if(MCAtom[IMTYPE].numb == 1)
    {
        // Initial configurations //
        double phi1  = 0.0;
        double cost1 = 1.0;
        double sint1 = sqrt(1.0 - cost1*cost1);

        double uvec1[NDIM];
        uvec1[0]     = sint1*cos(phi1);
        uvec1[1]     = sint1*sin(phi1);
        uvec1[2]     = cost1;

        int atom0    = 0;
        int type0    = MCType[atom0];
        int offset0  = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
        double sum_beads = 0.0;
        for (int it = 0; it < NumbRotTimes; it++)
        {
            int tm0      = offset0 + it/RotRatio;

            double cst   = 0.0;
            for (int id = 0; id < NDIM; id++)
            {
                cst    += MCCosine[id][tm0]*uvec1[id];
            }
            sum_beads += cst;
            scompxyz[0] = MCCosine[0][tm0];
            scompxyz[1] = MCCosine[1][tm0];
            scompxyz[2] = MCCosine[2][tm0];
        }
        scosTheta   = sum_beads;
    }

    cosTheta = scosTheta/NumbRotTimes;
    compxyz[0] = scompxyz[0]/(NumbAtoms*NumbRotTimes);
    compxyz[1] = scompxyz[1]/(NumbAtoms*NumbRotTimes);
    compxyz[2] = scompxyz[2]/(NumbAtoms*NumbRotTimes);
}

void GetCosThetaPIGS(double &cosTheta, double *compxyz)
{
    const char *_proc_=__func__; 

    // if user passed in a null pointer for array, bail out early!
    if (!compxyz)
        return;
    int it = (NumbRotTimes - 1)/2;

	double scosTheta;
	double scompxyz[NDIM];

	if(MCAtom[IMTYPE].numb > 1)
	{
        scosTheta    = 0.0;
    	for (int atom0 = 0; atom0 < (NumbAtoms-1); atom0++)
        {    
    	    for (int atom1 = (atom0+1); atom1 < NumbAtoms; atom1++)
    	    {
        	    int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
        	    int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

       		    int t0      = offset0 + it;
        	    int t1      = offset1 + it;

                double cst  = 0.0;
           	    for (int id = 0; id < NDIM; id++)
           	    {    
               	    cst    += MCCosine[id][t0]*MCCosine[id][t1];
           	    }
           	    scosTheta   += cst;
    		}     // LOOP OVER ATOM PAIRS
		}

		for (int id = 0; id<NDIM; id++)
		{
			double sum = 0.0;
    		for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
        	{    
           		int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
       			int t0      = offset0 + it;

				sum += MCCosine[id][t0];
			}
			scompxyz[id] = sum;
		}
	}
	if(MCAtom[IMTYPE].numb == 1)
	{
		// Initial configurations //
        double phi1  = 0.0;
        double cost1 = 1.0;
        double sint1 = sqrt(1.0 - cost1*cost1);

        double uvec1[NDIM];
        uvec1[0]     = sint1*cos(phi1);
        uvec1[1]     = sint1*sin(phi1);
        uvec1[2]     = cost1;

		int atom0    = 0;
     	int type0    = MCType[atom0];
       	int offset0  = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
        int tm0      = offset0 + it/RotRatio;

        double cst   = 0.0;
        for (int id = 0; id < NDIM; id++)
        {    
       	    cst    += MCCosine[id][tm0]*uvec1[id];
        }
		scosTheta   = cst;
		scompxyz[0] = MCCosine[0][tm0];
		scompxyz[1] = MCCosine[1][tm0];
		scompxyz[2] = MCCosine[2][tm0];
	}

	cosTheta = scosTheta;
	for (int id = 0; id < NDIM; id++)
	{
		compxyz[id] = scompxyz[id]/NumbAtoms;
	}
}

void GetCosThetaPIGSENT(double &cosTheta, double *compxyz)
{
    const char *_proc_=__func__; 

    // if user passed in a null pointer for array, bail out early!
    if (!compxyz)
        return;
#ifdef PIGSENTBOTH
	int atomStart = NumbAtoms/2;
#else
	int atomStart = 0;
#endif

    int it = (NumbRotTimes - 1)/2;

	double scosTheta_pair;
	double scompxyz_pair[NDIM];
	double scosTheta_sector;
	double scompxyz_sector[NDIM];

	if(MCAtom[IMTYPE].numb > 1)
	{
		int atomStart, atomEnd;
		for (int isector = 0; isector < 2; isector++)
		{
			if (isector == 0)
			{
				atomStart = 0;
				atomEnd   = NumbAtoms/2;
			}
			if (isector == 1)
			{
				atomStart = NumbAtoms/2;
				atomEnd   = NumbAtoms;
			}

        	scosTheta_pair= 0.0;
    		for (int atom0 = atomStart; atom0 < (atomEnd-1); atom0++)
        	{    
    	    	for (int atom1 = (atom0+1); atom1 < atomEnd; atom1++)
    	    	{
        	    	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
        	    	int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

       		    	int t0      = offset0 + it;
        	    	int t1      = offset1 + it;

                	double cst  = 0.0;
           	    	for (int id = 0; id < NDIM; id++)
           	    	{    
               	    	cst    += MCCosine[id][t0]*MCCosine[id][t1];
           	    	}
           	    	scosTheta_pair += cst;
    			}     // LOOP OVER ATOM PAIRS
			}

			scompxyz_pair[0] = 0.0;
			scompxyz_pair[1] = 0.0;
			scompxyz_pair[2] = 0.0;

    		for (int atom0 = atomStart; atom0 < atomEnd; atom0++)
        	{    
            	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
       			int t0      = offset0 + it;

				scompxyz_pair[0] += MCCosine[0][t0];
				scompxyz_pair[1] += MCCosine[1][t0];
				scompxyz_pair[2] += MCCosine[2][t0];
			}
		}
		scosTheta_sector += 0.5*scosTheta_pair;
		scompxyz_sector[0] += 0.5*scompxyz_pair[0];
        scompxyz_sector[1] += 0.5*scompxyz_pair[1];
        scompxyz_sector[2] += 0.5*scompxyz_pair[2];
	}
	cosTheta = scosTheta_sector;
	compxyz[0] = scompxyz_sector[0]/NumbAtoms;
	compxyz[1] = scompxyz_sector[1]/NumbAtoms;
	compxyz[2] = scompxyz_sector[2]/NumbAtoms;
}

#ifdef DDCORR
void GetDipoleCorrelationPIMC(double *DipoleCorrXYZ, double *DipoleCorrX, double *DipoleCorrY, double *DipoleCorrZ, double *DipoleCorrXY)
{
    const char *_proc_=__func__; 

	if(MCAtom[IMTYPE].numb > 1)
	{
        double totalCorr, xCorr, yCorr, zCorr, xyCorr;

		int ii = 0;
    	for (int atom0 = 0; atom0 < (NumbAtoms - 1); atom0++)
        {    
    	    for (int atom1 = (atom0 + 1); atom1 < NumbAtoms; atom1++)
    	    {
            	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
            	int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

            	totalCorr   = 0.0;
				xyCorr      = 0.0;
				xCorr       = 0.0;
				yCorr       = 0.0;
				zCorr       = 0.0;

#pragma omp parallel for reduction(+: totalCorr, xyCorr, xCorr, yCorr, zCorr)
				for (int it = 0; it < NumbRotTimes; it++)
				{
       	    		int t0      = offset0 + it;
            		int t1      = offset1 + it;

					double totalCorrid = 0.0;
					double xyCorrid    = 0.0;
               		for (int id = 0; id < NDIM; id++)
               		{    
           	    		totalCorrid += MCCosine[id][t0]*MCCosine[id][t1];
						if (id <= 1)
						{
           	    			xyCorrid += MCCosine[id][t0]*MCCosine[id][t1];
						}
               		}

           	    	totalCorr += totalCorrid;
           	    	xyCorr    += xyCorrid;
               		xCorr     += MCCosine[0][t0]*MCCosine[0][t1];
               		yCorr     += MCCosine[1][t0]*MCCosine[1][t1];
               		zCorr     += MCCosine[2][t0]*MCCosine[2][t1];

               		DipoleCorrXYZ[ii] = totalCorr/NumbRotTimes;
               		DipoleCorrX[ii]   = xCorr/NumbRotTimes;
               		DipoleCorrY[ii]   = yCorr/NumbRotTimes;
               		DipoleCorrZ[ii]   = zCorr/NumbRotTimes;
               		DipoleCorrXY[ii]  = xyCorr/NumbRotTimes;
				}
				ii++;
    		}
		}
	}
}
#endif

#ifdef DDCORR
void GetDipoleCorrelationPIGS(double *DipoleCorrXYZ, double *DipoleCorrX, double *DipoleCorrY, double *DipoleCorrZ, double *DipoleCorrXY)
{
    const char *_proc_=__func__; 

    int it = (NumbRotTimes - 1)/2;
	if(MCAtom[IMTYPE].numb > 1)
	{
        double totalCorr, xCorr, yCorr, zCorr, xyCorr;

		int ii = 0;
    	for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
        {    
    	    for (int atom1 = atom0; atom1 < NumbAtoms; atom1++)
    	    {
            	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
            	int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

       	    	int t0      = offset0 + it;
            	int t1      = offset1 + it;

               	totalCorr   = 0.0;
               	for (int id = 0; id < NDIM; id++)
               	{    
           	    	totalCorr += MCCosine[id][t0]*MCCosine[id][t1];
               	}

               	xCorr   = MCCosine[0][t0]*MCCosine[0][t1];
               	yCorr   = MCCosine[1][t0]*MCCosine[1][t1];
               	zCorr   = MCCosine[2][t0]*MCCosine[2][t1];

				xyCorr      = 0.0;
               	for (int id = 0; id < (NDIM-1); id++)
				{
           	    	xyCorr += MCCosine[id][t0]*MCCosine[id][t1];
				}
               	DipoleCorrXYZ[ii] = totalCorr;
               	DipoleCorrX[ii]   = xCorr;
               	DipoleCorrY[ii]   = yCorr;
               	DipoleCorrZ[ii]   = zCorr;
               	DipoleCorrXY[ii]  = xyCorr;
				ii++;
    		}
		}
	}
}
#endif

#ifdef DDCORR
void GetDipoleCorrelationPIGSENT(double *DipoleCorrXYZ, double *DipoleCorrX, double *DipoleCorrY, double *DipoleCorrZ, double *DipoleCorrXY)
{
    const char *_proc_=__func__; 
    int it = (NumbRotTimes - 1)/2;
	if(MCAtom[IMTYPE].numb > 1)
	{
        double totalCorr, xCorr, yCorr, zCorr, xyCorr;

		int atomStart, atomEnd;
		for (int isector = 0; isector < 2; isector++)
		{
			if (isector == 0)
			{
				atomStart = 0;
				atomEnd   = NumbAtoms/2;
			}
			if (isector == 1)
			{
				atomStart = NumbAtoms/2;
				atomEnd   = NumbAtoms;
			}

			int ii = 0;
    		for (int atom0 = atomStart; atom0 < atomEnd; atom0++)
        	{    
    	    	for (int atom1 = atom0; atom1 < atomEnd; atom1++)
    	    	{
            		int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
            		int offset1 = MCAtom[IMTYPE].offset + NumbRotTimes*atom1;

       	    		int t0      = offset0 + it;
            		int t1      = offset1 + it;

               		totalCorr   = 0.0;
               		for (int id = 0; id < NDIM; id++)
               		{    
           	    		totalCorr += MCCosine[id][t0]*MCCosine[id][t1];
               		}

               		xCorr   = MCCosine[0][t0]*MCCosine[0][t1];
               		yCorr   = MCCosine[1][t0]*MCCosine[1][t1];
               		zCorr   = MCCosine[2][t0]*MCCosine[2][t1];

					xyCorr      = 0.0;
               		for (int id = 0; id < (NDIM-1); id++)
					{
           	    		xyCorr += MCCosine[id][t0]*MCCosine[id][t1];
					}	
               		DipoleCorrXYZ[ii] += 0.5*totalCorr;
               		DipoleCorrX[ii]   += 0.5*xCorr;
               		DipoleCorrY[ii]   += 0.5*yCorr;
               		DipoleCorrZ[ii]   += 0.5*zCorr;
               		DipoleCorrXY[ii]  += 0.5*xyCorr;
					ii++;
    			}
			}
		}
	}
}
#endif

double *GetPhiEntanglement()
{
    const char *_proc_=__func__; 

    double *phiInstant = new double [2*NumbAtoms];
	int BeadMminus1 = (((NumbRotTimes - 1)/2) - 1); 
    int BeadM       = ((NumbRotTimes - 1)/2);

	for (int it = BeadMminus1; it <= BeadM; it++) 
	{
		for (int atom = 0; atom < NumbAtoms; atom++)
   		{
			int kk = atom + (it - BeadMminus1)*NumbAtoms;
    		int offset      = MCAtom[IMTYPE].offset + (NumbRotTimes*atom);
       		int tt          = offset + it;
       		phiInstant[kk] = MCAngles[PHI][tt];
		}
    }
    return phiInstant;
}

double *GetCosThetaEntanglement()
{
    const char *_proc_=__func__; 

    double *cosTheta = new double [2*NumbAtoms*NDIM];
	int BeadMminus1 = (((NumbRotTimes - 1)/2) - 1); 
    int BeadM       = ((NumbRotTimes - 1)/2);

	for (int id = 0; id < NDIM; id++)
   	{
		for (int it = BeadMminus1; it <= BeadM; it++) 
		{
            int jj = (it - BeadMminus1) + 2*id;
			for (int atom = 0; atom < NumbAtoms; atom++)
   			{
				int kk = atom + jj*NumbAtoms;
	    		int offset      = MCAtom[IMTYPE].offset + (NumbRotTimes*atom);
        		int tt          = offset + it;
        		cosTheta[kk] = MCCosine[id][tt];
			}
		}
    }
    return cosTheta;
}

double *GetProdUvec12()
{
    const char *_proc_=__func__; 
	int type = IMTYPE;

    double *ProdUvec12 = new double [2*NumbRotTimes*NumbAtoms];

	for (int atom = 0; atom < NumbAtoms; atom++)
   	{
		for (int it0 = 0; it0 <NumbRotTimes; it0++) 
		{
			int it1;
			if (it0 == (NumbRotTimes - 1)) it1 = 0;
			else it1 = it0+1;

	   		int offset      = MCAtom[IMTYPE].offset + (NumbRotTimes*atom);
       		int t0          = offset + it0;
       		int t1          = offset + it1;
			
			double p0 = 0.0;
			double p1 = 0.0;
			for (int id = 0; id < NDIM; id++)
   			{
				p0 += MCCosine[id][t0]*MCCosine[id][t1];
			}
			int kk = it0 + atom*NumbRotTimes;
       		ProdUvec12[2*kk] = p0;
       		ProdUvec12[2*kk+1] = SRotDensDeriv(p0,type);
		}
    }
    return ProdUvec12;
}

double GetPhi(void)
{
    const char *_proc_=__func__;
    if (NumbAtoms <= 1) nrerror(_proc_," Only one rotor/atom/molecule");

    double phi;
    for (int atom0=0;atom0<(NumbAtoms-1);atom0++)
    for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
    {
        int type0   = MCType[atom0];
        int type1   = MCType[atom1];

        int offset0 = NumbTimes*atom0;
        int offset1 = NumbTimes*atom1;

        int it      = (NumbTimes - 1)/2;
        int t0      = offset0 + it;
        int t1      = offset1 + it;

		MCCoords[0][t0] = 0.0;
		MCCoords[1][t0] = 0.0;
		MCCoords[2][t0] = 0.0;
		MCCoords[0][t1] = 0.0;
		MCCoords[1][t1] = 0.0;
		MCCoords[2][t1] = Distance;
        double dr2  = 0.0;
        for (int id = 0; id < NDIM; id++)
        {
            double dr = (MCCoords[id][t0] - MCCoords[id][t1]);
            dr2    += dr*dr;
        }
        double r    = sqrt(dr2);

        if ( ((MCAtom[type0].molecule == 4) && (MCAtom[type1].molecule == 4)) && (MCAtom[IMTYPE].numb > 1) )
        {
            int tm0=offset0 + it/RotRatio;
            int tm1=offset1 + it/RotRatio;

            double b1[NDIM];
            double b2[NDIM];
            double b3[NDIM];
            for (int id=0;id<NDIM;id++)
            {
                b1[id] = MCCosine[id][tm0];
                b2[id] = (MCCoords[id][t1] - MCCoords[id][t0])/r;
                b3[id] = MCCosine[id][tm1];
            }
            VectorNormalisation(b1);
            VectorNormalisation(b2);
            VectorNormalisation(b3);

            //Calculation of dihedral angle 
            double n1[3];
            double n2[3];
            double mm[3];

            CrossProduct(b2, b1, n1);
            CrossProduct(b2, b3, n2);
            CrossProduct(b2, n2, mm);

            double xx = DotProduct(n1, n2);
            double yy = DotProduct(n1, mm);

            phi       = atan2(yy, xx);
            if (phi<0.0) phi += 2.0*M_PI;

        }  // LOOP OVER TIME SLICES
    }     // LOOP OVER ATOM PAIRS
    return phi;
}


double GetPotEnergyEntanglement(int atom0, int atom1)
{
	const char *_proc_=__func__;

    int it      = (NumbRotTimes - 1)/2;

    int offset0 = NumbRotTimes*atom0;
    int offset1 = NumbRotTimes*atom1;
    int t0      = offset0 + it;
    int t1      = offset1 + it;

    double spot;
	double Eulang0[NDIM], Eulang1[NDIM];
   	Eulang0[PHI] = MCAngles[PHI][t0];
   	Eulang0[CTH] = acos(MCAngles[CTH][t0]);
   	Eulang0[CHI] = 0.0;
   	Eulang1[PHI] = MCAngles[PHI][t1];
   	Eulang1[CTH] = acos(MCAngles[CTH][t1]);
   	Eulang1[CHI] = 0.0;
    spot = 0.5*PotFunc(atom0, atom1, Eulang0, Eulang1, it);
    return spot;
}

double GetEstimNM(void)
{
    int atom0, atom1;
    int type        = IMTYPE;

   	int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	int particleA1Max = particleA1Min + NumbParticle - 1;
   	int particleA2Min = particleA1Max + 1;
   	int particleA2Max = particleA2Min + NumbParticle - 1;

    double spot    = 0.0;

    for (int atom0 = particleA1Min; atom0 <= particleA1Max; atom0++)
    {
        for (int atom1 = (particleA2Max+1); atom1 < NumbAtoms; atom1++)
        {
            spot      += GetPotEnergyEntanglement(atom0, atom1);
		}
    }

    for (int atom0 = particleA2Min; atom0 <= particleA2Max; atom0++)
    {
    	for (int atom1 = 0; atom1 < particleA1Min; atom1++)
    	{
        	spot      += GetPotEnergyEntanglement(atom0, atom1);
    	}
	}
    double potEstimNM = exp(-MCRotTau*spot);

    int it0  = (((NumbRotTimes - 1)/2)-1);
    int it1  = ((NumbRotTimes - 1)/2);

    double dens1 = 1.0;
    for (int atom0 = particleA1Min; atom0 <= particleA1Max; atom0++)
	{
    	int atom1 = particleA2Max - (atom0 - particleA1Min);
    	int offset0 = NumbRotTimes*atom0;
    	int offset1 = NumbRotTimes*atom1;

    	int t1M1 = offset0 + it0;
    	int t1M = offset1 + it1;

    	double p0   = 0.0;
    	for (int id = 0;id<NDIM;id++)
   	 	{
        	p0 += (MCCosine[id][t1M1]*MCCosine[id][t1M]);
    	}
    	dens1 *= SRotDens(p0,type);
	}

    double dens2 = 1.0;
    for (int atom0 = particleA2Min; atom0 <= particleA2Max; atom0++)
	{
    	int atom1 = particleA1Max - (atom0 - particleA2Min);

    	int offset0 = NumbRotTimes*atom0;
    	int offset1 = NumbRotTimes*atom1;

    	int t1M1 = offset0 + it0;
    	int t1M = offset1 + it1;

    	double p0   = 0.0;
    	for (int id = 0;id<NDIM;id++)
   	 	{
        	p0 += (MCCosine[id][t1M1]*MCCosine[id][t1M]);
    	}
    	dens2 *= SRotDens(p0,type);
	}
  	double estimNM = dens1*dens2*potEstimNM;
    return estimNM;
}

double GetEstimDM(void)
{
    int type       = IMTYPE;

   	int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	int particleA1Max = particleA1Min + NumbParticle - 1;
   	int particleA2Min = particleA1Max + 1;
   	int particleA2Max = particleA2Min + NumbParticle - 1;

    double spot    = 0.0;

    for (int atom0 = particleA1Min; atom0 <= particleA1Max; atom0++)
	{
    	for (int atom1 = 0; atom1 < particleA1Min; atom1++)
    	{
        	spot      += GetPotEnergyEntanglement(atom0, atom1);
    	}
	}

    for (int atom0 = particleA2Min; atom0 <= particleA2Max; atom0++)
	{
    	for (int atom1 = (particleA2Max+1); atom1 < NumbAtoms; atom1++)
    	{
        	spot      += GetPotEnergyEntanglement(atom0, atom1);
    	}
	}
    double potEstimDM = exp(-MCRotTau*spot);

    double dens = 1.0;
    for (int atom0 = particleA1Min; atom0 <= particleA2Max; atom0++)
    {
        int it0 = (((NumbRotTimes - 1)/2)-1);
        int it1 = ((NumbRotTimes - 1)/2);

        int offset0 = NumbRotTimes*atom0;

        int t0 = offset0 + it0;
        int t1 = offset0 + it1;

        double p0   = 0.0;
        for (int id = 0;id<NDIM;id++)
        {
            p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
        }
        dens *= SRotDens(p0,type);
    }
    double estimDM = dens*potEstimDM;
    return estimDM;
}

double GetPotEnergy(void)
// should be compatible with PotEnergy() from mc_piqmc.cc
{
   const char *_proc_=__func__; //  GetPotEnergy_Densities()  

#ifdef DEBUG_WORM
   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");
#endif

// double dr[NDIM];
   double spot = 0.0;

   for (int atom0=0;atom0<(NumbAtoms-1);atom0++)      
   for (int atom1=(atom0+1);atom1<NumbAtoms;atom1++)
   {
      int type0   = MCType[atom0];
      int type1   = MCType[atom1];

      int offset0 = NumbTimes*atom0;
      int offset1 = NumbTimes*atom1;

      double spot_pair=0.0;

      #pragma omp parallel for reduction(+: spot_pair)
      for (int it=0;it<NumbTimes;it++) 	    
      {  
         double dr[NDIM];
         int t0 = offset0 + it;
         int t1 = offset1 + it;

         double dr2 = 0.0;  		 
         for (int id=0;id<NDIM;id++)
         {
            dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

            if (MINIMAGE)
            dr[id] -= (BoxSize*rint(dr[id]/BoxSize));

            dr2    += (dr[id]*dr[id]);
         }
   	 
//#ifdef _CUTOFF_	     
//       if (dr2<dljcutoff2)
//#endif
         double r = sqrt(dr2);

//----------- [ATOM - MOLECULE] ----------------------

         if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
         {
         //  type 1 is a molecule 
     
             int sgn   = 1;             // set to -1 to correct the orientaion of dr

             int tm    = offset1 + it/RotRatio;
//           int tm    = offset1 + floor((double)it/(double)RotRatio);

             int typep = type1;         // define the type of the potential
             int typed = type0;         // define the type of the density

         //  type 0 is a molecule ?   
             if (MCAtom[type0].molecule == 1)  // does not work for two molecules
             {
                sgn   = -1;   
                tm    = offset0 + it/RotRatio;
                typep = type0; 
                typed = type1; 
             }

             double cost = 0.0;
             for (int id=0;id<NDIM;id++)    // n*dr = r*cos(theta) 
             cost += (MCCosine[id][tm]*dr[id]);   	 
	 
             cost /= r;                     // cos(theta)
             cost *= sgn;                   // correct the orientation 

//           bin_2Ddensity (r,cost,typed);  // densities 
             spot_pair += LPot2D(r,cost,typep);  // potential energy 
        }
//----------- [ATOM - NON-LINEAR MOLECULE] ----------------------
        else if (((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && (MCAtom[type0].molecule != MCAtom[type1].molecule) ) // 3D interaction, no density is calculated now
        {

            int tm;
            int typed;

            double RCOM[3];
            double Rpt[3];
            double Eulang[3];
            double vpot3d;
            double radret;
            double theret;
            double chiret;
            double hatx[3];
            double haty[3];
            double hatz[3];
            int    ivcord = 0;
            if(MCAtom[type0].molecule == 2)
            {
//             determine type of atoms for bin_3Ddensity
               typed = type1;
               tm  = offset0 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  RCOM[id] = MCCoords[id][t0];
                  Rpt[id]  = MCCoords[id][t1];
               }
            }
            else
            {
//             determine type of atoms for bin_3Ddensity
               typed = type0;
               tm  = offset1 + it/RotRatio;
               for (int id=0;id<NDIM;id++)
               {
                  Rpt[id]  = MCCoords[id][t0];
                  RCOM[id] = MCCoords[id][t1];
               }
            }
            Eulang[PHI]=MCAngles[PHI][tm];
            Eulang[CTH]=acos(MCAngles[CTH][tm]);
            Eulang[CHI]=MCAngles[CHI][tm];

            if( ISPHER == 0)
            {
               vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);
            }
            else if( ISPHER == 1)
            {
               radret = r;
               vspher_(&radret,&vpot3d);
               theret = 0.0;
               chiret = 0.0;
            }

//          bin_3Ddensity (radret,theret,chiret,typed);  // accumulate density

            spot_pair += vpot3d;
        }
//---------[NON-LINEAR - NON-LINEAR from GG]
         else if ( ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
         {
          //   if ( (MCType[atom0] == IMTYPE) && (MCType[atom1] == IMTYPE) )
         //   {
         //     if ( (MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2) )
         //   {
        //     cout<<"Le if de GG: MCAtom[type0].molecule MCAtom[type1].molecule"<<MCAtom[type0].molecule<<" "<<MCAtom[type1].molecule<<endl;
             double com_1[3];
             double com_2[3];
             double Eulang_1[3];
             double Eulang_2[3];
             double E_2H2O;
             for (int id=0;id<NDIM;id++)
             {
                  com_1[id] = MCCoords[id][t0];
                  com_2[id] = MCCoords[id][t1];
             }
             int tm0=offset0 + it/RotRatio;
             int tm1=offset1 + it/RotRatio;
             Eulang_1[PHI]=MCAngles[PHI][tm0];
             Eulang_1[CTH]=acos(MCAngles[CTH][tm0]);
             Eulang_1[CHI]=MCAngles[CHI][tm0];
             Eulang_2[PHI]=MCAngles[PHI][tm1];
             Eulang_2[CTH]=acos(MCAngles[CTH][tm1]);
             Eulang_2[CHI]=MCAngles[CHI][tm1];
             caleng_(com_1, com_2, &E_2H2O,
                        Eulang_1, Eulang_2);
             spot_pair += E_2H2O;
          //  }  //
          //  }
         }
//------------- [ATOM - ATOM] ------------------------------- 
         else                // only one atom type    
         if ((type0 == type1) && MCAtom[type0].molecule == 0) // no "cross" densities 
         {
//          bin_1Ddensity (r,type1);    // densities 
            spot_pair += SPot1D(r,type1);    // potential energy
         }
      }  // LOOP OVER TIME SLICES
      spot += spot_pair;
   }     // LOOP OVER ATOM PAIRS

// cout<<"in GetPotDensity"<<" _gr1D[0][80]="<<_gr1D[0][80]<<" _gr1D_sum[0][80]="<<_gr1D_sum[0][80]<<endl;
   return (spot/(double)NumbTimes);
}

double GetKinEnergy(void)
{
#ifdef DEBUG_PIMC
	const char *_proc_=__func__; //  GetKinEnergy() 
#ifdef DEBUG_WORM 
   	if (Worm.exists)
   	nrerror(_proc_," Only for Z-configurations");
#endif
#endif

   	int    numb  = 0;       // atom number counter, grand canonical
   	double r2avr = 0.0;     // <r^2> 

   	for (int atom=0;atom<NumbAtoms;atom++)
   	{  
      	numb ++;            // grand canonical only
 
      	int type    = MCType[atom];
      	int offset0 = NumbTimes*atom;
      	int offset1;    
      
      	int gatom   = MCAtom[type].offset/NumbTimes;   

      	double sum = 0.0;
 
      	#pragma omp parallel for reduction(+: sum)
      	for (int it=0;it<NumbTimes;it++) 
      	{
          	int t0  = offset0 + it;
        
          	offset1 = offset0;
          	if ((MCAtom[type].stat ==  BOSE) && ((it+1) == NumbTimes))          
          	offset1 = NumbTimes*(gatom + PIndex[atom-gatom]);

	  		int t1  = offset1 + (it+1) % NumbTimes; // = offset1

          	for (int dim=0;dim<NDIM;dim++)
	  		{
             	double dr = MCCoords[dim][t0] - MCCoords[dim][t1];

             	if (MINIMAGE)
             	dr  -= (BoxSize*rint(dr/BoxSize));

             	sum += (dr*dr);
          	}    
       	} // END loop over time slices

       	r2avr += (sum/(4.0*MCBeta*MCAtom[type].lambda)); 

   	}    // END loop over atoms

#ifdef DEBUG_PIMC
   	if (numb != NumbAtoms)   // should be removed for grand canonical calculations            
   	nrerror(_proc_,"Wrong number of atoms");
#endif

// 	r2avr /= (double)numb;

   	double kin = (double)NumbTimes*Temperature*(0.5*(double)(NDIM*numb) - r2avr);

   	return kin;
}

double GetRotEnergy(void)
{
   int type = IMTYPE; 

   int offset = MCAtom[type].offset;
 
   int atom  = 0;                   // only one molecular impurtiy
   offset   = NumbTimes*atom;
   int gatom = offset/NumbTimes;    // the same offset for rot and trans degrees

   double srot = 0.0;
   ErotSQ=0.0; // a global variable
   Erot_termSQ=0.0; // a global variable

   for (int it0=0;it0<NumbRotTimes;it0++)
   {
      int t0 = offset +  it0;
      int t1 = offset + (it0 + 1) % NumbRotTimes;

      double p0 = 0.0;
      for (int id=0;id<NDIM;id++)
      p0 += (MCCosine[id][t0]*MCCosine[id][t1]);

      if(RotDenType == 0)
      {
         double rdens = SRotDens(p0,type);

         if  (fabs(rdens)>RZERO)               // need to find asymptotic for small rot dens
	   srot += (SRotDensDeriv(p0,type)/rdens);
         Erot_termSQ += (SRotDensDeriv(p0,type)/rdens)*(SRotDensDeriv(p0,type)/rdens);
         ErotSQ += SRotDensEsqrt(p0,type)/rdens;
      }
      else if(RotDenType == 1)
      {
         double rho,erot;
         rsline_(&X_Rot,&p0,&MCRotTau,&rho,&erot);
//       srot += erot/(double)NumbRotTimes;
         srot += rho;
      }
   }

   if(RotDenType == 1)
   {
      srot = srot/(double)NumbRotTimes;
      srot = srot/MCRotTau + 1.0/MCRotTau;
   }
   
   return (srot);	      
}

double GetRotPlanarEnergy(void)
{
    if(RotDenType != 0) {cerr<<"only SOS rho"<<endl; exit(0);}
  
    int type = IMTYPE; 
    int offset = MCAtom[type].offset;

    double ERotPlanar=0.0;
    // changed by PN below
    ErotSQ=0.0;  // a global variable
    Erot_termSQ=0.0;  // a global variable

    for (int atom = 0; atom<MCAtom[type].numb; atom++)                   // multi molecular impurtiy
    {
        offset    = (NumbTimes*atom);

        double srot = 0.0;
        double sesq = 0.0;
        double se_termsq=0.0;
       
#pragma omp parallel for reduction(+: srot,sesq,se_termsq)       
        for (int it0 = 0; it0 < NumbRotTimes; it0++)
	    {
	        int t0 = offset +  it0;
			int it1 = (it0+1);
			if (it1 == NumbRotTimes) it1 = 0;
	        int t1 = offset + it1;
	   
            double p0 = 0.0;
	        for (int id = 0; id < NDIM; id++)
			{
	        	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
			}

#ifdef TYPE0
	        double rdens = SRotDens(p0,type);

	        if (fabs(rdens) > RZERO)               // need to find asymptotic for small rot dens
	        srot += (SRotDensDeriv(p0,type)/rdens);
#endif
#ifdef TYPE1
	        srot += SRotDensDeriv(p0,type);
#endif
#ifdef IOWRITE
	        se_termsq += (SRotDensDeriv(p0,type)/rdens)*(SRotDensDeriv(p0,type)/rdens);
	        sesq += SRotDensEsqrt(p0,type)/rdens;
#endif
	    }
        //      srot = srot / ((double)(NumbRotTimes));
        // sesq = sesq / ((double)(NumbRotTimes)*(double)(NumbRotTimes));
        // se_termsq = se_termsq/ ((double)(NumbRotTimes)*(double)(NumbRotTimes));
        // TOBY question: should I divide by NumbRotTimes above?
      
        ERotPlanar += srot;
#ifdef IOWRITE
        ErotSQ += sesq;
        Erot_termSQ += se_termsq; 
#endif
    }
    return (ERotPlanar);	      
}

double GetRotE3D(void)
{
  int type = IMTYPE;

  int offset = MCAtom[type].offset;

  double ERot3D=0.0;

  ErotSQ=0.0;
  Erot_termSQ=0.0;

  for(int atom  = 0;atom<MCAtom[type].numb;atom++)                   // multi molecular impurtiy
    {
      offset   = NumbTimes*atom;
      int gatom = offset/NumbTimes;    // the same offset for rot and trans degrees
       
      double srot = 0.0;
      double sesq = 0.0;
      double se_termsq=0.0;

      int RNskip;
      if(RotDenType == 0)
	{
	  RNskip = 1;
	}
      else if(RotDenType == 1)
	{
	  RNskip = RNratio;
	}
   
#pragma omp parallel for reduction(+: srot,sesq) // TOBY question: why not include se_termsq in deduction?
      for (int it0=0;it0<NumbRotTimes;it0=it0+RNskip)
	{
	  int t0 = offset +  it0;
	  int t1 = offset + (it0 + RNskip) % NumbRotTimes;
	  //    Given the two sets of Euler angles at t0 and t1, Toby calculates srot and sesq
	  double rho;
	  double erot;
	  double esq;
	  int istop=0;
	  double Eulan1[3];
	  double Eulan2[3];
	  double Eulrel[3];
	  double therel,phirel,chirel;

	  Eulan1[0]=MCAngles[PHI][t0];
	  Eulan1[1]=acos(MCAngles[CTH][t0]);
	  Eulan1[2]=MCAngles[CHI][t0];
	  Eulan2[0]=MCAngles[PHI][t1];
	  Eulan2[1]=acos(MCAngles[CTH][t1]);
	  Eulan2[2]=MCAngles[CHI][t1];
	   
	  rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);
	  phirel=Eulrel[0];
	  therel=Eulrel[1];
	  chirel=Eulrel[2];	   
	  int bin_t    = (int)floor(therel/_delta_theta);
	  if ((bin_t<MC_BINST) && (bin_t>=0))
	    _relthe_sum[bin_t] +=(double)RNskip;	   
	  int bin_p = (int)floor(phirel/_delta_chi);
	  if ((bin_p<MC_BINSC) && (bin_p>=0))
	    _relphi_sum[bin_p] +=(double)RNskip;	   
	  int bin_c = (int)floor(chirel/_delta_chi);
	  if ((bin_c<MC_BINSC) && (bin_c>=0))
	    _relchi_sum[bin_c] +=(double)RNskip;	   
	  //    Rattle Shake rotational energy
	  if(RotDenType == 1 && RNratio == 1)
	    rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);
	  //    srot += erot;
	  if(RotDenType == 1 && RNratio == 1)
	    {
	      srot += rho;
	    }
	  else
	    srot += erot;	   
	  sesq += esq;
	  se_termsq += erot*erot;	   
	}       
      srot = srot / ((double)(NumbRotTimes/RNskip));
      sesq = sesq / ((double)(NumbRotTimes/RNskip)*(double)(NumbRotTimes/RNskip));
      se_termsq = se_termsq/ ((double)(NumbRotTimes/RNskip)*(double)(NumbRotTimes/RNskip));
      
      ERot3D += srot;
      ErotSQ += sesq;
      Erot_termSQ += se_termsq;       
      if(RotDenType == 1 && RNratio == 1)      // Rattle Shake rotational energy
	{
	  ERot3D = ERot3D/(4.0*(MCRotTau/WNO2K)*(MCRotTau/WNO2K)); 
	  ERot3D += 0.25*(X_Rot+Y_Rot+Z_Rot) + 1.5/(MCRotTau/WNO2K);	   
	  ERot3D = ERot3D/WNO2K;
	}  
    }  
  return (ERot3D);
}

/* reactive */
void GetRCF(void) // rotational correlation function //
{
	int type = IMTYPE; 
	for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
	{
   		int offset0 = MCAtom[type].offset + (NumbTimes*atom0);

		for (int it0 = 0; it0 < NumbRotTimes; it0++)
		{
			int t0 = offset0 +  it0;

			for (int itc = 0; itc < NumbRotTimes; itc++)  
			{
				int tc = offset0 + (it0 + itc) % NumbRotTimes;

				double p0  = 0.0;

				for (int id = 0; id < NDIM; id++)
				{
					p0  += (MCCosine[id][t0]*MCCosine[id][tc]);
				}
				_rcf     [0][itc]  += p0;   // block average
				_rcf_sum [0][itc]  += p0;   // total average

				for (int in = 1; in < NUMB_RCF; in++)       // rcf[0][] should be the same as rcf[1][]
				{
					double pleg = 1.0;

					//if (p0<PLONE)
					//pleg = gsl_sf_legendre_Pl(in,p0); // inefficient
					_rcf     [in][itc] += pleg;
					_rcf_sum [in][itc] += pleg;
				}
			} // END offsets
		}    // END average over the time origin

#ifdef ROTCORR
		for (int it0 = 0; it0 < NumbRotTimes; it0++)
		{
			int t0 = offset0 +  it0;

			for (int itc = 0; itc < NumbRotTimes; itc++)  
			{
				int tc = offset0 + (it0 + itc) % NumbRotTimes;

				double p0x = 0.0;
				double p0y = 0.0;

				for (int id = 0; id < NDIM; id++)
				{
					p0x += (MCCosinex[id][t0]*MCCosinex[id][tc]);
					p0y += (MCCosiney[id][t0]*MCCosiney[id][tc]);
				}
				_rcfijx[atom0][atom0][0][itc]     +=p0x;
				_rcfijx_sum[atom0][atom0][0][itc] +=p0x;
				_rcfijy[atom0][atom0][0][itc]     +=p0y;
				_rcfijy_sum[atom0][atom0][0][itc] +=p0y;
				_rcfijz[atom0][atom0][0][itc]     +=p0;
				_rcfijz_sum[atom0][atom0][0][itc] +=p0;
				_rcfx     [0][itc] += p0x;   // block average x
				_rcfx_sum [0][itc] += p0x;   // total average x
				_rcfy     [0][itc] += p0y;   // block average y
				_rcfy_sum [0][itc] += p0y;   // total average y

				for (int in = 1; in < NUMB_RCF; in++)       // rcf[0][] should be the same as rcf[1][]
				{
					double pleg = 1.0;

					//if (p0<PLONE)
					//pleg = gsl_sf_legendre_Pl(in,p0); // inefficient

					_rcfijx[atom0][atom0][in][itc]     +=pleg;
					_rcfijx_sum[atom0][atom0][in][itc] +=pleg;
					_rcfijy[atom0][atom0][in][itc]     +=pleg;
					_rcfijy_sum[atom0][atom0][in][itc] +=pleg;
					_rcfijz[atom0][atom0][in][itc]     +=pleg;
					_rcfijz_sum[atom0][atom0][in][itc] +=pleg;

					_rcfx     [in][itc] += pleg;
					_rcfx_sum [in][itc] += pleg;
					_rcfy     [in][itc] += pleg;
					_rcfy_sum [in][itc] += pleg;
				}
			} // END offsets
		}    // END average over the time origin
#endif

//I J RCF
#ifdef ROTCORR
		for (int atom1 = (atom0+1); atom1 < NumbAtoms; atom1++)
		{
			int offset1 = (NumbTimes*atom1);

			for (int it0 = 0; it0 < NumbRotTimes; it0++)
			{
				int t0 = offset0 + it0;
   				int t1 = offset1 + it0;

   				for (int itc = 0; itc < NumbRotTimes; itc++)  // offsets
   				{
       				int tc = offset0 + (it0 + itc) % NumbRotTimes;
       				int tc1= offset1 + (it0 + itc) % NumbRotTimes;

          			double p0ijz = 0.0;
          			double p0ijx = 0.0;
          			double p0ijy = 0.0;

          			for (int id=0;id<NDIM;id++)
          			{
          				p0ijz += (MCCosine[id][t0]*MCCosine[id][tc1]);
          				p0ijx += (MCCosinex[id][t0]*MCCosinex[id][tc1]);
          				p0ijy += (MCCosiney[id][t0]*MCCosiney[id][tc1]);
          			}
         			_rcfijz     [atom0][atom1][0][itc] += p0ijz;   // block average z
         			_rcfijz_sum [atom0][atom1][0][itc] += p0ijz;   // total average z
         			_rcfijx     [atom0][atom1][0][itc] += p0ijx;   // block average x
         			_rcfijx_sum [atom0][atom1][0][itc] += p0ijx;   // total average x
         			_rcfijy     [atom0][atom1][0][itc] += p0ijy;   // block average y
         			_rcfijy_sum [atom0][atom1][0][itc] += p0ijy;   // total average y
         			_rcfijz     [atom1][atom0][0][itc] += p0ijz;   // block average z
         			_rcfijz_sum [atom1][atom0][0][itc] += p0ijz;   // total average z
         			_rcfijx     [atom1][atom0][0][itc] += p0ijx;   // block average x
         			_rcfijx_sum [atom1][atom0][0][itc] += p0ijx;   // total average x
         			_rcfijy     [atom1][atom0][0][itc] += p0ijy;   // block average y
         			_rcfijy_sum [atom1][atom0][0][itc] += p0ijy;   // total average y

          			for (int in	= 1;in < NUMB_RCF; in++)       // rcf[0][] should be the same as rcf[1][]
          			{
              			double pleg = 1.0;

             			_rcfijz     [atom0][atom1][in][itc] += pleg;
             			_rcfijz_sum [atom0][atom1][in][itc] += pleg;
             			_rcfijx     [atom0][atom1][in][itc] += pleg;
             			_rcfijx_sum [atom0][atom1][in][itc] += pleg;
             			_rcfijy     [atom0][atom1][in][itc] += pleg;
             			_rcfijy_sum [atom0][atom1][in][itc] += pleg;
             			_rcfijz     [atom1][atom0][in][itc] += pleg;
             			_rcfijz_sum [atom1][atom0][in][itc] += pleg;
             			_rcfijx     [atom1][atom0][in][itc] += pleg;
             			_rcfijx_sum [atom1][atom0][in][itc] += pleg;
             			_rcfijy     [atom1][atom0][in][itc] += pleg;
						_rcfijy_sum [atom1][atom0][in][itc] += pleg;

					}
				} // END offsets
			}    // END average over the time origin
		}//END loop aver atom1
#endif
	}//END loop over atom0
#ifdef IOWRITE
   int atom  = 0;                    // only one molecular impurtiy
   offset   = NumbTimes*atom;
   int gatom = offset/NumbTimes;

	for (int it0=0;it0<NumbRotTimes;it0++)    
	{
      	int t0 = offset +  it0;

      	for (int itc=0;itc<NumbRotTimes;itc++)  // offsets
      	{	
          	int tc = offset + (it0 + itc) % NumbRotTimes;

          	double p0 = 0.0;
          	for (int id=0;id<NDIM;id++)
          	p0 += (MCCosine[id][t0]*MCCosine[id][tc]);

         	_rcf     [0][itc] += p0;   // block average
         	_rcf_sum [0][itc] += p0;   // total average  

	  		for (int in=1;in<NUMB_RCF;in++)       // rcf[0][] should be the same as rcf[1][]
          	{
              	double pleg = 1.0;
 
	      		if (p0<PLONE)
//	      		pleg = gsl_sf_legendre_Pl(in,p0); // inefficient 	  

             	_rcf     [in][itc] += pleg;     
             	_rcf_sum [in][itc] += pleg;     
	  		}
      	} // END offsets	
   	}    // END average over the time origin 
#endif
}

void SaveRCF(const char fname [], double acount, int mode)
//  save rotational correlation functions
//  
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
	fstream fid;
	string  frcf;
 
	frcf  = fname;
	if (mode == MC_TOTAL)    // accumulated averages
   	{
    	frcf  += IO_SUM;
   	}
	frcf  += IO_EXT_RCF;

  	double norm = acount*(double)NumbRotTimes;
  	double ** rcf_save;
	rcf_save  = _rcf;

  	if (mode == MC_TOTAL)    // accumulated averages
    {
     	rcf_save  = _rcf_sum;
    }

    fid.open(frcf.c_str(),ios::out); io_setout(fid);
    for (int it = 0; it < NumbRotTimes; it++)    // save <n(tau)n(0)>
    {
        fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;
        fid << setw(IO_WIDTH) << rcf_save[0][it % NumbRotTimes]/norm << BLANK;
        fid << endl;

    }
    fid << endl;  // gnuplot index : at list two blank lines
    fid << endl;
    fid << COMMENTS << endl;

    for (int it = 0; it < NumbRotTimes; it++)              // save <Pl(nn)>
    {
        fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;

        for (int ip = 1; ip < NUMB_RCF; ip++)
		{
        	fid << setw(IO_WIDTH) << rcf_save[ip][it % NumbRotTimes]/norm << BLANK;
		}

        fid << endl;
    }
    fid.close();

#ifdef ROTCORR
	fstream fidx, fidy;
	string  frcfx, frcfy;
 
	frcfx = fname;
  	frcfy = fname;
	frcfx += IO_x;
	frcfy += IO_y;

	if (mode == MC_TOTAL)    // accumulated averages
   	{
    	frcfx += IO_SUM;  
    	frcfy += IO_SUM;
   	}
  	frcfx += IO_EXT_RCF;
  	frcfy += IO_EXT_RCF;  

  	double ** rcfx_save;
  	double ** rcfy_save;

  	rcfx_save = _rcfx;
  	rcfy_save = _rcfy;

  	if (mode == MC_TOTAL)    // accumulated averages
    {
     	rcfx_save = _rcfx_sum;
     	rcfy_save = _rcfy_sum;
    }

    fidx.open(frcfx.c_str(),ios::out); io_setout(fidx);
    fidy.open(frcfy.c_str(),ios::out); io_setout(fidy);

    for (int it = 0; it < NumbRotTimes; it++)    // save <n(tau)n(0)>
    {
        fidx << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;
        fidx << setw(IO_WIDTH) << rcfx_save[0][it % NumbRotTimes]/norm << BLANK;
        fidx << endl;

        fidy << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;
        fidy << setw(IO_WIDTH) << rcfy_save[0][it % NumbRotTimes]/norm << BLANK;
        fidy << endl;
    }
    fidx << endl;  // gnuplot index : at list two blank lines
    fidx << endl;
    fidx << COMMENTS << endl;

    fidy << endl;  // gnuplot index : at list two blank lines
    fidy << endl;
    fidy << COMMENTS << endl;

    for (int it = 0; it < NumbRotTimes; it++)              // save <Pl(nn)>
    {
        fidx << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;
        fidy << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK;

        for (int ip = 1; ip < NUMB_RCF; ip++)
		{
        	fidx << setw(IO_WIDTH) << rcfx_save[ip][it % NumbRotTimes]/norm << BLANK;
        	fidy << setw(IO_WIDTH) << rcfy_save[ip][it % NumbRotTimes]/norm << BLANK;
		}

        fidx << endl;
        fidy << endl;
    }
    fidx.close();
    fidy.close();
#endif

#ifdef IOWRITE
	string frcf0x,frcf0y,frcf0z
	ofstream fid0x,fid0y,fid0z
	const char fnamex [] = "_ITACF_X";
	const char fnamey [] = "_ITACF_Y";
	const char fnamez [] = "_ITACF_Z";
	frcf0x =MCFileName + fnamex;
	frcf0y =MCFileName + fnamey;
	frcf0z =MCFileName + fnamez;

	if (mode == MC_TOTAL)    // accumulated averages
   	{
  		frcf0x += IO_SUM;
  		frcf0y += IO_SUM;
  		frcf0z += IO_SUM;
   	}
  	fid0x.open(frcf0x.c_str(),ios::app);  
  	fid0y.open(frcf0y.c_str(),ios::app);
  	fid0z.open(frcf0z.c_str(),ios::app);
    for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
	{
    	for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
     	{
     		fid0x << atom0 << BLANK << atom1 <<BLANK;
     		fid0y << atom0 << BLANK << atom1 <<BLANK;
     		fid0z << atom0 << BLANK << atom1 <<BLANK;
     		for (int it = 0; it < NumbRotTimes; it++) 
     		{
     			fid0x << _rcfijx[atom0][atom1][0][it % NumbRotTimes]/norm<<BLANK;
     			fid0y << _rcfijy[atom0][atom1][0][it % NumbRotTimes]/norm<<BLANK;
     			fid0z << _rcfijz[atom0][atom1][0][it % NumbRotTimes]/norm<<BLANK;
     		} 
		}
	}
    fid0x << endl;
    fid0y << endl;
    fid0z << endl;
  
  	fid0x.close();
  	fid0y.close();
  	fid0z.close();
#endif
}
/* reactive */

#ifdef IOWRITE
void SaveRCF(const char fname [], double acount, int mode)
//  save rotational correlation functions
//  
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
	fstream fid;
  	string  frcf;

  	frcf  = fname;

  	if (mode == MC_TOTAL)    // accumulated averages
  	frcf += IO_SUM; 

  	frcf += IO_EXT_RCF;

  	fid.open(frcf.c_str(),ios::out); io_setout(fid);

  	double norm = acount*(double)NumbRotTimes;

  	double ** rcf_save;

  	rcf_save = _rcf;
  	if (mode == MC_TOTAL)    // accumulated averages
  	rcf_save = _rcf_sum;

  	for (int it=0;it<=NumbRotTimes;it++)    // save <n(tau)n(0)>
  	{	  
     	fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK; 
     	fid << setw(IO_WIDTH) << rcf_save[0][it % NumbRotTimes]/norm << BLANK;
 
     	fid << endl;
  	}

  	fid << endl;  // gnuplot index : at list two blank lines
  	fid << endl;
  	fid << COMMENTS << endl;

  	for (int it=0;it<=NumbRotTimes;it++)              // save <Pl(nn)>
  	{	  
     	fid << setw(IO_WIDTH) << (double)it*MCRotTau << BLANK; 
 
     	for (int ip=1;ip<NUMB_RCF;ip++) 
     	fid << setw(IO_WIDTH) << rcf_save[ip][it % NumbRotTimes]/norm << BLANK; 
    
     	fid << endl; 
  	}

  	fid.close();
}
#endif
/* reactive */

double GetConfPoten_Densities(void)
// should be compatible with ConfPot() from mc_piqmc.cc
{
   const char *_proc_=__func__; //  GetPotEnergy_Densities()  

   if (Worm.exists)
   nrerror(_proc_," Only for Z-configurations");

   double spot = 0.0;  		
 
   for (int atom=0;atom<NumbAtoms;atom++) 	   
   {
      int offset = NumbTimes*atom;
      int type   = MCType[atom];

      for (int it=0;it<NumbTimes;it++) 	    
      { 
//       bool wline = true;                  // skip if the time slice between ira and masha

//       if (WORM && Worm.exists && (Worm.type == type))  
//       wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
          
//       if (wline)
         {
         double r2 = 0.0;  		
         for (int id=0;id<NDIM;id++)
         r2   += (MCCoords[id][offset+it]*MCCoords[id][offset+it]);

         spot += (MCAtom[type].mass*r2);

         bin_1Ddensity (sqrt(r2),type);    // densities 

         }
      }
   }

   return (0.5*HOMEGA*HOMEGA*spot/(double)NumbTimes);
}

void bin_2Ddensity(double r, double cost, int dtype)
// (r, cost)  ==  (radius, cos(theta))
//  dtype     ==   density type
{
   int bin_r = (int)floor((r-_min_radius)/_delta_radius);

   if ((bin_r<MC_BINSR) && (bin_r>=0))
   { 
      double theta = acos(cost);
      int bin_t    = (int)floor(theta/_delta_theta);

      if ((bin_t<MC_BINST) && (bin_t>=0))
      {
         _gr2D[dtype][bin_r][bin_t] += 1.0;  // block average
         _gr2D_sum[dtype][bin_r][bin_t] += 1.0; // total average //added by Hui Li
      }
   }
}

void bin_3Ddensity(double r, double theta, double chi, int dtype)
// dtype == density type
{

   int bin_r = (int)floor((r-_min_radius)/_delta_radius);
   if ((bin_r<MC_BINSR) && (bin_r>=0))
   {
      int bin_t    = (int)floor(theta/_delta_theta);
      if ((bin_t<MC_BINST) && (bin_t>=0))
      {
         int bin_c = (int)floor(chi/_delta_chi);
         if((bin_c < MC_BINSC) && (bin_c >= 0))
         {
           int ijk = ((bin_r * MC_BINST) + bin_t) * MC_BINSC + bin_c;
           _gr3D[dtype][ijk] += 1.0; // block average
           _gr3D_sum[dtype][ijk] += 1.0; // total average
         }
      }

   }

}

void bin_1Ddensity(double r,int dtype)
//  r         ==   radius
//  dtype     ==   density type
{
   int bin_r = (int)floor((r-_min_radius)/_delta_radius);

   if ((bin_r<MC_BINSR) && (bin_r>=0))
   {
     _gr1D[dtype][bin_r] += 1.0;
     _gr1D_sum[dtype][bin_r] += 1.0;
   }
}

#ifdef HISTOGRAM
void binxyz_1Ddensity(double xyz,int id)
{
   int bin_xyz = (int)floor((xyz-_min_xyz)/_delta_xyz);

   if ((bin_xyz<MC_BINSXYZ) && (bin_xyz>=0))
   {
     _gxyz1D[id][bin_xyz] += 1.0;
     _gxyz1D_sum[id][bin_xyz] += 1.0;
   }
}
#endif

void SaveGraSum(const char fname [], double acount)
// accumulate sum for inter-atomic distribution.  should be similar to the pair distribution in SaveDensities1D
{
  fstream fid;
  string fdens;

  fdens  = fname;
  fdens += IO_SUM;
  fdens += IO_EXT_GRA;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

#ifndef NEWDENSITY
  double norma = _delta_radius*acount*(double)(NumbTimes);
#else
  double norma = _delta_radius*acount;
#endif

  double r;

  for (int ir=0;ir<MC_BINSR;ir++)
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;

     fid<<setw(IO_WIDTH)<<(r*Units.length)<<BLANK;

     for (int id=0;id<NUMB_DENS1D;id++)
     {
//      double nfact = norma*(double)(MCAtom[id].numb*(MCAtom[id].numb-1));
//      the following scaling is to let the gra_sum to be normalized to one by integrating over dr, without any jacobian factor
//      norma = norma * (MCAtom[id].numb*(MCAtom[id].numb-1))/2.0;
#ifndef NEWDENSITY
        fid <<setw(IO_WIDTH)<<_gr1D_sum[id][ir]/(norma*(MCAtom[id].numb*(MCAtom[id].numb-1))/2.0)<<BLANK;   // gra_sum
#else
        fid <<setw(IO_WIDTH)<<_gr1D_sum[id][ir]/norma<<BLANK;   // gra_sum
#endif
     }

     fid<<endl;
  }
//cout<<"norma="<<norma<<" _gr1D_sum[0][80]="<<_gr1D_sum[0][80]<<" numb="<<<<endl;

  fid.close();
}

#ifdef HISTOGRAM
void SaveGxyzSum(const char fname [], double acount)
// accumulate sum for inter-atomic distribution.  should be similar to the pair distribution in SaveDensities1D
{
  fstream fid;
  string fdens;

  fdens  = fname;
  fdens += IO_SUM;
  fdens += IO_EXT_GXYZ;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double norma = _delta_xyz*acount;

  double xyz;

  for (int ixyz=0;ixyz<MC_BINSXYZ;ixyz++)
  {
     xyz   =  (double)ixyz*_delta_xyz;
     xyz  +=  (0.5*_delta_xyz);
     xyz  +=  _min_xyz;

     fid<<setw(IO_WIDTH)<<xyz<<BLANK;

     for (int id=0;id<NDIM;id++)
     {
        fid <<setw(IO_WIDTH)<<_gxyz1D_sum[id][ixyz]/norma<<BLANK;   // gra_sum
     }

     fid<<endl;
  }

  fid.close();
}
#endif

void SaveDensities1D(const char fname [], double acount)
// the density type corresponds to the atom type
{
  fstream fid;
  string fdens;

// ------- 1D radial (pair) distribution functions -----------------------------

  fdens  = fname;
  fdens += IO_EXT_GRA;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double volume = pow(BoxSize,(double)NDIM);     // 3D only  

#ifndef NEWDENSITY
  double norm0  = 2.0*M_PI*_delta_radius*acount  // normalization factor for 
                *(double)(NumbTimes)/volume;     // radial distribution functions 
#else
  double norm0 = _delta_radius*acount;
#endif

  double r,r2;

  for (int ir=0;ir<MC_BINSR;ir++) // normalization
  {	  
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =  r*r;
   
     fid<<setw(IO_WIDTH)<<(r*Units.length)<<BLANK; 
 
     for (int id=0;id<NUMB_DENS1D;id++)
     { 
        double nfact = norm0*(double)(MCAtom[id].numb*(MCAtom[id].numb-1));
#ifndef NEWDENSITY
        fid <<setw(IO_WIDTH)<<_gr1D[id][ir]/(r2*nfact)<<BLANK;   // gra
#else
        fid <<setw(IO_WIDTH)<<_gr1D[id][ir]/norm0<<BLANK;   // gra
#endif
     } 

     fid<<endl;
  }

  fid.close();

// CONVERT 2D TO 1D
  
  if (IMPURITY && MCAtom[IMTYPE].molecule == 1)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;
  fdens += IO_EXT_GRI;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm1 = _delta_radius*(double)NumbTimes     // the normalization for the RADIAL
                  *acount;                            // distribution around impurity

  norm1 *= (4.0*M_PI*pow(Units.length,(double)NDIM)); // 3D only

  for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =   r*r;
   
     fid << setw(IO_WIDTH) << (r*Units.length) << BLANK; 

     for (int id=0;id<NUMB_DENS1D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
     {
        double densr = 0.0;          
        for (int it=0;it<MC_BINST;it++)
        densr += _gr2D[id][ir][it]; 
  
        fid << setw(IO_WIDTH) << densr/(norm1*r2) << BLANK;  
     }
      
     fid << endl;
   } 

   fid.close();

//fid << endl;  // gnuplot index : at list two blank lines
//fid << endl;

// ------  1D angular distributions around an impurity --------------------------

   fdens  = fname;
   fdens += IO_EXT_GRT;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm2 = (double)NumbTimes*acount*_delta_theta;

   double theta;

   for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS1D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
      {
         double denst = 0.0;
         for (int ir=0;ir<MC_BINSR;ir++)  // total contribution
         denst += _gr2D[id][ir][it];

         fid << setw(IO_WIDTH) << denst/(norm2*(double)MCAtom[id].numb) << BLANK; 
      }

      fid << endl;
   }
 
   fid.close();
   } // END if (IMPURITY) convert 2D into 1D
}

void SaveRho1D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
// this subroutine is designed for nonlinear rotor only and should not be called for linear rotor under
// ALL circumstances
{
  fstream fid;
  string fdens;

  double ** dens_save;

  dens_save = _gr3D;
  if(mode == MC_TOTAL)
  dens_save = _gr3D_sum;


// convert 3D density to 1D and save
  if (IMPURITY && MCAtom[IMTYPE].molecule == 2)  //
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if(mode == MC_TOTAL) // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_GRI;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm1 = _delta_radius*(double)NumbTimes     // the normalization for the RADIAL
                  *acount;                            // distribution around impurity

  double r,r2;

  for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
  {
     r   =  (double)ir*_delta_radius;
     r  +=  (0.5*_delta_radius);
     r  +=  _min_radius;
     r2  =   r*r;

     fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;

     for (int id=0;id<NUMB_DENS3D;id++) // rho(r,theta,chi) -> rho(r), rho(theta) and rho(chi) convert 
     {
        double densr = 0.0;
        for (int it=0;it<MC_BINST;it++)
        for (int ic=0;ic<MC_BINSC;ic++)
        {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           densr += dens_save[id][ijk];
        }

//      fid << setw(IO_WIDTH) << densr/(norm1*r2) << BLANK;  
        fid << setw(IO_WIDTH) << densr/(norm1) << BLANK;

     }

     fid << endl;
   }

   fid.close();

//fid << endl;  // gnuplot index : at list two blank lines
//fid << endl;

// ------  1D theta distributions around an impurity --------------------------

   fdens  = fname;

   if(mode == MC_TOTAL) // accumulated averages
   fdens += IO_SUM;

   fdens += IO_EXT_GRT;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm2 = (double)NumbTimes*acount*_delta_theta*(180.0/M_PI); // the last factor converts density to per degree

   double theta;

   for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS3D;id++) // rho(r,cost) -> rho(r) and rho(t) convert 
      {
         double denst = 0.0;
         for (int ir=0;ir<MC_BINSR;ir++)  // total contribution
         for (int ic=0;ic<MC_BINSC;ic++)
         {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           denst += dens_save[id][ijk];
         }

         fid << setw(IO_WIDTH) << denst/(norm2*(double)MCAtom[id].numb) << BLANK;
      }

      fid << endl;
   }

   fid.close();

// ------  1D chi distributions around an impurity --------------------------
   fdens  = fname;

   if(mode == MC_TOTAL) // accumulated averages
   fdens += IO_SUM;

   fdens += IO_EXT_GRC;

   fid.open(fdens.c_str(),ios::out); io_setout(fid);

   double norm4 = (double)NumbTimes*acount*_delta_chi*(180./M_PI); // the last factor is to convert the density in the chi element of degree.

   double chi;

   for (int ic=0;ic<MC_BINSC;ic++)
   {
      chi  = (double)ic*_delta_chi;
      chi += (0.5*_delta_chi);

      fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;

      for (int id=0;id<NUMB_DENS3D;id++)
      {
          double densc = 0.0;
          for (int ir=0;ir<MC_BINSR;ir++)
          for (int it=0;it<MC_BINST;it++)
          {
           int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
           densc += dens_save[id][ijk];
          }

         fid << setw(IO_WIDTH) << densc/(norm4*(double)MCAtom[id].numb) << BLANK;
      }

      fid << endl;
   }

   fid.close();

// ------  Punch out relative euler angles between the adjacent rotational imaginary time slices  --------------------------
      if(mode == MC_TOTAL)
      {
         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_REP;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double norm5 = (double)NumbRotTimes*acount*_delta_chi*(180./M_PI);

         double phirel;

         for (int ip=0;ip<MC_BINSC;ip++)
         {
            phirel = (double)ip*_delta_chi;
            phirel += (0.5*_delta_chi);

            fid << setw(IO_WIDTH) << (phirel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relphi_sum[ip]/norm5 << endl;

         }

         fid.close();

         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_REC;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double chirel;

         for (int ic=0;ic<MC_BINSC;ic++)
         {
            chirel = (double)ic*_delta_chi;
            chirel += (0.5*_delta_chi);

            fid << setw(IO_WIDTH) << (chirel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relchi_sum[ic]/norm5 << endl;

         }

         fid.close();

         fdens = fname;
         fdens += IO_SUM;
         fdens += IO_EXT_RET;

         fid.open(fdens.c_str(),ios::out); io_setout(fid);

         double norm6 = (double)NumbRotTimes*acount*_delta_theta*(180.0/M_PI);

         double therel;

         for (int it=0;it<MC_BINST;it++)
         {
            therel = (double)it*_delta_theta;
            therel += (0.5*_delta_theta);

            fid << setw(IO_WIDTH) << (therel*180.0/M_PI) << BLANK;

            fid << setw(IO_WIDTH) << _relthe_sum[it]/norm6 << endl;

         }

         fid.close();

      }

   } // END if (IMPURITY) convert 3D into 1D
}

// added by Hui Li
void SaveDensities2D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;
         
  if (IMPURITY)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_DENS2D;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm3 = _delta_radius*_delta_theta*(double)NumbTimes     // the normalization for density 
                  *acount;                                        // distribution around impurity

  double *** dens_save;
  
  dens_save = _gr2D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr2D_sum;


  double theta,r,r2;

  for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
       {
         r   =  (double)ir*_delta_radius;
         r  +=  (0.5*_delta_radius);
         r  +=  _min_radius;
         r2  =   r*r;
    
         for (int id=0;id<NUMB_DENS2D;id++) // 
          {
           fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
           fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
           fid << setw(IO_WIDTH) << dens_save[id][ir][it]/(norm3) << BLANK;
          }

         fid << endl;
        }
   }
   fid.close();
  } //endif 
 } //end subroutine

// added by Toby Zeng
void SaveDensities3D(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;

  if (IMPURITY)  // ? separate from 1D case to avoid long if{}  ?
  {
// ------  1D radial distributions around an impurity --------------------------

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_DENS3D;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm5 = _delta_radius*_delta_theta*_delta_chi*(double)NumbTimes     // the normalization for density 
                  *acount;                                        // distribution around impurity

  double ** dens_save;

  dens_save = _gr3D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr3D_sum;


  double theta,r,r2,chi;

  for (int ir=0;ir<MC_BINSR;ir++)
  {
      r   =  (double)ir*_delta_radius;
      r  +=  (0.5*_delta_radius);
      r  +=  _min_radius;

      for (int it=0;it<MC_BINST;it++)
      {
         theta  = (double)it*_delta_theta;
         theta += (0.5*_delta_theta);

         for (int ic=0;ic<MC_BINSC;ic++)
         {
            chi  = (double)ic*_delta_chi;
            chi += (0.5*_delta_chi);

            for (int id=0;id<NUMB_DENS3D;id++)
            {
               fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
               fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
               fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;
               int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
               fid << setw(IO_WIDTH) << dens_save[id][ijk]/(norm5*r*r*sin(theta)) << BLANK;
            }
            fid<<endl;
         }

      }
  }

/*
  for (int it=0;it<MC_BINST;it++)  // angular distribution
   {
      theta  = (double)it*_delta_theta;
      theta += (0.5*_delta_theta);

      for (int ir=0;ir<MC_BINSR;ir++) // radial distribution
       {
         r   =  (double)ir*_delta_radius;
         r  +=  (0.5*_delta_radius);
         r  +=  _min_radius;
         r2  =   r*r;

         for (int id=0;id<NUMB_DENS2D;id++) // 
          {
           fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
           fid << setw(IO_WIDTH) << (r*Units.length) << BLANK;
           fid << setw(IO_WIDTH) << dens_save[id][ir][it]/(norm5) << BLANK;
          }

         fid << endl;
        }
   }
*/
   fid.close();
  } //endif 
 } //end subroutine

// added by Toby Zeng
void IOxyzAng(int tstatus, const char file_name[])
{
   const char *_proc_=__func__;    // "IOxyz"; 

   string fdens;

//---------------- Open  ------------

   ios::openmode mode;

   switch (tstatus)
   {
      case IOWrite: mode = ios::out;  break;
      case IORead : mode = ios::in;   break;
      default     :
      nrerror (_proc_,IO_ERR_WMODE);  break;
   }

   fdens = file_name;

   if(IOWrite)
   fdens += IO_EXT_XYZ;

   fstream fid(fdens.c_str(),mode);

   if (!fid.good())
   _io_error(_proc_,IO_ERR_FOPEN,file_name);

   io_setout(fid);

//---------------- Read/Write ------------

// stringstream stype;
   string       sbuff,stype;

   int offset;

   int type = 0;
   int atom = 0;  // first atom # will be 1, NOT 0 
   switch (tstatus)
   {
      case IOWrite:
         fid<<MaxnTimes<<" ";                    // total number of "atoms" 
//       permutation table
         if(BOSONS)
         {
            for (int atom=0; atom<MCAtom[BSTYPE].numb;atom++)
            {
             fid <<"  "<<PIndex[atom]<<BLANK;
            }
         }
         fid << endl;

         fid<<COMMENTS<<BLANK<<IO_COM_XYZ<<endl;  // comments

         for (int it=0;it<MaxnTimes;it++)
         {
            if (it==MCAtom[type+1].offset) {type++; atom=0;}    // new atom type
            if ((it-MCAtom[type].offset)%NumbTimes==0) atom++;  // new atom

            fid<<MCAtom[type].type<<atom;      // atom label
//          fid<<setw(5)<<stype<<BLANK;

            for (int id=0;id<NDIM;id++)
            {
            fid<<setw(IO_WIDTH)<<MCCoords[id][it]<<BLANK;
            fid<<setw(IO_WIDTH)<<MCAngles[id][it]<<BLANK;
//          Toby replaces the above line by
//          fid<<setw(IO_WIDTH)<<MCAngles[id][it]<<BLANK;
//          by doing that, Toby stores the three Euler angles, not the unit vector of the two angles orientation
            }
            fid<<endl;
         }

         break;

      case IORead:
         fid>>MaxnTimes;
         getline(fid,sbuff);    // skip a comment line

         for (int type=0;type<NumbTypes;type++)
         for (int atom=0;atom<MCAtom[type].numb;atom++)
         for (int it=0;it<NumbTimes;it++)
         {
            fid>>sbuff;         // skip an atom type

            offset=MCAtom[type].offset+NumbTimes*atom;
            for (int id=0;id<NDIM;id++)
            {
             fid>>MCCoords[id][offset+it];
             fid>>MCCosine[id][offset+it];
//           Toby replaces the above line by
//           fid>>MCAngles[id][offset+it];
//           by doing that, Toby reads the three Euler angles, not the unit vector of the two angle orientation
            }
         }

         break;

      default :
         nrerror (_proc_,IO_ERR_WMODE);
         break;
   }

   fid.close();
}

// added by Toby Zeng
void SaveRhoThetaChi(const char fname [], double acount, int mode)
// the density type corresponds to the atom type
//
//  mode:  MC_TOTAL - accumulated averages
//  mode:  MC_BLOCK - block averages
//
{
  fstream fid;
  string fdens;

  fdens  = fname;

  if (mode == MC_TOTAL)    // accumulated averages
  fdens += IO_SUM;

  fdens += IO_EXT_GTC;

  fid.open(fdens.c_str(),ios::out); io_setout(fid);

  double  norm6 = _delta_theta*_delta_chi*(double)NumbTimes     // the normalization for density 
                  *acount*(180.0*180.0/(M_PI*M_PI));            // distribution around impurity

  double ** dens_save;

  dens_save = _gr3D;
  if (mode == MC_TOTAL)    // accumulated averages density
  dens_save = _gr3D_sum;


  double theta,r,r2,chi;

  for (int it=0;it<MC_BINST;it++)
  {
     theta  = (double)it*_delta_theta;
     theta += (0.5*_delta_theta);

     for (int ic=0;ic<MC_BINSC;ic++)
     {
        chi  = (double)ic*_delta_chi;
        chi += (0.5*_delta_chi);

        fid << setw(IO_WIDTH) << (theta*180.0/M_PI) << BLANK;
        fid << setw(IO_WIDTH) << (chi*180.0/M_PI) << BLANK;

        for (int id=0;id<NUMB_DENS3D;id++)
        {

           double denstc = 0.0;

           for (int ir=0;ir<MC_BINSR;ir++)
           {
              int ijk = (ir*MC_BINST + it)*MC_BINSC + ic;
              denstc = denstc + dens_save[id][ijk];
//            fid << setw(IO_WIDTH) << dens_save[id][ijk]/(norm6) << BLANK;
           }
           fid << setw(IO_WIDTH) << denstc/(norm6*(double)MCAtom[id].numb) << BLANK;
        }
        fid<<endl;
     }

     fid<<endl;

  }

   fid.close();
 } //end subroutine

void GetExchangeLength(void)
{
   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
  _pflags[atom] = 0;

   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   if (_pflags[atom] == 0)
   { 
//    _pflags[patom] = 1;           // do not need if try all atoms in order
 
       int clenght = 0;             // start a new cycle
       int patom   = PIndex[atom]; 
	
       while (patom != atom)
       {
         _pflags[patom] = 1; 
          patom = PIndex[patom];
          clenght++;
       }

      _ploops[clenght] += 1.0;	 
   }
}

void SaveExchangeLength (const char fname [], double acount, long int blocknumb)
{
  const char *_proc_=__func__;    //  SaveExchangeLength() 
 
//------  open file ---------------

  fstream fid;
  string  fperm;

  fperm  = fname;
  fperm += IO_EXT_PRL;

  fid.open(fperm.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fperm.c_str());

//--------------------------------

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

  double excited = 0.0;  
  double ground  = 0.0;   

  for (int clength=0;clength<MCAtom[BSTYPE].numb;clength++)
  {
//   check <norm> below
     double norm  = (double)(clength+1)/(acount*(double)MCAtom[BSTYPE].numb);

     if (clength<=GSLOOP_MAX)  
     excited  += (_ploops[clength]*norm);
     else
     ground   += (_ploops[clength]*norm);
  }

  fid << setw(IO_WIDTH) << ground  << BLANK;   
  fid << setw(IO_WIDTH) << excited << BLANK;  
 
  fid << setw(IO_WIDTH) << (ground+excited) << BLANK; // norm check

// ----------------------------------------------

  for (int clength=0;clength<MCAtom[BSTYPE].numb;clength++)
  {
//   check <norm> above 
     double norm  = (double)(clength+1)/(acount*(double)MCAtom[BSTYPE].numb);

     fid << setw(IO_WIDTH) << _ploops[clength]*norm << BLANK;   
  }
  fid << endl;

//added by Hui Li, test whether it is changed
  for (int atom=0; atom<MCAtom[BSTYPE].numb;atom++)
  {
   fid << setw(IO_WIDTH)<<PIndex[atom]<<BLANK;  
  }

// print PrintXYZprl in the prl file, to specify whether the instantaneous XYZ and PRL are recorded
   fid << PrintXYZprl;

   fid << endl;
// end of test

  fid.close();
}

void GetAreaEstimators(void)
{
   double dr0[NDIM];
   double dr1[NDIM];

   double n_perp[NDIM];
   double n_parl[NDIM];

   double area[NDIM];

   double rn0[NDIM];     // r x n cross product
   double rn1[NDIM];

   double area_perp;
   double area_parl;

   double inert_perp;
   double inert_parl;

// define the ref point for the area estimator
// (i)   instanteneous COM
// (ii)  instanteneous COM of a dopant molecule

// uncomment one section below

/*
//(i)

// double tmass  = 0.0;   // define total mass globally (canonical only)?

// define mass globally change the order of loops

// for (int type=0;type<NumbTypes;type++)
// tmass += (MCAtom[type].mass * (double)MCAtom[type].numb);

   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   {	
      double      tmass  = 0.0;   // define total mass globally (canonical only)?
      newcoords[dim][it] = 0.0;   // center of mass
 
      for (int atom=0;atom<NumbAtoms;atom++)
      {
         int    type = MCType[atom];
         double mass = MCAtom[type].mass;
 
         newcoords[dim][it] += (mass*MCCoords[dim][atom*NumbTimes + it]);	 
         tmass              +=  mass;
      }

      newcoords[dim][it] /= tmass; 
   }
*/

//(ii) only one molecule in the system
   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   newcoords[dim][it] = MCCoords[dim][MCAtom[IMTYPE].offset + it];	 
   int offset = MCAtom[BSTYPE].offset;

   area_perp  = 0.0;
   area_parl  = 0.0;

   inert_perp = 0.0;
   inert_parl = 0.0;

   for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
   for (int it0=0;it0<NumbTimes;it0++)
   {
      int it1 = (it0 + 1) % NumbTimes;

      int pt0 = offset + NumbTimes*atom;     // offset only
      int pt1 = pt0;
 
      if (it1!= (it0 + 1))          
      pt1 = offset + NumbTimes*PIndex[atom]; // offset only  
 
      pt0 += it0;
      pt1 += it1;

      for (int dim=0;dim<NDIM;dim++)         // COM adjustment 
      {
         dr0[dim] = MCCoords[dim][pt0] - newcoords[dim][it0];
         dr1[dim] = MCCoords[dim][pt1] - newcoords[dim][it1];
      }
// ------------- orientations ----------------------------

      int it_rot = it0/RotRatio;

      n_parl[AXIS_X] = MCCosine[AXIS_X][MCAtom[IMTYPE].offset + it_rot];
      n_parl[AXIS_Y] = MCCosine[AXIS_Y][MCAtom[IMTYPE].offset + it_rot];
      n_parl[AXIS_Z] = MCCosine[AXIS_Z][MCAtom[IMTYPE].offset + it_rot];
/*
      n_parl[AXIS_X] = 0.0;
      n_parl[AXIS_Y] = 0.0;
      n_parl[AXIS_Z] = 1.0;
*/
      double static zero = 10e-4; 

      double tg = 0.0;
      double st = 1.0;

      if (fabs(n_parl[AXIS_X]) > zero) // treat separetely n_x=n_y=0 limit ?
      {
         tg = n_parl[AXIS_Y]/n_parl[AXIS_X]; 
         st = sqrt(1.0 + tg*tg);
      }
 
      n_perp[AXIS_X] =  tg /st;
      n_perp[AXIS_Y] = -1.0/st;
      n_perp[AXIS_Z] =  0.0;

/* 
      n_perp[AXIS_X] =  0.0;
      n_perp[AXIS_Y] = -1.0;
      n_perp[AXIS_Z] =  0.0;
*/

//-------------------------------------------------

      area[AXIS_X]  = 0.5*(dr0[AXIS_Y]*dr1[AXIS_Z] - dr0[AXIS_Z]*dr1[AXIS_Y]);
      area[AXIS_Y]  = 0.5*(dr0[AXIS_Z]*dr1[AXIS_X] - dr0[AXIS_X]*dr1[AXIS_Z]);
      area[AXIS_Z]  = 0.5*(dr0[AXIS_X]*dr1[AXIS_Y] - dr0[AXIS_Y]*dr1[AXIS_X]);

      for (int dim=0;dim<NDIM;dim++) // projections
      {
         area_perp += (n_perp[dim]*area[dim]);  
         area_parl += (n_parl[dim]*area[dim]);   	 
      } 

// -- moment of inertia  -----------------------

      rn0[AXIS_X] = n_perp[AXIS_Y]*dr0[AXIS_Z] - n_perp[AXIS_Z]*dr0[AXIS_Y];
      rn0[AXIS_Y] = n_perp[AXIS_Z]*dr0[AXIS_X] - n_perp[AXIS_X]*dr0[AXIS_Z];
      rn0[AXIS_Z] = n_perp[AXIS_X]*dr0[AXIS_Y] - n_perp[AXIS_Y]*dr0[AXIS_X];

      rn1[AXIS_X] = n_perp[AXIS_Y]*dr1[AXIS_Z] - n_perp[AXIS_Z]*dr1[AXIS_Y];
      rn1[AXIS_Y] = n_perp[AXIS_Z]*dr1[AXIS_X] - n_perp[AXIS_X]*dr1[AXIS_Z];
      rn1[AXIS_Z] = n_perp[AXIS_X]*dr1[AXIS_Y] - n_perp[AXIS_Y]*dr1[AXIS_X];

      for (int dim=0;dim<NDIM;dim++) // the perpendicular component
      inert_perp += (rn0[dim]*rn1[dim]);
 
      rn0[AXIS_X] = n_parl[AXIS_Y]*dr0[AXIS_Z] - n_parl[AXIS_Z]*dr0[AXIS_Y];
      rn0[AXIS_Y] = n_parl[AXIS_Z]*dr0[AXIS_X] - n_parl[AXIS_X]*dr0[AXIS_Z];
      rn0[AXIS_Z] = n_parl[AXIS_X]*dr0[AXIS_Y] - n_parl[AXIS_Y]*dr0[AXIS_X];

      rn1[AXIS_X] = n_parl[AXIS_Y]*dr1[AXIS_Z] - n_parl[AXIS_Z]*dr1[AXIS_Y];
      rn1[AXIS_Y] = n_parl[AXIS_Z]*dr1[AXIS_X] - n_parl[AXIS_X]*dr1[AXIS_Z];
      rn1[AXIS_Z] = n_parl[AXIS_X]*dr1[AXIS_Y] - n_parl[AXIS_Y]*dr1[AXIS_X];

      for (int dim=0;dim<NDIM;dim++)  // the parallel component
      inert_parl += (rn0[dim]*rn1[dim]);
   }

  _areas[PERP] += area_perp;   
  _areas[PARL] += area_parl; 
 
  _area2[PERP] += (area_perp*area_perp);   
  _area2[PARL] += (area_parl*area_parl);  

  _inert[PERP] += (inert_perp/(double)NumbTimes);
  _inert[PARL] += (inert_parl/(double)NumbTimes);
}

void GetAreaEstim3D(int iframe)
{

// iframe specifies the frame on which the superfluid response is projected 0: Space-fixed frame; 1: dopant-fixed frame

   double dr0[NDIM];     // the position vector of a bead with respect to the COM of dopant
   double dr1[NDIM];     // the position vector of the bead one slice after

   double inert3D[NDIM*NDIM];

   double rn0[NDIM];     // dr0 x n cross product
   double rn1[NDIM];     // dr1 x n cross product

   double area[NDIM];    // the vector area of dr0 cross dr1

   double area_proj[NDIM]; // projection of the area vector onto the three principle axes of the non-linear dopant

// define the ref point for the area estimator
// (i)   instanteneous COM
// (ii)  instanteneous COM of a dopant molecule

// uncomment one section below

//(i)

// double tmass  = 0.0;   // define total mass globally (canonical only)?

// define mass globally change the order of loops

// for (int type=0;type<NumbTypes;type++)
// tmass += (MCAtom[type].mass * (double)MCAtom[type].numb);
   if(iframe == 0)
   {
///*
   for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
   for (int it=0;it<NumbTimes;it++)
   {    
      double      tmass  = 0.0;   // define total mass globally (canonical only)?
      newcoords[dim][it] = 0.0;   // center of mass
 
      for (int atom=0;atom<NumbAtoms;atom++)
      {
         int    type = MCType[atom];
         double mass = MCAtom[type].mass;
 
         newcoords[dim][it] += (mass*MCCoords[dim][atom*NumbTimes + it]);        
         tmass              +=  mass;
      }

      newcoords[dim][it] /= tmass; 
   }
//*/
   }
//(ii) only one molecule in the system

   if(iframe == 1)
   {
      for (int dim=0;dim<NDIM;dim++) // instantaneous center of mass
      for (int it=0;it<NumbTimes;it++)
      newcoords[dim][it] = MCCoords[dim][MCAtom[IMTYPE].offset + it];
   }


   int offset = MCAtom[BSTYPE].offset;

   double bmass = MCAtom[BSTYPE].mass;

// clean area_proj and inert3D
   for (int id=0;id<NDIM;id++)
   {
      area_proj[id]=0.0;
      for (int jd=0;jd<NDIM;jd++)
      {
         int ij = id *NDIM + jd;
         inert3D[ij]=0.0;
      }
   }

// private variables that accumulate results in each cpu
   double area_projx=0.0,area_projy=0.0,area_projz=0.0,icxx=0.0,icxy=0.0,icxz=0.0,icyx=0.0,icyy=0.0,icyz=0.0,iczx=0.0,iczy=0.0,iczz=0.0;

// compared to the original GetAreaEstimators, here Toby switches the order of loop over it0 and loop over atom.
// This is to minimize the call of vcord
   #pragma omp parallel for reduction(+: area_projx,area_projy,area_projz,icxx,icxy,icxz,icyx,icyy,icyz,iczx,iczy,iczz) private(area,rn0,rn1,dr0,dr1)
   for (int it0=0;it0<NumbTimes;it0++)
   {

      int it_rot = it0/RotRatio + MCAtom[IMTYPE].offset; // only one molecular impurity

      double RCOM[3];
      double Rpt[3];
      double Eulang[3];
      double vpot3d;
      double radret;
      double theret;
      double chiret;
      double hatx[3];
      double haty[3];
      double hatz[3];
      int    ivcord = 1;

      for (int id=0;id<NDIM;id++)
      {
//       Rpt[id]  = MCCoords[id][it_rot];
//       RCOM[id] = MCCoords[id][it_rot];
//       orientation of molecule should not dipend on the position of mass of centre
         Rpt[id]  = 0.0;
         RCOM[id] = 0.0;
      }

      Eulang[PHI]=MCAngles[PHI][it_rot];
      Eulang[CTH]=acos(MCAngles[CTH][it_rot]);
      Eulang[CHI]=MCAngles[CHI][it_rot];

      vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

/*
// get the orientation of the dopant in the next time slice
      int it_rot2 = ((it0 + 1) % NumbTimes)/RotRatio;

      it_rot2 += MCAtom[IMTYPE].offset;

      double hatx2[3];
      double haty2[3];
      double hatz2[3];

      Eulang[PHI]=MCAngles[PHI][it_rot2];
      Eulang[CTH]=acos(MCAngles[CTH][it_rot2]);
      Eulang[CHI]=MCAngles[CHI][it_rot2];

      vcord_(Eulang,RCOM,Rpt,vtable,&vpot3d,&radret,&theret,&chiret,hatx2,haty2,hatz2,&ivcord);
*/

/*
//    check hatz with MCCosine
      for(int id=0;id<NDIM;id++)
      {
         if(fabs(hatz[id] - MCCosine[id][it_rot] > 0.00005))
         cout<<id<<" "<<it_rot<<" "<<hatz[id]<<" "<<MCCosine[id][it_rot]<<endl;
      }
*/
      
//    SFF axes
      if(iframe == 0)
      {
         hatx[0]=1.0,hatx[1]=0.0,hatx[2]=0.0;
         haty[0]=0.0,haty[1]=1.0,haty[2]=0.0;
         hatz[0]=0.0,hatz[1]=0.0,hatz[2]=1.0;
      }

      for (int atom=0;atom<MCAtom[BSTYPE].numb;atom++)
      {

         int it1 = (it0 + 1) % NumbTimes;

         int pt0 = offset + NumbTimes*atom;     // offset only
         int pt1 = pt0;

         if (it1!= (it0 + 1))
         pt1 = offset + NumbTimes*PIndex[atom]; // offset only

         pt0 += it0;
         pt1 += it1;

         for (int dim=0;dim<NDIM;dim++)         // COM adjustment 
         {
            dr0[dim] = MCCoords[dim][pt0] - newcoords[dim][it0];
            dr1[dim] = MCCoords[dim][pt1] - newcoords[dim][it1];
         }

//       get the x0, y0, z0, x1, y1, and z1 in their respective dopant-fixed frame
/*
         double x0=0.0;
         double y0=0.0;
         double z0=0.0;
         double x1=0.0;
         double y1=0.0;
         double z1=0.0;
         for (int dim=0;dim<NDIM;dim++)
         {
            x0 += hatx[dim]*dr0[dim];
            y0 += haty[dim]*dr0[dim];
            z0 += hatz[dim]*dr0[dim];
            x1 += hatx2[dim]*dr1[dim];
            y1 += haty2[dim]*dr1[dim];
            z1 += hatz2[dim]*dr1[dim];
         }

         area_proj[AXIS_X] += 0.5*(y0*z1 - z0*y1);
         area_proj[AXIS_Y] += 0.5*(z0*x1 - x0*z1);
         area_proj[AXIS_Z] += 0.5*(x0*y1 - y0*x1);
*/
//----------- cross product of the two adjacent position vectors with respect to the dopant COM to tet area --------------------------------------

///*
         area[AXIS_X]  = 0.5*(dr0[AXIS_Y]*dr1[AXIS_Z] - dr0[AXIS_Z]*dr1[AXIS_Y]);
         area[AXIS_Y]  = 0.5*(dr0[AXIS_Z]*dr1[AXIS_X] - dr0[AXIS_X]*dr1[AXIS_Z]);
         area[AXIS_Z]  = 0.5*(dr0[AXIS_X]*dr1[AXIS_Y] - dr0[AXIS_Y]*dr1[AXIS_X]);

//       project the vector area onto the three principal axes of the dopant
         for (int id=0;id<NDIM;id++)
         {
            area_projx += area[id] * hatx[id];
            area_projy += area[id] * haty[id];
            area_projz += area[id] * hatz[id];
         }
//*/

///*     calculate the classical moment of inertia
         for (int id=0;id<NDIM;id++)
         {
//          diagonal elements first
            double *hat_dum;
            if(id == AXIS_X)
            hat_dum = hatx;

            if(id == AXIS_Y)
            hat_dum = haty;

            if(id == AXIS_Z)
            hat_dum = hatz;

//          cross product dr0 x hat(x,y,z)
            rn0[AXIS_X] = dr0[AXIS_Y]*hat_dum[AXIS_Z] - dr0[AXIS_Z]*hat_dum[AXIS_Y];
            rn0[AXIS_Y] = dr0[AXIS_Z]*hat_dum[AXIS_X] - dr0[AXIS_X]*hat_dum[AXIS_Z];
            rn0[AXIS_Z] = dr0[AXIS_X]*hat_dum[AXIS_Y] - dr0[AXIS_Y]*hat_dum[AXIS_X];

//          cross product dr1 x hat(x,y,z)
            rn1[AXIS_X] = dr1[AXIS_Y]*hat_dum[AXIS_Z] - dr1[AXIS_Z]*hat_dum[AXIS_Y];
            rn1[AXIS_Y] = dr1[AXIS_Z]*hat_dum[AXIS_X] - dr1[AXIS_X]*hat_dum[AXIS_Z];
            rn1[AXIS_Z] = dr1[AXIS_X]*hat_dum[AXIS_Y] - dr1[AXIS_Y]*hat_dum[AXIS_X];

//          accumulate the dot product of rn0 and rn1 into the diagonal element of classical moment of inertia
            double sum=0.0;
            for (int dim=0;dim<NDIM;dim++)
            sum += rn0[dim]*rn1[dim]*bmass;
//          inert3D[id*NDIM+id] += rn0[dim]*rn1[dim]*bmass;

            if(id == AXIS_X)
            icxx += sum;

            if(id == AXIS_Y)
            icyy += sum;

            if(id == AXIS_Z)
            iczz += sum;

//          off diagonal elements below
//          projection of dr0 on hat_dum (id)
            double dr0_id = 0.0;
            for (int dim=0;dim<NDIM;dim++)
            dr0_id +=hat_dum[dim]*dr0[dim];

            for (int jd=0;jd<NDIM;jd++)
            if( jd != id )
            {
               double *hat_dum2;
               if(jd == AXIS_X)
               hat_dum2 = hatx;

               if(jd == AXIS_Y)
               hat_dum2 = haty;

               if(jd == AXIS_Z)
               hat_dum2 = hatz;

               double dr1_jd = 0.0;
               for (int dim=0;dim<NDIM;dim++)
               dr1_jd +=hat_dum2[dim]*dr1[dim];

               //inert3D[id*NDIM + jd] += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_X && jd == AXIS_Y)
               icxy += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_X && jd == AXIS_Z)
               icxz += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Y && jd == AXIS_X)
               icyx += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Y && jd == AXIS_Z)
               icyz += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Z && jd == AXIS_X)
               iczx += -bmass * dr0_id * dr1_jd;

               if(id == AXIS_Z && jd == AXIS_Y)
               iczy += -bmass * dr0_id * dr1_jd;

            }

         }
//*/

      }
   }

// transfer the parallelly accumulated results in the corresponding arrays
   area_proj[AXIS_X] = area_projx;
   area_proj[AXIS_Y] = area_projy;
   area_proj[AXIS_Z] = area_projz;
   inert3D[AXIS_X*NDIM+AXIS_X] = icxx;
   inert3D[AXIS_X*NDIM+AXIS_Y] = icxy;
   inert3D[AXIS_X*NDIM+AXIS_Z] = icxz;
   inert3D[AXIS_Y*NDIM+AXIS_X] = icyx;
   inert3D[AXIS_Y*NDIM+AXIS_Y] = icyy;
   inert3D[AXIS_Y*NDIM+AXIS_Z] = icyz;
   inert3D[AXIS_Z*NDIM+AXIS_X] = iczx;
   inert3D[AXIS_Z*NDIM+AXIS_Y] = iczy;
   inert3D[AXIS_Z*NDIM+AXIS_Z] = iczz;

// scale the inert3D by 1/NumbTimes and add it to the block accumulation
   for (int id=0;id<NDIM*NDIM;id++)
   {
      if(iframe == 1)
      _inert3DMFF[id] += inert3D[id]/(double)NumbTimes;

      if(iframe == 0)
      _inert3DSFF[id] += inert3D[id]/(double)NumbTimes;

   }

// add area_proj to the block accumulation
   int ind = 0;
   for (int id=0;id<NDIM;id++)
   {
//    cout<<area_proj[id]<<BLANK;
      for (int jd=0;jd<=id;jd++)
      {
         if(iframe == 1)
         _areas3DMFF[ind] += area_proj[id]*area_proj[jd];

         if(iframe == 0)
         _areas3DSFF[ind] += area_proj[id]*area_proj[jd];


         ind ++;
      }
   }
// cout<<endl;

}

void SaveAreaEstimators (const char fname [], double acount, long int blocknumb)
{
  const char *_proc_=__func__;    //  SaveAreaDensities() 
 
//------  open file ---------------

  fstream fid;
  string  fsuper;

  fsuper  = fname;
  fsuper += IO_EXT_SUP;

  fid.open(fsuper.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fsuper.c_str());

//--------------------------------
//  shift the center of mass 
//--------------------------------

  double mass   = MCAtom[BSTYPE].mass; 
  double lambda = MCAtom[BSTYPE].lambda; 

  double norm   = 2.0/(MCBeta*lambda);

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

// super density (normalization: inertia does not include mass)

  fid << setw(IO_WIDTH) << _area2[PERP]*norm/_inert[PERP] << BLANK; 
  fid << setw(IO_WIDTH) << _area2[PARL]*norm/_inert[PARL] << BLANK;

// moment of inertia

  double units = Units.mass*Units.length*Units.length;

  fid << setw(IO_WIDTH) << _inert[PERP]*units*mass/acount << BLANK;
  fid << setw(IO_WIDTH) << _inert[PARL]*units*mass/acount << BLANK;

// area (normalized and averaged)

  fid << setw(IO_WIDTH) << _areas[PERP]*sqrt(norm/_inert[PERP])/acount << BLANK; 
  fid << setw(IO_WIDTH) << _areas[PARL]*sqrt(norm/_inert[PARL])/acount << BLANK;
 
  fid << endl;
  fid.close();
}

void super_reset (void)
{
   _areas[PERP] = 0.0;
   _areas[PARL] = 0.0;

   _area2[PERP] = 0.0;
   _area2[PARL] = 0.0;

   _inert[PERP] = 0.0;
   _inert[PARL] = 0.0;

   for (int id=0;id<6;id++)
   {
      _areas3DMFF[id] = 0.0;
      _areas3DSFF[id] = 0.0;
   }

   for (int id=0;id<NDIM*NDIM;id++)
   {
      _inert3DMFF[id]=0.0;
      _inert3DSFF[id]=0.0;
   }

}

void SaveAreaEstim3D (const char fname [], double acount, long int blocknumb, int iframe)
{

  const char *_proc_=__func__;    //  SaveAreaDensities()

//------  open file ---------------

  fstream fid;
  string  fsuper;

  fsuper  = fname;

  if(iframe == 1)
  fsuper += IO_EXT_MFFSUP3D;

  if(iframe == 0)
  fsuper += IO_EXT_SFFSUP3D;

  fid.open(fsuper.c_str(),ios::app | ios::out); io_setout(fid);

  if (!fid.is_open())
 _io_error(_proc_,IO_ERR_FOPEN,fsuper.c_str());

  double bmass = MCAtom[BSTYPE].mass;
  double lambda = MCAtom[BSTYPE].lambda;

  double norm = 2.0*bmass/(MCBeta * lambda);  // norm = 4m^2/(hbar^2 * beta)

  fid << setw(IO_WIDTH_BLOCK) << blocknumb << BLANK;     // block number

// 9 components of classical moment of inertia
  for(int id=0;id<NDIM*NDIM;id++)
  {
     if(iframe == 1)
     fid << setw(IO_WIDTH) << _inert3DMFF[id]/(double)acount<< BLANK;

     if(iframe == 0)
     fid << setw(IO_WIDTH) << _inert3DSFF[id]/(double)acount<< BLANK;
  }

// 6 components of 4m^2/(hbar^2*beta) <A_iA_j>
   for(int id=0;id<6;id++)
   {
      if(iframe == 1)
      fid << setw(IO_WIDTH) << _areas3DMFF[id]*norm/(double)acount<< BLANK;

      if(iframe == 0)
      fid << setw(IO_WIDTH) << _areas3DSFF[id]*norm/(double)acount<< BLANK;
   }

// 3 compponents of area vector
// for(int id=0;id<NDIM;id++)
// fid << setw(IO_WIDTH) << _area3D[id]/(double)acount <<BLANK;

// fid << setw(IO_WIDTH) << norm <<BLANK;

  fid << endl;
  fid.close();

}

void GetPermutation()
{

   const char *_proc_=__func__;    //  SaveAreaDensities()

   fstream fid;
   string fname;
   fname = FPERMU;

   fid.open(fname.c_str(),ios::app | ios::out); io_setout(fid);

   if (!fid.is_open())
   _io_error(_proc_,IO_ERR_FOPEN,fname.c_str());

   for(int iat=0;iat<MCAtom[BSTYPE].numb;iat++)
   {
      fid<<PIndex[iat]<<" ";
   }

   fid<<endl;
   fid.close();
}
//Below the routines are used for the calculation of dot product and cross product between two vectors ---  added by Tapas Sahoo
void VectorNormalisation(double *v)
{
    double lenght = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    v[0] /= lenght;
    v[1] /= lenght;
    v[2] /= lenght;
}

double DotProduct(double *v, double *w)
{
    return (v[0] * w[0] + v[1] * w[1] + v[2] * w[2]);
}

void CrossProduct(double *v, double *w, double *cross)
{
    cross[0] = w[1] * v[2] - w[2] * v[1];
    cross[1] = w[2] * v[0] - w[0] * v[2];
    cross[2] = w[0] * v[1] - w[1] * v[0];
}

// PotFunc() for analytic potential

double PotFunc(int atom0, int atom1, const double *Eulang0, const double *Eulang1, int it)
{
	int offset0 = NumbRotTimes*atom0;
	int offset1 = NumbRotTimes*atom1;
   	int t0 = offset0 + it;
   	int t1 = offset1 + it;

	double DipoleMomentInAU = DipoleMoment/AuToDebye; // DipoleMoment in Debye

	double R12[NDIM];
	double dr2 = 0.0;
	double DipoleMoment0[NDIM], DipoleMoment1[NDIM];
    for (int id = 0; id < NDIM; id++)
	{
		DipoleMoment0[id] = 0.0;
		DipoleMoment1[id] = 0.0;
        R12[id]  = (MCCoords[id][t1] - MCCoords[id][t0]);
		R12[id] /= BOHRRADIUS;
        dr2     += (R12[id]*R12[id]);
	}

	double RCOM = sqrt(dr2);	
	double R2   = RCOM*RCOM;
	double R5   = R2*R2*RCOM;

#ifdef SHORTFORM
	double R3         = R2*RCOM; 
    double PreFactor  = DipoleMomentInAU*DipoleMomentInAU/R3;
    double potential = PreFactor*(sin(Eulang0[CTH])*sin(Eulang1[CTH])*cos(Eulang0[PHI] - Eulang1[PHI]) - 2.0*cos(Eulang0[CTH])*cos(Eulang1[CTH]));
#else
	double dm[NDIM];
	dm[0] = 0.0;
	dm[1] = 0.0;
	dm[2] = DipoleMomentInAU;

	double RotMat0[NDIM*NDIM];
	for (int i = 0; i < (NDIM*NDIM); i++) RotMat0[i] = 0.0;
	UnitVectors(Eulang0, RotMat0);

	double RotMat1[NDIM*NDIM];
	for (int i = 0; i < (NDIM*NDIM); i++) RotMat1[i] = 0.0;
	UnitVectors(Eulang1, RotMat1);

	double delta[NDIM*NDIM];
    for (int i = 0; i < NDIM; i++)
	{
		for (int j = 0; j < NDIM; j++)
		{
			int jj = j + i*NDIM;
			DipoleMoment0[i] += RotMat0[jj]*dm[j];
			DipoleMoment1[i] += RotMat1[jj]*dm[j];

			if (i == j)
            {
                delta[jj] = 1.0;
            }
            else
            {
                delta[jj] = 0.0;
            }
		}
	}

	double potential = 0.0;
	for (int i = 0; i < NDIM; i++)
	{
		for (int j = 0; j < NDIM; j++)
		{
			int jj = j + i*NDIM;
        	potential += - DipoleMoment0[i]*DipoleMoment1[j]*(3.0*R12[i]*R12[j] - R2*delta[jj])/R5;
		}
	}
#endif
    double PotReturn = potential*AuToKelvin;
#ifdef POTZERO
	PotReturn = 0.0;
#endif
    return PotReturn;
}

void UnitVectors(const double *Eulang, double *RotMat)
{
	double theta = Eulang[CTH];
	double phi   = Eulang[PHI];
	double chi   = Eulang[CHI];

	double cp = cos(phi);
    double sp = sin(phi);
    double ct = cos(theta);
    double st = sin(theta);
    double ck = cos(chi);
    double sk = sin(chi);

    RotMat[0] = cp*ct*ck-sp*sk;
    RotMat[1] = -cp*ct*sk-sp*ck;
    RotMat[2] = cp*st;
	RotMat[3] = sp*ct*ck+cp*sk;
    RotMat[4] = -sp*ct*sk+cp*ck;
    RotMat[5] = sp*st;
    RotMat[6] = -st*ck;
    RotMat[7] = st*sk;
    RotMat[8] = ct;
}

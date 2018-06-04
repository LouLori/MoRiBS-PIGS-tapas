// LIMITATION :  only one atom type and one molecule type at this point

#include "mc_input.h"
#include "mc_confg.h"
#include "mc_const.h"
#include "mc_setup.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_piqmc.h"
#include "mc_qworm.h"
#include "mc_estim.h"
//#include <mpi.h>
#include <omp.h>
#include <string>

#include <math.h>
#include "rngstream.h"
#include "omprng.h"

// counters

double ** MCTotal;   // MC counters (total number of moves)
double ** MCAccep;   // MC counters (number of accepted moves)

// counters for parallel 3d rotation move.  They are supposed to be declared for all CPUs
double MCRotChunkAcp;    // total number of rotational moves for one chunk loop
double MCRotChunkTot;    // total accept number of rotational moves for one chunk loop
double MCRotTot; // the sum of all MCRotChunkTot from all CPUs
double MCRotAcp; // the sum of all MCRotChunkAcp from all CPUs

extern "C" void rotden_(double *Eulan1,double *Eulan2,double *Eulrel,double *rho,double *erot,double *esq,double *rhoprp, double *erotpr,double *erotsq,int *istop); // external fortran subroutine by Toby

extern "C" void vcord_(double *Eulang, double *RCOM, double *Rpt, double *vtable, int *Rgrd, int *THgrd, int *CHgrd, double *Rvmax, double *Rvmin, double *Rvstep, double *vpot3d, double *radret, double *theret, double *chiret, double *hatx, double *haty, double *hatz, int *ivcord);

extern "C" void rsrot_(double *Eulan1,double *Eulan2,double *X_Rot,double *Y_Rot,double *Z_Rot,double *tau,int *iodevn,double *eoff,double *rho,double *erot);

extern "C" void rsline_(double *X_Rot,double *p0,double *tau,double *rho,double *erot);

extern "C" void vspher_(double *r,double *vpot);

extern "C" void reflec_(double *coord,double *rcom,double *hatx,double *haty,double *hatz);

extern "C" void rflmfx_(double *RCOM,double *hatx,double *haty,double *hatz,double *Eulang);
extern "C" void rflmfy_(double *RCOM,double *hatx,double *haty,double *hatz,double *Eulang);
extern "C" void rflmfz_(double *RCOM,double *hatx,double *haty,double *hatz,double *Eulang);

// GG ---> potentiel H2O ---- H2O
extern "C" void caleng_(double *com_1, double *com_2, double *E_2H2O, double *Eulang_1, double *Eulang_2);
// Hinde potential for H2 - H2
extern "C" void vh2h2_(double *rd, double *r1, double *r2, double *t1, double *t2, double *phi, double *potl);

int PrintYrfl; // integer flag for printing reflected coordinates
int PrintXrfl; // integer flag for printing reflected coordinates
int PrintZrfl; // integer flag for printing reflected coordinates

void MCMolecularMove(int type)
{
	int numb = MCAtom[type].numb;  

  	double disp[NDIM];

  	for (int atom = 0; atom < numb; atom++)
  	{
    	int offset = MCAtom[type].offset + NumbTimes*atom;
    	int gatom  = offset/NumbTimes;

    	for (int id = 0; id < NDIM; id++)	  // MOVE
    	disp[id] = MCAtom[type].mcstep*(rnd1()-0.5);

    	for (int id = 0; id < NDIM; id++)	  // MOVE
    	{
    		#pragma omp parallel for
    		for (int it = 0; it < NumbTimes; it++)
    		{ 
       			newcoords[id][offset+it]  =  MCCoords[id][offset+it];
       			newcoords[id][offset+it] +=  disp[id];
    		}
    	}

    	double deltav = 0.0;         // ACCEPT/REJECT
    
    	deltav += (PotEnergy(gatom,newcoords)-PotEnergy(gatom,MCCoords));

    	bool Accepted = false;

    	if (deltav<0.0)             Accepted = true;
    	else if
    	(exp(-deltav*MCTau)>rnd2()) Accepted = true;

    	MCTotal[type][MCMOLEC] += 1.0;  
      
    	if (Accepted)
    	{
       		MCAccep[type][MCMOLEC] += 1.0; 

       		for (int id = 0; id < NDIM; id++)       // save accepted configuration	
       		{
       			#pragma omp parallel for
       			for (int it = 0; it < NumbTimes; it++)
       			MCCoords[id][offset+it] = newcoords[id][offset+it];
       		}
    	}	     
  	}   
}

#ifdef GAUSSIANMOVE
void MCMolecularMoveGauss(int type)
{
	int numb = MCAtom[type].numb;  

	GetRandomCoords();
  	for (int atom = 0; atom < numb; atom++)
  	{
    	int offset = MCAtom[type].offset + NumbTimes*atom;
    	int gatom  = offset/NumbTimes;

    	for (int id = 0; id < NDIM; id++)	  // MOVE
    	{
    		//#pragma omp parallel for
    		for (int it = 0; it < NumbTimes; it++)
    		{ 
   				MCCoords[id][offset+it] = gausscoords[id][offset+it];

       		}
    	}
  	}
}
#endif

void MCMolecularMoveExchange(int type)
// particular atom type
{
#ifdef DEBUG_PIMC
   const char *_proc_ = __func__;  // MCMolecularMoveExchange(int)
 
   if (type != BSTYPE)
   nrerror(_proc_,"Wrong atom type");
#endif

  bool Accepted;
  double disp[NDIM];

  int numb = MCAtom[type].numb; 
 
  for (int atom=0;atom<numb;atom++)
 _pflags[atom] = 0;

  for (int atom=0;atom<numb;atom++)
  if (_pflags[atom] == 0)           // start a new cycle
  {
    int amax = 0;                   // number of atoms in the current
    atom_list[amax] = atom;         // exchange loop  
    amax ++;

    int patom   = PIndex[atom]; 
	
    while (patom != atom)         // count all atoms involved  
    {                             // in current exchange loop  
       atom_list[amax] = patom; 
       amax ++;
      _pflags[patom] = 1; 

       patom = PIndex[patom];
    }

    for (int id=0;id<NDIM;id++)      // MOVE
    disp[id] = MCAtom[type].mcstep*(rnd1()-0.5);

    for (int ia=0;ia<amax;ia++)     // loop over atoms in the current
    {                               // exchange loop
       int offset = MCAtom[type].offset + NumbTimes*atom_list[ia];
       int gatom  = offset/NumbTimes;


       for (int id=0;id<NDIM;id++)   // MOVE
       {
       #pragma omp parallel for
       for (int it=0;it<NumbTimes;it++)
       { 
          newcoords[id][offset+it]  =  MCCoords[id][offset+it];
          newcoords[id][offset+it] +=  disp[id];
       }
       }

       double deltav = 0.0;         // ACCEPT/REJECT
    
       deltav += (PotEnergy(gatom,newcoords)-PotEnergy(gatom,MCCoords));

       Accepted = false;

       if (deltav<0.0)            Accepted = true;
       else if
      (exp(-deltav*MCTau)>rnd2()) Accepted = true;

       if (!Accepted) break;

     } // END sum over atoms in the current exchange loop

     MCTotal[type][MCMOLEC] += 1.0;  
      
     if (Accepted)
     {
        MCAccep[type][MCMOLEC] += 1.0;
 
        for (int ia=0;ia<amax;ia++)       // loop over atoms in the current
        {                             // exchange loop
            int offset = MCAtom[type].offset + NumbTimes*atom_list[ia];
 
            for (int id=0;id<NDIM;id++)       // save accepted configuration	
            {
            #pragma omp parallel for
            for (int it=0;it<NumbTimes;it++)
            MCCoords[id][offset+it] = newcoords[id][offset+it];
            }
        } 
     }
   }    // END sum over atoms (fixed atom type)
}

void MCBisectionMove(int type, int time)  // multilevel Metropolis
{
	int numb = MCAtom[type].numb;

   	double mclambda = MCAtom[type].lambda;    
   	int    mclevels = MCAtom[type].levels;  // number of levels
   	int    seg_size = MCAtom[type].mlsegm;  // segmen size  

   	for (int atom=0;atom<numb;atom++)         // one atom to move only
   	{

      	int offset = MCAtom[type].offset + NumbTimes*atom;
      	int gatom  = offset/NumbTimes;

// initialize the end points

      	int pit = (time+seg_size) % NumbTimes;  // periodicity in time 	       	
      	for (int id=0;id<NDIM;id++)             
      	{	  
         	newcoords[id][offset + time] = MCCoords[id][offset + time];
         	newcoords[id][offset + pit]  = MCCoords[id][offset + pit];
      	}

      	double bnorm = 1.0/(mclambda*MCTau);  // variance for gaussian sampling 

      	bool Accepted; 

      	int t0,t1,t2;

      	double pot0 = 0.0;  // potential, current  level
      	double pot1 = 0.0;  // potential, previous level
       	
      	for (int level=0;level<mclevels;level++) // loop over bisection levels
      	{	                                                          
         	int level_seg_size = (int)pow(2.0,(mclevels-level));

         	double bkin_norm = bnorm/(double) level_seg_size;
         	double bpot_norm = MCTau*(double)(level_seg_size/2);
	   
         	pot1 = pot0;   // swap level potentials
         	pot0 = 0.0;
	   
         	t2 = 0;
         	do             // loop over middle points
         	{
            	t0 =  t2;                   // left point
            	t2 =  t0 + level_seg_size;  // right point	
            	t1 = (t0 + t2)/2;           // middle point

            	int pt0 = (time + t0) % NumbTimes;	
            	int pt1 = (time + t1) % NumbTimes;	
            	int pt2 = (time + t2) % NumbTimes;	

//  change the offset if exchange
 
            	for (int id=0;id<NDIM;id++)
            	{  	   
               		newcoords[id][offset+pt1]  = 0.5*(newcoords[id][offset+pt0]+newcoords[id][offset+pt2]);
               		newcoords[id][offset+pt1] += gauss(bkin_norm);
            	} 
//---------------------------- the end point approximation    

            	pot0 += (PotEnergy(gatom,newcoords,pt1) - PotEnergy(gatom,MCCoords,pt1));
  
            	if (t0!=0)                // skip the contributions of the end points
            	pot0 += (PotEnergy(gatom,newcoords,pt0) - PotEnergy(gatom,MCCoords,pt0));
         	}   	      
         	while (t2<seg_size);        // end the loop over middle points 

// inefficient version

         	double deltav = (pot0-2.0*pot1);  // rho(0,1;tau) 
         	deltav *= bpot_norm;
 
         	Accepted = false;
       
         	if (deltav<0.0)               Accepted = true;
         	else if (exp(-deltav)>rnd3()) Accepted = true;

         	if (!Accepted) break;

     	}  // END loop over levels        

     	MCTotal[type][MCMULTI] += 1.0;
     
     	if (Accepted)     
     	{
         	MCAccep[type][MCMULTI] += 1.0;
 
         	for (int id=0;id<NDIM;id++)                // save new coordinates
         	for (int it=time;it<=(time+seg_size);it++)    
         	{  
            	int pit = it % NumbTimes;                // periodicity in time 	       	
            	MCCoords[id][offset+pit] = newcoords[id][offset+pit];
         	}                                                           
      	}	     
//-----------------------------------------------------------------------  
//      END bisection 
//-----------------------------------------------------------------------
  	}  // END loop over time slices/atoms
}

void MCBisectionMovePIGS(int type, int time)  // multilevel Metropolis
{
	int numb = MCAtom[type].numb;

   	double mclambda = MCAtom[type].lambda;    
   	int    mclevels = MCAtom[type].levels;  // number of levels
   	int    seg_size = MCAtom[type].mlsegm;  // segmen size  

   	for (int atom=0;atom<numb;atom++)         // one atom to move only
   	{

      	int offset = MCAtom[type].offset + NumbTimes*atom;
      	int gatom  = offset/NumbTimes;

// initialize the end points

      	int pit = (time+seg_size) % NumbTimes;  // periodicity in time 	       	
      	for (int id=0;id<NDIM;id++)             
      	{	  
         	newcoords[id][offset + time] = MCCoords[id][offset + time];
         	newcoords[id][offset + pit]  = MCCoords[id][offset + pit];
      	}

      	double bnorm = 1.0/(mclambda*MCTau);  // variance for gaussian sampling 

      	bool Accepted; 

      	int t0,t1,t2;

      	double pot0 = 0.0;  // potential, current  level
      	double pot1 = 0.0;  // potential, previous level
       	
      	for (int level=0;level<mclevels;level++) // loop over bisection levels
      	{	                                                          
         	int level_seg_size = (int)pow(2.0,(mclevels-level));

         	double bkin_norm = bnorm/(double) level_seg_size;
         	double bpot_norm = MCTau*(double)(level_seg_size/2);
	   
         	pot1 = pot0;   // swap level potentials
         	pot0 = 0.0;
	   
         	t2 = 0;
         	do             // loop over middle points
         	{
            	t0 =  t2;                   // left point
            	t2 =  t0 + level_seg_size;  // right point	
            	t1 = (t0 + t2)/2;           // middle point

            	int pt0 = (time + t0) % NumbTimes;	
            	int pt1 = (time + t1) % NumbTimes;	
            	int pt2 = (time + t2) % NumbTimes;	

//  change the offset if exchange
 
            	for (int id=0;id<NDIM;id++)
            	{  	   
               		newcoords[id][offset+pt1]  = 0.5*(newcoords[id][offset+pt0]+newcoords[id][offset+pt2]);
               		newcoords[id][offset+pt1] += gauss(bkin_norm);
            	} 
//---------------------------- the end point approximation    

            	pot0 += (PotEnergy(gatom,newcoords,pt1) - PotEnergy(gatom,MCCoords,pt1));
  
            	if (t0!=0)                // skip the contributions of the end points
            	pot0 += (PotEnergy(gatom,newcoords,pt0) - PotEnergy(gatom,MCCoords,pt0));
         	}   	      
         	while (t2<seg_size);        // end the loop over middle points 

// inefficient version

         	double deltav = (pot0-2.0*pot1);  // rho(0,1;tau) 
         	deltav *= bpot_norm;
 
         	Accepted = false;
       
         	if (deltav<0.0)               Accepted = true;
         	else if (exp(-deltav)>rnd3()) Accepted = true;

         	if (!Accepted) break;

     	}  // END loop over levels        

     	MCTotal[type][MCMULTI] += 1.0;
     
     	if (Accepted)     
     	{
         	MCAccep[type][MCMULTI] += 1.0;
 
         	for (int id=0;id<NDIM;id++)                // save new coordinates
         	for (int it=time;it<=(time+seg_size);it++)    
         	{  
            	int pit = it % NumbTimes;                // periodicity in time 	       	
            	MCCoords[id][offset+pit] = newcoords[id][offset+pit];
         	}                                                           
      	}	     
//-----------------------------------------------------------------------  
//      END bisection 
//-----------------------------------------------------------------------
  	}  // END loop over time slices/atoms
}

void MCBisectionMoveExchange(int type, int time0)  // multilevel Metropolis
{
   int numb = MCAtom[type].numb;

   double mclambda = MCAtom[type].lambda;    
   int    mclevels = MCAtom[type].levels;  // number of levels
   int    seg_size = MCAtom[type].mlsegm;  // segment size  

   for (int atom=0;atom<numb;atom++)       // one atom to move only
   {
      int offset0 = MCAtom[type].offset + NumbTimes*atom;
      int offset1 = offset0;

      int time1 = time0 + seg_size;  // the end of the segment 
      int timep = time1 % NumbTimes;

      if (timep != time1)
      offset1 = MCAtom[type].offset + NumbTimes*PIndex[atom];

      for (int id=0;id<NDIM;id++)   
      {  
         newcoords[id][offset0 + time0] = MCCoords[id][offset0 + time0];
         newcoords[id][offset1 + timep] = MCCoords[id][offset1 + timep];
      }                                                      

      double bnorm = 1.0/(mclambda*MCTau);  // variance for gaussian sampling 

      bool Accepted; 

      int t0,t1,t2;

      double pot0 = 0.0;  // potential, current  level
      double pot1 = 0.0;  // potential, previous level
       	
      for (int level=0;level<mclevels;level++) // loop over bisection levels
      {	                                                          
         int level_seg_size = (int)pow(2.0,(mclevels-level));

         double bkin_norm = bnorm/(double) level_seg_size;
         double bpot_norm = MCTau*(double)(level_seg_size/2);
	   
         pot1 = pot0;    // swap level potentials
         pot0 = 0.0;
	   
         t2 = 0;
         do              // loop over middle points
         {
            t0 =  t2;                       // left point
            t2 =  t0 + level_seg_size;      // right point	
            t1 = (t0 + t2)/2;               // middle point

            int pt0 = (time0 + t0) % NumbTimes;	
            int pt1 = (time0 + t1) % NumbTimes;	
            int pt2 = (time0 + t2) % NumbTimes;	

            int off0 = offset0; 
            int off1 = offset0; 
            int off2 = offset0; 
 
            if (pt0 != (time0 + t0))
            off0 = offset1;
 
            if (pt1 != (time0 + t1))
            off1 = offset1;

            if (pt2 != (time0 + t2))
            off2 = offset1;

//  change the offset if exchange
 
            for (int id=0;id<NDIM;id++)
            {  	   
               newcoords[id][off1+pt1]  = 0.5*(newcoords[id][off0+pt0] + newcoords[id][off2+pt2]);
               newcoords[id][off1+pt1] += gauss(bkin_norm);
            } 
//---------------------------- the end point approximation
     
            int gatom0  = off0/NumbTimes;
            int gatom1  = off1/NumbTimes;

            pot0 += (PotEnergy(gatom0,newcoords,pt1) - PotEnergy(gatom0,MCCoords,pt1));
  
            if (t0!=0)                // skip the contributions of the end points
            pot0 += (PotEnergy(gatom0,newcoords,pt0) - PotEnergy(gatom0,MCCoords,pt0));
         }   	      
         while (t2<seg_size);        // end the loop over middle points 

// inefficient version

         double deltav = (pot0-2.0*pot1);  // rho(0,1;tau) 
         deltav *= bpot_norm;
 
         Accepted = false;
       
         if (deltav<0.0)               Accepted = true;
         else if (exp(-deltav)>rnd3()) Accepted = true;

         if (!Accepted) break;

     }  // END loop over levels        

     MCTotal[type][MCMULTI] += 1.0;
     
     if (Accepted)     
     {
         MCAccep[type][MCMULTI] += 1.0;
 
         for (int id=0;id<NDIM;id++)                // save new coordinates
         for (int it=time0;it<=time1;it++)    
         {  
            int pit = it % NumbTimes;               // periodicity in time
 	
            int offset = offset0;
            if (pit != it)
            offset = offset1;

            MCCoords[id][offset+pit] = newcoords[id][offset+pit];
         }                                                           
      }	     
//-----------------------------------------------------------------------  
//      END bisection 
//-----------------------------------------------------------------------
  }  // END loop over time slices/atoms
}

void MCRotationsMove(int type) // update all time slices for rotational degrees of freedom
{
#ifdef DEBUG_PIMC
	const char *_proc_=__func__;    //  MCRotationsMove() 
   	if (type != IMTYPE)
   	nrerror(_proc_,"Wrong impurity type");

   	if (NDIM != 3)
   	nrerror(_proc_,"Rotational sampling for 3D systems only");
#endif

   	double step   = MCAtom[type].rtstep; 
   	double MCRotChunkTot = 0.0;
   	double MCRotChunkAcp = 0.0;

   	RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"
   	double rand1,rand2,rand3;
	int offset, gatom;

#pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3,offset,gatom)
	for (int itrot=0;itrot<NumbRotTimes;itrot += 2)
	{
		for(int atom0=0;atom0<MCAtom[type].numb;atom0++)
		{
			offset = MCAtom[type].offset+(NumbTimes*atom0);   // the same offset for rotational
			gatom  = offset/NumbTimes;    // and translational degrees of freedom
			rand1=runif(Rng);
			rand2=runif(Rng);
			rand3=runif(Rng);
#ifdef PIMCTYPE
			MCRotLinStepPIMC(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
//
#ifdef PIGSTYPE
			MCRotLinStepPIGS(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
//
#ifdef PIGSENTTYPE
#ifdef SWAPTOUNSWAP
			MCRotLinStepSwap(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp,Distribution);
#else
			MCRotLinStepSwapBroken(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
#endif
		}
	}

	MCTotal[type][MCROTAT] += MCRotChunkTot;
	MCAccep[type][MCROTAT] += MCRotChunkAcp;

	MCRotChunkTot = 0;
	MCRotChunkAcp = 0;

#pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3,offset,gatom)
	for (int itrot = 1; itrot < NumbRotTimes; itrot += 2)
	{
		for(int atom0=0;atom0<MCAtom[type].numb;atom0++)
		{
			offset = MCAtom[type].offset+(NumbTimes*atom0);   // the same offset for rotational
			gatom  = offset/NumbTimes;    // and translational degrees of freedom
 			rand1=runif(Rng);
			rand2=runif(Rng);
			rand3=runif(Rng);
#ifdef PIMCTYPE
			MCRotLinStepPIMC(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
//
#ifdef PIGSTYPE
			MCRotLinStepPIGS(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
//
#ifdef PIGSENTTYPE
#ifdef SWAPTOUNSWAP
			MCRotLinStepSwap(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp,Distribution);
#else
			MCRotLinStepSwapBroken(itrot,offset,gatom,type,step,rand1,rand2,rand3,MCRotChunkTot,MCRotChunkAcp);
#endif
#endif
		}
	}

	MCTotal[type][MCROTAT] += MCRotChunkTot;
	MCAccep[type][MCROTAT] += MCRotChunkAcp;

#ifdef SWAPTOUNSWAP
    //double rand4 = (double)rand() / ((double)RAND_MAX + 1);
    double rand4 = runif(Rng);
    MCSwap(rand4, Distribution);
    if (Distribution == "Swap") 
	{
		MCAccepSwap += 1;
		iSwap = 1;
		iUnSwap = 0;
	}
    else
	{
		MCAccepUnSwap += 1;
		iSwap = 0;
		iUnSwap = 1;
	}
#endif
}

// Perform MC Rigid Body Rotations

void MCRotLinStepPIMC(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp)
{
	int it0 = (it1 - 1);
	int it2 = (it1 + 1);

   	if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   	if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   	int t0 = offset + it0;
   	int t1 = offset + it1;
   	int t2 = offset + it2;

   	double cost = MCAngles[CTH][t1];
   	double phi  = MCAngles[PHI][t1];
	double EulangOld[NDIM], EulangNew[NDIM];
	EulangOld[PHI] = phi;
	EulangOld[CTH] = acos(cost);
	EulangOld[CHI] = 0.0;

   	cost += (step*(rand1-0.5));
   	phi  += (step*(rand2-0.5));

   	if (cost >  1.0)
   	{
      	cost = 2.0 - cost;
   	}

   	if (cost < -1.0)
   	{
       	cost = -2.0 - cost;
   	}

	if (abs(cost) > 2.0) 
	{
        cout<<"Upper or lower limit of cost is excided " << cost<<endl;
		exit(0);
	}

	EulangNew[PHI] = phi;
	EulangNew[CTH] = acos(cost);
	EulangNew[CHI] = 0.0;

   	double sint = sqrt(1.0 - cost*cost);

   	newcoords[AXIS_X][t1] = sint*cos(phi);
   	newcoords[AXIS_Y][t1] = sint*sin(phi);
   	newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// 	the old density

   	double p0 = 0.0;
   	double p1 = 0.0;

   	for (int id=0;id<NDIM;id++)
   	{
      	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      	p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_old;
   	double rho1,rho2,erot;
// 	If it1 = 0 (the first bead), dens_new = SRotDens(p1,type)
// 	if it1 = (NumbRotTimes-1) is the last bead, dens_new = SRotDens(p0,type)

	if(RotDenType == 0)
	{
        dens_old = SRotDens(p0,type)*SRotDens(p1,type);
	}
    else if(RotDenType == 1)
    {
        rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
        rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
        dens_old = rho1+rho2; //dens_old
    }

   if (fabs(dens_old)<RZERO) dens_old = 0.0;
#ifndef NEGATIVEDENSITY
   if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
   if (dens_old<0.0) dens_old=fabs(dens_old);
#endif

   double pot_old  = 0.0;

   int itr0 = it1  * RotRatio;     // interval to average over
   int itr1 = itr0 + RotRatio;     // translational time slices

   	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
   		//pot_old  += (PotRotEnergyPIMC(gatom,MCCosine,it));
   		pot_old  += (PotRotEnergyPIMC(gatom,EulangOld,it));
	}

// the new density 

   p0 = 0.0;
   p1 = 0.0;


   for (int id=0;id<NDIM;id++)
   {
       p0 += (MCCosine[id][t0]*newcoords[id][t1]);
       p1 += (newcoords[id][t1]*MCCosine[id][t2]);
   }

   double dens_new;

	if(RotDenType == 0)
	{
        dens_new = SRotDens(p0,type)*SRotDens(p1,type);
	}
	else if(RotDenType == 1)
	{
		rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
		rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
		dens_new = rho1 + rho2;
	}

	if (fabs(dens_new)<RZERO) dens_new = 0.0;
#ifndef NEGATIVEDENSITY
	if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
	if (dens_new<0.0) dens_new=fabs(dens_new);
#endif

	double pot_new  = 0.0;

	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		//pot_new  += (PotRotEnergyPIMC(gatom,newcoords,it));
		pot_new  += (PotRotEnergyPIMC(gatom,EulangNew,it));
	}

	double rd;

	if(RotDenType == 0)
	{
		if (dens_old>RZERO)
			rd = dens_new/dens_old;
		else rd = 1.0;

		rd *= exp(- MCTau*(pot_new-pot_old));
	}
	else if(RotDenType == 1)
	{
		rd = dens_new - dens_old - MCTau*(pot_new-pot_old);
		//rd = exp(rd);
	}

	bool Accepted = false;
	if(RotDenType == 0)
	{
		if (rd>1.0)         Accepted = true;
		//else if (rd>rnd7()) Accepted = true;
		else if (rd>rand3) Accepted = true;
	}
	else if (RotDenType == 1)
	{
		if (rd > 0.0)   Accepted = true;
		//else if (rd > log(rnd7())) Accepted = true;
    	else if (rd > log(rand3)) Accepted = true;
	}

	MCRotChunkTot += 1.0;

	if (Accepted)
	{
		MCRotChunkAcp += 1.0;

		MCAngles[CTH][t1] = cost;
		MCAngles[PHI][t1] = phi;

		for (int id=0;id<NDIM;id++)
    		MCCosine[id][t1] = newcoords[id][t1];
	}

}

void MCRotLinStepPIGS(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp)
{
	int it0 = (it1 - 1);
	int it2 = (it1 + 1);

   	if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   	if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   	int t0 = offset + it0;
   	int t1 = offset + it1;
   	int t2 = offset + it2;

   	double cost = MCAngles[CTH][t1];
   	double phi  = MCAngles[PHI][t1];
	double EulangOld[NDIM], EulangNew[NDIM];
	EulangOld[PHI] = phi;
	EulangOld[CTH] = acos(MCAngles[CTH][t1]);
	EulangOld[CHI] = 0.0;

   	cost += (step*(rand1-0.5));
   	phi  += (step*(rand2-0.5));

   	if (cost >  1.0)
   	{
      	cost = 2.0 - cost;
   	}

   	if (cost < -1.0)
   	{
       	cost = -2.0 - cost;
   	}

	if (abs(cost) > 2.0) 
	{
        cout<<"Upper or lower limit of cost is excided " << cost<<endl;
		exit(0);
	}

	EulangNew[PHI] = phi;
	EulangNew[CTH] = acos(cost);
	EulangNew[CHI] = 0.0;

   	double sint = sqrt(1.0 - cost*cost);

   	newcoords[AXIS_X][t1] = sint*cos(phi);
   	newcoords[AXIS_Y][t1] = sint*sin(phi);
   	newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// 	the old density

   	double p0 = 0.0;
   	double p1 = 0.0;

   	for (int id=0;id<NDIM;id++)
   	{
      	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      	p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_old;
   	double rho1,rho2,erot;
// 	If it1 = 0 (the first bead), dens_new = SRotDens(p1,type)
// 	if it1 = (NumbRotTimes-1) is the last bead, dens_new = SRotDens(p0,type)

	if(RotDenType == 0)
	{
        if (it1 == 0 || it1 == (NumbRotTimes - 1))
        {
            if (it1 == 0)
            {
                dens_old = SRotDens(p1, type);
            }
            else
            {
                dens_old = SRotDens(p0, type);
            }
        }
        else
        {
            dens_old = SRotDens(p0,type)*SRotDens(p1,type);
        }
	}
    else if(RotDenType == 1)
    {
        rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
        rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
        dens_old = rho1+rho2;
    }

   if (fabs(dens_old)<RZERO) dens_old = 0.0;
#ifndef NEGATIVEDENSITY
   if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
   if (dens_old<0.0) dens_old=fabs(dens_old);
#endif

   double pot_old  = 0.0;

   int itr0 = it1  * RotRatio;     // interval to average over
   int itr1 = itr0 + RotRatio;     // translational time slices

   	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
   		//pot_old  += (PotRotEnergyPIGS(gatom,MCCosine,it));
   		pot_old  += (PotRotEnergyPIGS(gatom,EulangOld,it));
	}

// the new density 

   p0 = 0.0;
   p1 = 0.0;


   for (int id=0;id<NDIM;id++)
   {
       p0 += (MCCosine[id][t0]*newcoords[id][t1]);
       p1 += (newcoords[id][t1]*MCCosine[id][t2]);
   }

   double dens_new;

	if(RotDenType == 0)
	{
        if ((it1 == 0) || (it1 == (NumbRotTimes - 1)))
        {
            if (it1 == 0)
            {
                dens_new = SRotDens(p1, type);
            }
            else
            {
                dens_new = SRotDens(p0, type);
            }
        }
        else
        {
            dens_new = SRotDens(p0,type)*SRotDens(p1,type);
        }
	}
	else if(RotDenType == 1)
	{
		rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
		rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
		dens_new = rho1 + rho2;
	}

	if (fabs(dens_new)<RZERO) dens_new = 0.0;
#ifndef NEGATIVEDENSITY
	if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
	if (dens_new<0.0) dens_new=fabs(dens_new);
#endif

	double pot_new  = 0.0;

	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		//pot_new  += (PotRotEnergyPIGS(gatom,newcoords,it));
		pot_new  += (PotRotEnergyPIGS(gatom,EulangNew,it));
	}

	double rd;

	if(RotDenType == 0)
	{
		if (dens_old>RZERO)
			rd = dens_new/dens_old;
		else rd = 1.0;

		rd *= exp(- MCTau*(pot_new-pot_old));
	}
	else if(RotDenType == 1)
	{
		rd = dens_new - dens_old - MCTau*(pot_new-pot_old);
		//rd = exp(rd);
	}

	bool Accepted = false;
	if(RotDenType == 0)
	{
		if (rd>1.0)         Accepted = true;
		//else if (rd>rnd7()) Accepted = true;
		else if (rd>rand3) Accepted = true;
	}
	else if (RotDenType == 1)
	{
		if (rd > 0.0)   Accepted = true;
		//else if (rd > log(rnd7())) Accepted = true;
    	else if (rd > log(rand3)) Accepted = true;
	}

	MCRotChunkTot += 1.0;

	if (Accepted)
	{
		MCRotChunkAcp += 1.0;

		MCAngles[CTH][t1] = cost;
		MCAngles[PHI][t1] = phi;

		for (int id=0;id<NDIM;id++)
    		MCCosine[id][t1] = newcoords[id][t1];
	}

}

double PotRotEnergyPIGS(int atom0, double *Eulang0, int it)   
{
	int type0   =  MCType[atom0];
	double spot;

    double weight;
	weight = 1.0;
    if (it == 0 || it == (NumbRotTimes - 1)) weight = 0.5;

	if ( (MCAtom[type0].molecule == 4) && (MCAtom[type0].numb > 1) )
	{
	    int offset0 =  atom0*NumbRotTimes;
        int t0  = offset0 + it;

        spot = 0.0;
        for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
        if (atom1 != atom0)                    
        {
            int offset1 = atom1*NumbRotTimes;
            int t1  = offset1 + it;

	        string stype = MCAtom[type0].type;
			/*
			if (stype == H2)
	        {
				double cosine[NDIM][NumbAtoms*NumbRotTimes];
				cosine[0][t0] = sin(Eulang0[CTH])*cos(Eulang0[PHI]);
				cosine[1][t0] = sin(Eulang0[CTH])*sin(Eulang0[PHI]);
				cosine[2][t0] = cos(Eulang0[CTH]);
                double s1 = 0.0;
                double s2 = 0.0;
                double dr2 = 0.0;
				double dr[NDIM];

                for (int id=0;id<NDIM;id++)
                {
                    dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);
                    dr2    += (dr[id]*dr[id]);
                    double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*cosine[id][t0];
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
                    b1[id] = cosine[id][t0];
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
                spot += weight*potl*CMRECIP2KL;
			}  //stype
			*/

		    if (stype == HF )
            {
				double Eulang1[NDIM];
				Eulang1[PHI] = MCAngles[PHI][t1];
        		Eulang1[CTH] = acos(MCAngles[CTH][t1]);
        		Eulang1[CHI] = 0.0;
        		spot += weight*PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            }  //stype
        } //loop over atom1 (molecules)
    }

	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
	{
        double E12 = -2.0*DipoleMomentAU2*cos(Eulang0[CTH])/(RR*RR*RR);
        spot        = weight*E12*AuToKelvin;
    }
    double spot_cage;
#ifdef CAGEPOT
    double cost = cos(Eulang0[CTH]);
    double phi = Eulang0[PHI];
    if (phi < 0.0) phi = 2.0*M_PI + phi;
    phi = fmod(phi,2.0*M_PI);
    spot_cage = weight*LPot2DRotDOF(cost,phi,type0);
#else
    spot_cage = 0.0;
#endif
	double spotReturn = (spot + spot_cage);
    return spotReturn;
}

void MCRotLinStepSwap(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp, string Distribution)
{
	int it0 = (it1 - 1);
   	int it2 = (it1 + 1);

   	if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   	if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   	int t0 = offset + it0;
   	int t1 = offset + it1;
   	int t2 = offset + it2;

   	double cost = MCAngles[CTH][t1];
   	double phi  = MCAngles[PHI][t1];

	double EulangOld[NDIM], EulangNew[NDIM];
	EulangOld[PHI] = phi;
	EulangOld[CTH] = acos(cost);
	EulangOld[CHI] = 0.0;

   	cost += (step*(rand1-0.5));
   	phi  += (step*(rand2-0.5));

   	if (cost >  1.0)
   	{
      	cost = 2.0 - cost;
   	}

   	if (cost < -1.0)
   	{
       	cost = -2.0 - cost;
   	}

	if (abs(cost) > 2.0) 
	{
        cout<<"Upper or lower limit of cost is excided " << cost<<endl;
		exit(0);
	}

	EulangNew[PHI] = phi;
	EulangNew[CTH] = acos(cost);
	EulangNew[CHI] = 0.0;

   	double sint = sqrt(1.0 - cost*cost);

   	newcoords[AXIS_X][t1] = sint*cos(phi);
   	newcoords[AXIS_Y][t1] = sint*sin(phi);
   	newcoords[AXIS_Z][t1] = cost;

// the old density

   	double p0 = 0.0;
   	double p1 = 0.0;

   	for (int id=0;id<NDIM;id++)
   	{
      	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      	p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_old;
// If it1 = 0 (the first bead), dens_new = SRotDens(p1,type)
// if it1 = (NumbRotTimes-1) is the last bead, dens_new = SRotDens(p0,type)
   	int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	int particleA1Max = particleA1Min + NumbParticle - 1;
   	int particleA2Min = particleA1Max + 1;
   	int particleA2Max = particleA2Min + NumbParticle - 1;

	if (it1 == 0 || it1 == (NumbRotTimes - 1))
	{
		if (it1 == 0)
		{	
			dens_old = SRotDens(p1, type);
		}
		else
		{
			dens_old = SRotDens(p0, type);
		}
	}
	else
	{
		dens_old = SRotDens(p0,type)*SRotDens(p1,type);
	}
	if (Distribution == "Swap")
	{
       	if ((gatom >= particleA1Min) && (gatom <= particleA1Max))
       	{
       		int gatomSwap = particleA2Max - (gatom - particleA1Min);
           	int offsetSwap = NumbRotTimes*gatomSwap;

         	if (it1 == (((NumbRotTimes - 1)/2) - 1))
           	{
           		int tSwap = offsetSwap + it2;
				double pSwap = 0.0;
           		for (int id = 0; id < NDIM; id++)
           		{
               		pSwap += (MCCosine[id][t1]*MCCosine[id][tSwap]);
           		}
               	dens_old = SRotDens(p0,type)*SRotDens(pSwap,type);
           	}
           	if (it1 == ((NumbRotTimes - 1)/2))
           	{
           		int tSwap = offsetSwap + it0;
				double pSwap = 0.0;
           		for (int id = 0; id < NDIM; id++)
           		{
               		pSwap += (MCCosine[id][tSwap]*MCCosine[id][t1]);
           		}
               	dens_old = SRotDens(pSwap,type)*SRotDens(p1,type);
           	}
       	}	

       	if ((gatom >= particleA2Min) && (gatom <= particleA2Max))
       	{
           	int gatomSwap = particleA1Max - (gatom - particleA2Min);
           	int offsetSwap = NumbRotTimes*gatomSwap;

           	if (it1 == (((NumbRotTimes - 1)/2) - 1))
           	{
           		int tSwap = offsetSwap + it2;
				double pSwap = 0.0;
           		for (int id = 0; id < NDIM; id++)
           		{
               		pSwap += (MCCosine[id][t1]*MCCosine[id][tSwap]);
           		}
               	dens_old = SRotDens(p0,type)*SRotDens(pSwap,type);
           	}
           	if (it1 == ((NumbRotTimes - 1)/2))
           	{
           		int tSwap = offsetSwap + it0;
				double pSwap = 0.0;
           		for (int id = 0; id < NDIM; id++)
           		{
               		pSwap += (MCCosine[id][tSwap]*MCCosine[id][t1]);
           		}
               	dens_old = SRotDens(pSwap,type)*SRotDens(p1,type);
           	}
       	}
	}
//
   	if (fabs(dens_old)<RZERO) dens_old = 0.0;
#ifndef NEGATIVEDENSITY
   	if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
   	if (dens_old<0.0) dens_old=fabs(dens_old);
#endif
//
   	double pot_old  = 0.0;

   	int itr0 = it1  * RotRatio;     // interval to average over
   	int itr1 = itr0 + RotRatio;     // translational time slices

   	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		pot_old  += (PotRotEnergySwap(gatom,EulangOld,it,Distribution));
	}

// the new density 

   	p0 = 0.0;
   	p1 = 0.0;


   	for (int id=0;id<NDIM;id++)
   	{
       	p0 += (MCCosine[id][t0]*newcoords[id][t1]);
       	p1 += (newcoords[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_new;

	if ((it1 == 0) || (it1 == (NumbRotTimes - 1)))
	{
		if (it1 == 0)
		{
			dens_new = SRotDens(p1, type);
		}
		else
		{
			dens_new = SRotDens(p0, type);
		}
	}
	else
	{
		dens_new = SRotDens(p0,type)*SRotDens(p1,type);
	}
	if (Distribution == "Swap")
    {
       	if ((gatom >= particleA1Min) && (gatom <= particleA1Max))
       	{
           	int gatomSwap = particleA2Max - (gatom - particleA1Min);
           	int offsetSwap = NumbRotTimes*gatomSwap;

           	if (it1 == (((NumbRotTimes - 1)/2) - 1))
           	{
               	int tSwap = offsetSwap + it2;
               	double pSwap = 0.0;
               	for (int id = 0; id < NDIM; id++)
               	{
                   	pSwap += (newcoords[id][t1]*MCCosine[id][tSwap]);
               	}
               	dens_new = SRotDens(p0,type)*SRotDens(pSwap,type);
           	}
           	if (it1 == ((NumbRotTimes - 1)/2))
           	{
               	int tSwap = offsetSwap + it0;
               	double pSwap = 0.0;
               	for (int id = 0; id < NDIM; id++)
               	{
                   	pSwap += (MCCosine[id][tSwap]*newcoords[id][t1]);
               	}
               	dens_new = SRotDens(pSwap,type)*SRotDens(p1,type);
           	}
       	}

    	if ((gatom >= particleA2Min) && (gatom <= particleA2Max))
   	    {
           	int gatomSwap = particleA1Max - (gatom - particleA2Min);
           	int offsetSwap = NumbRotTimes*gatomSwap;

           	if (it1 == (((NumbRotTimes - 1)/2) - 1))
           	{
               	int tSwap = offsetSwap + it2;
               	double pSwap = 0.0;
               	for (int id = 0; id < NDIM; id++)
               	{
                   	pSwap += (newcoords[id][t1]*MCCosine[id][tSwap]);
               	}
               	dens_new = SRotDens(p0,type)*SRotDens(pSwap,type);
           	}
           	if (it1 == ((NumbRotTimes - 1)/2))
           	{
               	int tSwap = offsetSwap + it0;
               	double pSwap = 0.0;
               	for (int id = 0; id < NDIM; id++)
               	{
                   	pSwap += (MCCosine[id][tSwap]*newcoords[id][t1]);
               	}
               	dens_new = SRotDens(pSwap,type)*SRotDens(p1,type);
           	}
       	}
	}

//
	if (fabs(dens_new)<RZERO) dens_new = 0.0;
#ifndef NEGATIVEDENSITY
	if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
	if (dens_new<0.0) dens_new=fabs(dens_new);
#endif
//
	double pot_new  = 0.0;
	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		pot_new  += (PotRotEnergySwap(gatom,EulangNew,it,Distribution));
	}

	double rd;
	if (dens_old>RZERO)
		rd = dens_new/dens_old;
	else rd = 1.0;

	rd *= exp(- MCTau*(pot_new-pot_old));

	bool Accepted = false;
	if (rd>1.0)         Accepted = true;
	else if (rd>rand3) Accepted = true;

	MCRotChunkTot += 1.0;

	if (Accepted)
	{
		MCRotChunkAcp += 1.0;

		MCAngles[CTH][t1] = cost;
		MCAngles[PHI][t1] = phi;

		for (int id=0;id<NDIM;id++)
    		MCCosine[id][t1] = newcoords[id][t1];
	}
}

//double PotRotEnergySwap(int atom0, double **cosine, int it, int Distribution)   
double PotRotEnergySwap(int atom0, const double *Eulang0, int it, string Distribution)   
{
	int type0   =  MCType[atom0];
	double spot, spotSwap;

    double weight, weight1;
	weight = 1.0;

    if (it == 0 || it == (NumbRotTimes - 1)) weight = 0.5;

	if ( (MCAtom[type0].molecule == 4) && (MCAtom[type0].numb > 1) )
	{
	    int offset0 =  atom0*NumbRotTimes;
        int t0  = offset0 + it;

		int atom1Init, NumbAtoms1;
   	    int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	    int particleA1Max = particleA1Min + NumbParticle - 1;
   	    int particleA2Min = particleA1Max + 1;
   	    int particleA2Max = particleA2Min + NumbParticle - 1;

        if (atom0 <= particleA1Max)
        {
            atom1Init  = 0;
            NumbAtoms1 = (particleA1Max+1);
        }
        else
        {
            atom1Init  = particleA2Min;
            NumbAtoms1 = NumbAtoms;
        }

        spot = 0.0;
        for (int atom1 = atom1Init; atom1 < NumbAtoms1; atom1++)
        if (atom1 != atom0)                    
        {
            int offset1 = atom1*NumbRotTimes;
            int t1  = offset1 + it;

	        string stype = MCAtom[type0].type;
		    if (stype == HF )
            {
				double Eulang1[NDIM];
				Eulang1[PHI] = MCAngles[PHI][t1];
                Eulang1[CTH] = acos(MCAngles[CTH][t1]);
                Eulang1[CHI] = 0.0;

				weight1 = 1.0;
				if (Distribution == "Swap")
				{
					if (((atom0 < particleA1Min) || (atom0 > particleA2Max)) && ((atom1 >= particleA1Min) && (atom1 <= particleA2Max)))
            		{
           				if (it == ((NumbRotTimes - 1)/2))
              			{
                   			weight1 = 0.5;
						}
            		} 
            		if (((atom0 >= particleA1Min) && (atom0 <= particleA2Max)) && ((atom1 < particleA1Min) || (atom1 > particleA2Max)))
            		{
            			if (it == ((NumbRotTimes - 1)/2))
            	   		{
        	           		weight1 = 0.5;
						}
	            	} 
				}
				spot += weight*weight1*PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            }  //stype
        } //loop over atom1 (molecules)
		spotSwap = 0.0;
		if (Distribution == "Swap")
		{
			if (it == ((NumbRotTimes - 1)/2))
			{
				int atomSwapMin, atomSwapMax;
				if (atom0 < particleA1Min)
    	   		{
					atomSwapMin = particleA2Min;
                	atomSwapMax = (particleA2Max + 1);
				}
				if (atom0 > particleA2Max)
    	   		{
					atomSwapMin = particleA1Min;
                	atomSwapMax = (particleA1Max + 1);
				}
	    		if ((atom0 >= particleA1Min) && (atom0 <= particleA1Max))
       	    	{
			    	atomSwapMin = (particleA2Max + 1);
                	atomSwapMax = NumbAtoms;
		    	}
     			if ((atom0 >= particleA2Min) && (atom0 <= particleA2Max))
         		{
	     			atomSwapMin = 0;
                	atomSwapMax = particleA1Min;
		    	}
				spotSwap = 0.0;
     			for (int atomSwap = atomSwapMin; atomSwap < atomSwapMax; atomSwap++)
	    		{
                	int offsetSwap = NumbRotTimes*atomSwap;
                	int tSwap  = offsetSwap + it;
	        		string stype = MCAtom[type0].type;
                	if (stype == HF )
                	{
						double EulangSwap[NDIM];
						EulangSwap[PHI] = MCAngles[PHI][tSwap];
                		EulangSwap[CTH] = acos(MCAngles[CTH][tSwap]);
                		EulangSwap[CHI] = 0.0;
                    	spotSwap += 0.5*PotFunc(atom0, atomSwap, Eulang0, EulangSwap, it);
                	}  //stype
		    	}
        	}
		}
    }

    double spot_cage;
#ifdef CAGEPOT
	double cost = cos(Eulang0[CTH]);
	double phi = Eulang0[PHI];
    if (phi < 0.0) phi = 2.0*M_PI + phi;
    phi = fmod(phi,2.0*M_PI);
    spot_cage = weight*LPot2DRotDOF(cost,phi,type0);
#else
    spot_cage = 0.0;
#endif
	double spotReturn = spot + spotSwap + spot_cage;
    return spotReturn;
}

void MCSwap(double rand4, string &Distribution)
{

    double rd;

    if (Distribution == "unSwap")
    {
        rd = GetEstimNM()/GetEstimDM();
    }

    if (Distribution == "Swap")
    {
        rd = GetEstimDM()/GetEstimNM();
    }

    bool Accepted = false;
    if (rd>1.0)         Accepted = true;
    else if (rd>rand4) Accepted = true;

    string DistributionInit = Distribution;
    if (Accepted)
    {
        if (DistributionInit == "unSwap") Distribution = "Swap";
        if (DistributionInit == "Swap" ) Distribution = "unSwap";
    }
}

/*
void MCRotationsMove(int type) // update all time slices for rotational degrees of freedom
{
#ifdef DEBUG_PIMC
   const char *_proc_=__func__;    //  MCRotationsMove() 
   if (type != IMTYPE)
   nrerror(_proc_,"Wrong impurity type");

   if (NDIM != 3)
   nrerror(_proc_,"Rotational sampling for 3D systems only");
#endif

   double step   = MCAtom[type].rtstep; 
   int    offset = MCAtom[type].offset;
 
   int atom0  = 0;                   // only one molecular impurtiy
   offset    += (NumbTimes*atom0);   // the same offset for rotational
   int gatom  = offset/NumbTimes;    // and translational degrees of freedom

   for (int it1=0;it1<NumbRotTimes;it1++)
   {
      int it0 = (it1 - 1);
      int it2 = (it1 + 1);
 
      if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
      if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0
      
      int t0 = offset + it0;
      int t1 = offset + it1;
      int t2 = offset + it2;

      double n1[NDIM];

      double cost = MCAngles[CTH][t1];
      double phi  = MCAngles[PHI][t1];

      cost += (step*(rnd1()-0.5));
      phi  += (step*(rnd1()-0.5));

      if (cost >  1.0)
      {
         cost = 2.0 - cost;    
//       phi  = phi + M_PI;
      }		 
      
      if (cost < -1.0)
      {
          cost = -2.0 - cost;    
//        phi  = phi  + M_PI;
      }		  

      double sint = sqrt(1.0 - cost*cost);

      newcoords[AXIS_X][t1] = sint*cos(phi);
      newcoords[AXIS_Y][t1] = sint*sin(phi);
      newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// the old density

      double p0 = 0.0;
      double p1 = 0.0;
 
      for (int id=0;id<NDIM;id++)
      {
         p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
         p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
      }

      double dens_old;
      double rho1,rho2,erot;

      if(RotDenType == 0)
      {
         dens_old = SRotDens(p0,type)*SRotDens(p1,type);
      }
      else if(RotDenType == 1)
      {
         rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
         rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
         dens_old = rho1+rho2;
      }

      if (fabs(dens_old)<RZERO) dens_old = 0.0;
      if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");

      double pot_old  = 0.0;

      int itr0 = it1  * RotRatio;     // interval to average over
      int itr1 = itr0 + RotRatio;     // translational time slices
 
      for (int it=itr0;it<itr1;it++)  // average over tr time slices
      pot_old  += (PotRotEnergy(gatom,MCCosine,it));

//   the new density 
    
      p0 = 0.0;
      p1 = 0.0;
 

      for (int id=0;id<NDIM;id++)
      {
          p0 += (MCCosine [id][t0]*newcoords[id][t1]);
          p1 += (newcoords[id][t1]*MCCosine [id][t2]);
      }

      double dens_new;

      if(RotDenType == 0)
      {
         dens_new = SRotDens(p0,type)*SRotDens(p1,type);
      }
      else if(RotDenType == 1)
      {
         rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
         rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
         dens_new = rho1 + rho2;
      }

      if (fabs(dens_new)<RZERO) dens_new = 0.0;
      if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");

      double pot_new  = 0.0;

      for (int it=itr0;it<itr1;it++)  // average over tr time slices
      pot_new  += (PotRotEnergy(gatom,newcoords,it));
     
      double rd;

      if(RotDenType == 0)
      {
         if (dens_old>RZERO)
          rd = dens_new/dens_old;
         else rd = 1.0;

         rd *= exp(- MCTau*(pot_new-pot_old));   
      }
      else if(RotDenType == 1)
      {
         rd = dens_new - dens_old - MCTau*(pot_new-pot_old);
//       rd = exp(rd);
      }
     

      bool Accepted = false;
      if(RotDenType == 0)
      {
         if (rd>1.0)         Accepted = true;
         else if (rd>rnd7()) Accepted = true;
      }
      else if (RotDenType == 1)
      {
         if (rd > 0.0)   Accepted = true;
         else if (rd > log(rnd7())) Accepted = true;
      }

      MCTotal[type][MCROTAT] += 1.0;  
      
      if (Accepted)
      {
         MCAccep[type][MCROTAT] += 1.0;

         MCAngles[CTH][t1] = cost;
         MCAngles[PHI][t1] = phi;
  
         for (int id=0;id<NDIM;id++)
         MCCosine [id][t1] = newcoords[id][t1];
      }	      
   } // end of the loop over time slices
}
*/

void MCRotations3D(int type) // update all time slices for rotational degrees of freedom
{
#ifdef DEBUG_PIMC
   const char *_proc_=__func__;    //  MCRotationsMove() 
   if (type != IMTYPE)
   nrerror(_proc_,"Wrong impurity type");

   if (NDIM != 3)
   nrerror(_proc_,"Rotational sampling for 3D systems only");
#endif


   double step   = MCAtom[type].rtstep; 
 
// for(int atom0=0;atom0<MCAtom[type].numb;atom0++)
// {
//    int offset = MCAtom[type].offset+(NumbTimes*atom0);   // the same offset for rotational
//    int gatom  = offset/NumbTimes;    // and translational degrees of freedom


//    serial code
/*
      MCRotChunkTot = 0;
      MCRotChunkAcp = 0;

      for (int it1=0;it1<NumbRotTimes;it1++)
      {
        MCRot3Dstep(it1,offset,gatom,type,step,MCRotChunkTot,MCRotChunkAcp);
      }



   

      MCTotal[type][MCROTAT] += MCRotChunkTot;
      MCAccep[type][MCROTAT] += MCRotChunkAcp;
*/

// openmp code
   MCRotChunkTot = 0;
   MCRotChunkAcp = 0;
// randomseed();   //set seed according to clock
   RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"
   double rand1,rand2,rand3,rand4;

   #pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3,rand4)
   for (int itrot=0;itrot<NumbRotTimes;itrot=itrot+2)
   {
      for(int atom0=0;atom0<MCAtom[type].numb;atom0++)
      {
         int offset = MCAtom[type].offset+(NumbTimes*atom0);   // the same offset for rotational
         int gatom  = offset/NumbTimes;    // and translational degrees of freedom
         rand1=runif(Rng);
         rand2=runif(Rng);
         rand3=runif(Rng);
         rand4=runif(Rng);
         MCRot3Dstep(itrot,offset,gatom,type,step,rand1,rand2,rand3,rand4,IROTSYM,NFOLD_ROT,MCRotChunkTot,MCRotChunkAcp);
      }
   }

   MCTotal[type][MCROTAT] += MCRotChunkTot;
   MCAccep[type][MCROTAT] += MCRotChunkAcp;

   MCRotChunkTot = 0;
   MCRotChunkAcp = 0;

   #pragma omp parallel for reduction(+: MCRotChunkTot,MCRotChunkAcp) private(rand1,rand2,rand3,rand4)
   for (int itrot=1;itrot<NumbRotTimes;itrot=itrot+2)
   {
      for(int atom0=0;atom0<MCAtom[type].numb;atom0++)
      {
         int offset = MCAtom[type].offset+(NumbTimes*atom0);   // the same offset for rotational
         int gatom  = offset/NumbTimes;    // and translational degrees of freedom
         rand1=runif(Rng);
         rand2=runif(Rng);
         rand3=runif(Rng);
         rand4=runif(Rng);
         MCRot3Dstep(itrot,offset,gatom,type,step,rand1,rand2,rand3,rand4,IROTSYM,NFOLD_ROT,MCRotChunkTot,MCRotChunkAcp);
      }
   }

   MCTotal[type][MCROTAT] += MCRotChunkTot;
   MCAccep[type][MCROTAT] += MCRotChunkAcp;

// }

}

#ifdef IOWRITE
void MCRotLinStep(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp)
{
	int it0 = (it1 - 1);
	int it2 = (it1 + 1);

   	if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   	if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   	int t0 = offset + it0;
   	int t1 = offset + it1;
   	int t2 = offset + it2;

#ifndef PROPOSED
   	double cost = MCAngles[CTH][t1];
   	double phi  = MCAngles[PHI][t1];
	double EulangOld[NDIM], EulangNew[NDIM];
	EulangOld[PHI] = phi;
	EulangOld[CTH] = acos(cost);
	EulangOld[CHI] = 0.0;
#endif

#ifdef PROPOSED
	double PreDistribution[NCOST*NPHI];
    double sum = 0.0;
	for (int itp = 0; itp < NCOST*NPHI; itp++)
    {
		double sintProposed = sqrt(1.0 - costProposed[itp]*costProposed[itp]);
        tempcoords[0][t1] = sintProposed*cos(phiProposed[itp]);
        tempcoords[1][t1] = sintProposed*sin(phiProposed[itp]);
        tempcoords[2][t1] = costProposed[itp];

		double p0 = 0.0;
   		double p1 = 0.0;

   		for (int id=0;id<NDIM;id++)
   		{
   			p0 += (MCCosine[id][t0]*tempcoords[id][t1]);
   			p1 += (tempcoords[id][t1]*MCCosine[id][t2]);
   		}

		double weight;
      	if (it1 == 0 || it1 == (NumbRotTimes - 1))
        {
			weight = 0.5;
            if (it1 == 0)
            {
                PreDistribution[itp] = SRotDens(p1, type);
            }
            else
            {
                PreDistribution[itp] = SRotDens(p0, type);
            }
        }
        else
        {
			weight = 1.0;
            PreDistribution[itp] = SRotDens(p0,type)*SRotDens(p1,type);
        }
        
    	double E12 = -2.0*DipoleMomentAU2*costProposed[itp]*AuToKelvin*weight/(RR*RR*RR);
        PreDistribution[itp] *= exp(-MCRotTau*E12);
		//var+=MCRotTau*E12;
        sum += PreDistribution[itp];
	}
    //PreDistribution[itp]=exp(-var);
	for (int itp = 0; itp < NCOST*NPHI; itp++)
    {
        PreDistribution[itp] /= sum;
	}
   	iChooseNew = myRand(PreDistribution, rand1);
	if (iChooseNew == -1) exit(0);
   	double cost = costProposed[iChooseNew];
   	double phi  = phiProposed[iChooseNew];
	MCAngles[CTH][t1] = cost;
	MCAngles[PHI][t1] = phi;

	double sint = sqrt(1.0 - cost*cost);
	MCCosine[0][t1] = sint*cos(phi);
   	MCCosine[1][t1] = sint*sin(phi);
   	MCCosine[2][t1] = cost;
#endif
#ifndef PROPOSED
   	cost += (step*(rand1-0.5));
   	phi  += (step*(rand2-0.5));

   	if (cost >  1.0)
   	{
      	cost = 2.0 - cost;
   	}

   	if (cost < -1.0)
   	{
       	cost = -2.0 - cost;
   	}

	if (abs(cost) > 2.0) 
	{
        cout<<"Upper or lower limit of cost is excided " << cost<<endl;
		exit(0);
	}

	EulangNew[PHI] = phi;
	EulangNew[CTH] = acos(cost);
	EulangNew[CHI] = 0.0;

   	double sint = sqrt(1.0 - cost*cost);

   	newcoords[AXIS_X][t1] = sint*cos(phi);
   	newcoords[AXIS_Y][t1] = sint*sin(phi);
   	newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// 	the old density

   	double p0 = 0.0;
   	double p1 = 0.0;

   	for (int id=0;id<NDIM;id++)
   	{
      	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      	p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_old;
   	double rho1,rho2,erot;
// 	If it1 = 0 (the first bead), dens_new = SRotDens(p1,type)
// 	if it1 = (NumbRotTimes-1) is the last bead, dens_new = SRotDens(p0,type)

	if(RotDenType == 0)
	{
#ifdef PIGSTYPE
        if (it1 == 0 || it1 == (NumbRotTimes - 1))
        {
            if (it1 == 0)
            {
                dens_old = SRotDens(p1, type);
            }
            else
            {
                dens_old = SRotDens(p0, type);
            }
        }
        else
        {
            dens_old = SRotDens(p0,type)*SRotDens(p1,type);
        }
#else
        dens_old = SRotDens(p0,type)*SRotDens(p1,type);
#endif
	}
    else if(RotDenType == 1)
    {
        rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
        rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
        dens_old = rho1+rho2;
    }

   if (fabs(dens_old)<RZERO) dens_old = 0.0;
#ifndef NEGATIVEDENSITY
   if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
// tapas's temporary treatment for negative rho of two linear rotors
   if (dens_old<0.0) dens_old=fabs(dens_old);
#endif

   double pot_old  = 0.0;

   int itr0 = it1  * RotRatio;     // interval to average over
   int itr1 = itr0 + RotRatio;     // translational time slices

   	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
   		//pot_old  += (PotRotEnergy(gatom,MCCosine,it));
   		pot_old  += (PotRotEnergy(gatom,EulangOld,it));
	}

// the new density 

   p0 = 0.0;
   p1 = 0.0;


   for (int id=0;id<NDIM;id++)
   {
       p0 += (MCCosine[id][t0]*newcoords[id][t1]);
       p1 += (newcoords[id][t1]*MCCosine[id][t2]);
   }

   double dens_new;

	if(RotDenType == 0)
	{
#ifdef PIGSTYPE
        if ((it1 == 0) || (it1 == (NumbRotTimes - 1)))
        {
            if (it1 == 0)
            {
                dens_new = SRotDens(p1, type);
            }
            else
            {
                dens_new = SRotDens(p0, type);
            }
        }
        else
        {
            dens_new = SRotDens(p0,type)*SRotDens(p1,type);
        }
#else
        dens_new = SRotDens(p0,type)*SRotDens(p1,type);
#endif
	}
	else if(RotDenType == 1)
	{
		rsline_(&X_Rot,&p0,&MCRotTau,&rho1,&erot);
		rsline_(&X_Rot,&p1,&MCRotTau,&rho2,&erot);
		dens_new = rho1 + rho2;
	}

	if (fabs(dens_new)<RZERO) dens_new = 0.0;
#ifndef NEGATIVEDENSITY
	if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
	if (dens_new<0.0) dens_new=fabs(dens_new);
#endif

	double pot_new  = 0.0;

	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		//pot_new  += (PotRotEnergy(gatom,newcoords,it));
		pot_new  += (PotRotEnergy(gatom,EulangNew,it));
	}

	double rd;

	if(RotDenType == 0)
	{
		if (dens_old>RZERO)
			rd = dens_new/dens_old;
		else rd = 1.0;

		rd *= exp(- MCTau*(pot_new-pot_old));
	}
	else if(RotDenType == 1)
	{
		rd = dens_new - dens_old - MCTau*(pot_new-pot_old);
		//rd = exp(rd);
	}

	bool Accepted = false;
	if(RotDenType == 0)
	{
		if (rd>1.0)         Accepted = true;
		//else if (rd>rnd7()) Accepted = true;
		else if (rd>rand3) Accepted = true;
	}
	else if (RotDenType == 1)
	{
		if (rd > 0.0)   Accepted = true;
		//else if (rd > log(rnd7())) Accepted = true;
    	else if (rd > log(rand3)) Accepted = true;
	}

	MCRotChunkTot += 1.0;

	if (Accepted)
	{
		MCRotChunkAcp += 1.0;

		MCAngles[CTH][t1] = cost;
		MCAngles[PHI][t1] = phi;

		for (int id=0;id<NDIM;id++)
    		MCCosine[id][t1] = newcoords[id][t1];
	}
#endif
}
#endif

void MCRotLinStepSwapBroken(int it1,int offset,int gatom,int type,double step,double rand1,double rand2,double rand3,double &MCRotChunkTot,double &MCRotChunkAcp)
{
	int it0 = (it1 - 1);
	int it2 = (it1 + 1);

   	if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
   	if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0

   	int t0 = offset + it0;
   	int t1 = offset + it1;
   	int t2 = offset + it2;

   	double cost = MCAngles[CTH][t1];
   	double phi  = MCAngles[PHI][t1];

	double EulangOld[NDIM], EulangNew[NDIM];
	EulangOld[PHI] = phi;
	EulangOld[CTH] = acos(cost);
	EulangOld[CHI] = 0.0;

   	cost += (step*(rand1-0.5));
   	phi  += (step*(rand2-0.5));

   	if (cost >  1.0)
   	{
      	cost = 2.0 - cost;
   	}

   	if (cost < -1.0)
   	{
       	cost = -2.0 - cost;
   	}

	if (abs(cost) > 2.0) 
	{
        cout<<"Upper or lower limit of cost is excided " << cost<<endl;
		exit(0);
	}

	EulangNew[PHI] = phi;
	EulangNew[CTH] = acos(cost);
	EulangNew[CHI] = 0.0;

   	double sint = sqrt(1.0 - cost*cost);

   	newcoords[AXIS_X][t1] = sint*cos(phi);
   	newcoords[AXIS_Y][t1] = sint*sin(phi);
   	newcoords[AXIS_Z][t1] = cost;

//----------------------------------------------

// 	the old density

   	double p0 = 0.0;
   	double p1 = 0.0;

   	for (int id=0;id<NDIM;id++)
   	{
      	p0 += (MCCosine[id][t0]*MCCosine[id][t1]);
      	p1 += (MCCosine[id][t1]*MCCosine[id][t2]);
   	}

   	double dens_old;
// 	If it1 = 0 (the first bead), dens_new = SRotDens(p1,type)
// 	if it1 = (NumbRotTimes-1) is the last bead, dens_new = SRotDens(p0,type)
   	int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	int particleA1Max = particleA1Min + NumbParticle - 1;
   	int particleA2Min = particleA1Max + 1;
   	int particleA2Max = particleA2Min + NumbParticle - 1;

    if (it1 == 0 || it1 == (NumbRotTimes - 1))
    {
        if (it1 == 0)
        {
            dens_old = SRotDens(p1, type);
        }
        else
        {
            dens_old = SRotDens(p0, type);
        }
    }
    else
    {
        dens_old = SRotDens(p0,type)*SRotDens(p1,type);
    }
#ifdef BROKENPATH
    if ((gatom >= particleA1Min) && (gatom <= particleA2Max))
    {
        if (it1 == (((NumbRotTimes - 1)/2) - 1))
        {
            dens_old = SRotDens(p0,type);
        }
        if (it1 == ((NumbRotTimes - 1)/2))
        {
            dens_old = SRotDens(p1,type);
        }
    }
#endif
#ifdef SWAP
    if ((gatom >= particleA1Min) && (gatom <= particleA1Max))
    {
       	int gatomSwap = particleA2Max - (gatom - particleA1Min);
        int offsetSwap = NumbRotTimes*gatomSwap;

        if (it1 == (((NumbRotTimes - 1)/2) - 1))
        {
          	int tSwap = offsetSwap + it2;
			double pSwap = 0.0;
           	for (int id = 0; id < NDIM; id++)
           	{
               	pSwap += (MCCosine[id][t1]*MCCosine[id][tSwap]);
           	}
            dens_old = SRotDens(p0,type)*SRotDens(pSwap,type);
        }
        if (it1 == ((NumbRotTimes - 1)/2))
        {
          	int tSwap = offsetSwap + it0;
			double pSwap = 0.0;
           	for (int id = 0; id < NDIM; id++)
           	{
               	pSwap += (MCCosine[id][tSwap]*MCCosine[id][t1]);
           	}
            dens_old = SRotDens(pSwap,type)*SRotDens(p1,type);
        }
    }

    if ((gatom >= particleA2Min) && (gatom <= particleA2Max))
    {
        int gatomSwap = particleA1Max - (gatom - particleA2Min);
        int offsetSwap = NumbRotTimes*gatomSwap;

        if (it1 == (((NumbRotTimes - 1)/2) - 1))
        {
          	int tSwap = offsetSwap + it2;
			double pSwap = 0.0;
           	for (int id = 0; id < NDIM; id++)
           	{
               	pSwap += (MCCosine[id][t1]*MCCosine[id][tSwap]);
           	}
            dens_old = SRotDens(p0,type)*SRotDens(pSwap,type);
        }
        if (it1 == ((NumbRotTimes - 1)/2))
        {
          	int tSwap = offsetSwap + it0;
			double pSwap = 0.0;
           	for (int id = 0; id < NDIM; id++)
           	{
               	pSwap += (MCCosine[id][tSwap]*MCCosine[id][t1]);
           	}
            dens_old = SRotDens(pSwap,type)*SRotDens(p1,type);
        }
	}
#endif

   if (fabs(dens_old)<RZERO) dens_old = 0.0;
#ifndef NEGATIVEDENSITY
   if (dens_old<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
// tapas's temporary treatment for negative rho of two linear rotors
   if (dens_old<0.0) dens_old=fabs(dens_old);
#endif

   double pot_old  = 0.0;

   int itr0 = it1  * RotRatio;     // interval to average over
   int itr1 = itr0 + RotRatio;     // translational time slices

   	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
   		//pot_old  += (PotRotEnergy(gatom,MCCosine,it));
   		pot_old  += (PotRotEnergySwapBroken(gatom,EulangOld,it));
	}

// the new density 

   p0 = 0.0;
   p1 = 0.0;


   for (int id=0;id<NDIM;id++)
   {
       p0 += (MCCosine[id][t0]*newcoords[id][t1]);
       p1 += (newcoords[id][t1]*MCCosine[id][t2]);
   }

   double dens_new;

   if ((it1 == 0) || (it1 == (NumbRotTimes - 1)))
   {
       if (it1 == 0)
       {
           dens_new = SRotDens(p1, type);
       }
       else
       {
           dens_new = SRotDens(p0, type);
       }
    }
    else
    {
        dens_new = SRotDens(p0,type)*SRotDens(p1,type);
    }
#ifdef BROKENPATH
    if ((gatom >= particleA1Min) && (gatom <= particleA2Max))
    {
        if (it1 == (((NumbRotTimes - 1)/2) - 1))
        {
            dens_new = SRotDens(p0,type);
        }
        if (it1 == ((NumbRotTimes - 1)/2))
        {
            dens_new = SRotDens(p1,type);
        }
    }
#endif
#ifdef SWAP
        if ((gatom >= particleA1Min) && (gatom <= particleA1Max))
        {
            int gatomSwap = particleA2Max - (gatom - particleA1Min);
            int offsetSwap = NumbRotTimes*gatomSwap;

            if (it1 == (((NumbRotTimes - 1)/2) - 1))
            {
                int tSwap = offsetSwap + it2;
                double pSwap = 0.0;
                for (int id = 0; id < NDIM; id++)
                {
                    pSwap += (newcoords[id][t1]*MCCosine[id][tSwap]);
                }
                dens_new = SRotDens(p0,type)*SRotDens(pSwap,type);
            }
            if (it1 == ((NumbRotTimes - 1)/2))
            {
                int tSwap = offsetSwap + it0;
                double pSwap = 0.0;
                for (int id = 0; id < NDIM; id++)
                {
                    pSwap += (MCCosine[id][tSwap]*newcoords[id][t1]);
                }
                dens_new = SRotDens(pSwap,type)*SRotDens(p1,type);
            }
        }

        if ((gatom >= particleA2Min) && (gatom <= particleA2Max))
        {
            int gatomSwap = particleA1Max - (gatom - particleA2Min);
            int offsetSwap = NumbRotTimes*gatomSwap;

            if (it1 == (((NumbRotTimes - 1)/2) - 1))
            {
                int tSwap = offsetSwap + it2;
                double pSwap = 0.0;
                for (int id = 0; id < NDIM; id++)
                {
                    pSwap += (newcoords[id][t1]*MCCosine[id][tSwap]);
                }
                dens_new = SRotDens(p0,type)*SRotDens(pSwap,type);
            }
            if (it1 == ((NumbRotTimes - 1)/2))
            {
                int tSwap = offsetSwap + it0;
                double pSwap = 0.0;
                for (int id = 0; id < NDIM; id++)
                {
                    pSwap += (MCCosine[id][tSwap]*newcoords[id][t1]);
                }
                dens_new = SRotDens(pSwap,type)*SRotDens(p1,type);
            }
        }
#endif

	if (fabs(dens_new)<RZERO) dens_new = 0.0;
#ifndef NEGATIVEDENSITY
	if (dens_new<0.0 && RotDenType == 0) nrerror("Rotational Moves: ","Negative rot density");
#else
	if (dens_new<0.0) dens_new=fabs(dens_new);
#endif

	double pot_new  = 0.0;

	for (int it=itr0;it<itr1;it++)  // average over tr time slices
	{
		//pot_new  += (PotRotEnergy(gatom,newcoords,it));
		pot_new  += (PotRotEnergySwapBroken(gatom,EulangNew,it));
	}

	double rd;
	if (dens_old>RZERO)
		rd = dens_new/dens_old;
	else rd = 1.0;

	rd *= exp(- MCTau*(pot_new-pot_old));

	bool Accepted = false;
	if (rd>1.0)         Accepted = true;
	else if (rd>rand3) Accepted = true;

	MCRotChunkTot += 1.0;

	if (Accepted)
	{
		MCRotChunkAcp += 1.0;

		MCAngles[CTH][t1] = cost;
		MCAngles[PHI][t1] = phi;

		for (int id=0;id<NDIM;id++)
    		MCCosine[id][t1] = newcoords[id][t1];
	}
}

//double PotRotEnergy(int atom0, double **cosine, int it)   
double PotRotEnergySwapBroken(int atom0, double *Eulang0, int it)   
{
	int type0   =  MCType[atom0];
	double spot, spotSwap;

    double weight, weight1;
	weight = 1.0;

    if (it == 0 || it == (NumbRotTimes - 1)) weight = 0.5;

	if ( (MCAtom[type0].molecule == 4) && (MCAtom[type0].numb > 1) )
	{
	    int offset0 =  atom0*NumbRotTimes;
        int t0  = offset0 + it;

		int atom1Init, NumbAtoms1;
   	    int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	    int particleA1Max = particleA1Min + NumbParticle - 1;
   	    int particleA2Min = particleA1Max + 1;
   	    int particleA2Max = particleA2Min + NumbParticle - 1;

        if (atom0 <= particleA1Max)
        {
            atom1Init  = 0;
            NumbAtoms1 = (particleA1Max+1);
        }
        else
        {
            atom1Init  = particleA2Min;
            NumbAtoms1 = NumbAtoms;
        }

        spot = 0.0;
        for (int atom1 = atom1Init; atom1 < NumbAtoms1; atom1++)
        if (atom1 != atom0)                    
        {
            int offset1 = atom1*NumbRotTimes;
            int t1  = offset1 + it;

	        string stype = MCAtom[type0].type;
		    if (stype == HF )
            {
				double Eulang1[NDIM];
				Eulang1[PHI] = MCAngles[PHI][t1];
                Eulang1[CTH] = acos(MCAngles[CTH][t1]);
                Eulang1[CHI] = 0.0;

				weight1 = 1.0;
				if (((atom0 < particleA1Min) || (atom0 > particleA2Max)) && ((atom1 >= particleA1Min) && (atom1 <= particleA2Max)))
            	{
           			if (it == ((NumbRotTimes - 1)/2))
              		{
                   		weight1 = 0.5;
					}
            	} 
            	if (((atom0 >= particleA1Min) && (atom0 <= particleA2Max)) && ((atom1 < particleA1Min) || (atom1 > particleA2Max)))
            	{
            		if (it == ((NumbRotTimes - 1)/2))
               		{
                   		weight1 = 0.5;
					}
            	} 
        		spot += weight*weight1*PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            }  //stype
        } //loop over atom1 (molecules)
		spotSwap = 0.0;
#ifdef SWAP
		if (it == ((NumbRotTimes - 1)/2))
		{
			int atomSwapMin, atomSwapMax;
			if (atom0 < particleA1Min)
    	   	{
				atomSwapMin = particleA2Min;
                atomSwapMax = (particleA2Max + 1);
			}
			if (atom0 > particleA2Max)
    	   	{
				atomSwapMin = particleA1Min;
                atomSwapMax = (particleA1Max + 1);
			}
	    	if ((atom0 >= particleA1Min) && (atom0 <= particleA1Max))
       	    {
			    atomSwapMin = (particleA2Max + 1);
                atomSwapMax = NumbAtoms;
		    }
     		if ((atom0 >= particleA2Min) && (atom0 <= particleA2Max))
         	{
	     		atomSwapMin = 0;
                atomSwapMax = particleA1Min;
		    }
			spotSwap = 0.0;
     		for (int atomSwap = atomSwapMin; atomSwap < atomSwapMax; atomSwap++)
	    	{
                int offsetSwap = NumbRotTimes*atomSwap;
                int tSwap  = offsetSwap + it;
	        	string stype = MCAtom[type0].type;
                if (stype == HF )
                {
					double EulangSwap[NDIM];
					EulangSwap[PHI] = MCAngles[PHI][tSwap];
        			EulangSwap[CTH] = acos(MCAngles[CTH][tSwap]);
        			EulangSwap[CHI] = 0.0;
        			spotSwap += 0.5*PotFunc(atom0, atomSwap, Eulang0, EulangSwap, it);
                }  //stype
		    }
        }
#endif
    }
	double spotReturn = (spot + spotSwap);
    return spotReturn;
}
void MCRot3Dstep(int it1, int offset, int gatom, int type, double step,double rand1,double rand2,double rand3,double rand4,int IROTSYM, int NFOLD_ROT,double &MCRotChunkTot,double &MCRotChunkAcp)
{
      int it0 = (it1 - 1);
      int it2 = (it1 + 1);
 
      if (it0<0)             it0 += NumbRotTimes; // NumbRotTimes - 1
      if (it2>=NumbRotTimes) it2 -= NumbRotTimes; // 0
      
      int t0 = offset + it0;
      int t1 = offset + it1;
      int t2 = offset + it2;

      double cost = MCAngles[CTH][t1];
      double phi  = MCAngles[PHI][t1];
      double chi  = MCAngles[CHI][t1];

//    cost += (step*(rnd1()-0.5));
      cost += (step*(rand1-0.5));

//    Toby change:
//    cout<<"before random change "<<phi<<" "<<chi<<" "<<endl;
//    phi += 2.0*M_PI*(step*(rnd1()-0.5));
//    chi += 2.0*M_PI*(step*(rnd1()-0.5));
      phi += 2.0*M_PI*(step*(rand2-0.5));
      chi += 2.0*M_PI*(step*(rand3-0.5));
/*
//    axial symmetry of the molecule controlled by IROTSYM and NFOLD_ROT. use rand1 to judge whether rotate or not
      if( IROTSYM == 1 )
      {
//       if(rand1 < 1.0/3.0)
//       chi += 2.0*M_PI/(double)NFOLD_ROT;

         if(rand1 < 2.0/3.0 && rand1 >= 1.0/3.0)
         chi += 2.0*M_PI/(double)NFOLD_ROT;

         if(rand1 >= 2.0/3.0 )
         chi -= 2.0*M_PI/(double)NFOLD_ROT;
      }
*/

//    get to the positive values of phi and chi
      if(phi<0.0) phi = 2.0*M_PI + phi;
      if(chi<0.0) chi = 2.0*M_PI + chi;
//    Toby needs to recover the [0:2*Pi] range for phi and chi
      phi = fmod(phi,2.0*M_PI);
      chi = fmod(chi,2.0*M_PI);

      if (cost >  1.0)
      {
         cost = 2.0 - cost;    
//       phi  = phi + M_PI;
      }		 
      
      if (cost < -1.0)
      {
          cost = -2.0 - cost;    
//        phi  = phi  + M_PI;
      }		  

      double sint = sqrt(1.0 - cost*cost);

      newcoords[PHI][t1] = phi;
      newcoords[CHI][t1] = chi;
      newcoords[CTH][t1] = cost;

//----------------------------------------------

// the old density

      double rho = 0.0;
      double erot = 0.0;
      double esq  = 0.0;
      double Eulan1[3];
      double Eulan2[3];
      double Eulrel[3];
      int istop=0;
      Eulan1[0]=MCAngles[PHI][t0];
      Eulan1[1]=acos(MCAngles[CTH][t0]);
      Eulan1[2]=MCAngles[CHI][t0];
      Eulan2[0]=MCAngles[PHI][t1];
      Eulan2[1]=acos(MCAngles[CTH][t1]);
      Eulan2[2]=MCAngles[CHI][t1];

      if(RotDenType == 0)
      {
         rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);

         if(istop == 1)
         {
          cerr<<"large matrix test error"<<endl;
          exit(0);
         }
      }
      else if(RotDenType == 1)
      {
         rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);

      }
      double dens_old = rho;
      double rhoold = rho;

      Eulan1[0]=MCAngles[PHI][t1];
      Eulan1[1]=acos(MCAngles[CTH][t1]);
      Eulan1[2]=MCAngles[CHI][t1];
      Eulan2[0]=MCAngles[PHI][t2];
      Eulan2[1]=acos(MCAngles[CTH][t2]);
      Eulan2[2]=MCAngles[CHI][t2];

      istop=0;
      if(RotDenType == 0)
      {
         rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);
         if(istop == 1)
         {
          cerr<<"large matrix test error"<<endl;
          exit(0);
         }
      }
      else if(RotDenType == 1)
      {
         rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);
      }
// PN pigs
// ask Toby if rhoold is for PQC only?
#ifdef PIGSTYPE
	if (t1 == 0 || t1 == (NumbRotTimes-1))
   		if (t1 == 0)
		{
       		dens_old = rho;
        	rhoold   = rho;
		}
    	else 
		{
       		dens_old = dens_old; 
        	rhoold   =  rhoold;
		}
   	else 
	{
#endif
   	dens_old         = dens_old*rho;
   	rhoold           = rhoold + rho;
#ifdef PIGSTYPE
   	}
#endif

      if (fabs(dens_old)<RZERO) dens_old = 0.0;
//    if (dens_old<0.0) nrerror("Rotational Moves: ","Negative rot density");
//    toby's temporary treatment for negative rho of paraH2O
      if (dens_old<0.0) dens_old=fabs(dens_old);

      double pot_old  = 0.0;

      int itr0 = it1  * RotRatio;     // interval to average over
      int itr1 = itr0 + RotRatio;     // translational time slices

#ifndef FREE_ROTORS
      for (int it=itr0;it<itr1;it++)  // average over tr time slices
      pot_old  += (PotRotE3D(gatom,Eulan1,it));
#endif
//    Toby: pot_old can be calculated with MCAngles

//   the new density 

      Eulan1[0]=MCAngles[PHI][t0];
      Eulan1[1]=acos(MCAngles[CTH][t0]);
      Eulan1[2]=MCAngles[CHI][t0];
      Eulan2[0]=newcoords[PHI][t1];
      Eulan2[1]=acos(newcoords[CTH][t1]);
      Eulan2[2]=newcoords[CHI][t1];

      istop=0;
      if(RotDenType == 0)
      {
         rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);
         if(istop == 1)
         {
          cerr<<"large matrix test error"<<endl;
          exit(0);
         }
      }
      else if(RotDenType == 1)
      {
         rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);
      }

      double dens_new = rho;
      double rhonew = rho;

      Eulan1[0]=newcoords[PHI][t1];
      Eulan1[1]=acos(newcoords[CTH][t1]);
      Eulan1[2]=newcoords[CHI][t1];
      Eulan2[0]=MCAngles[PHI][t2];
      Eulan2[1]=acos(MCAngles[CTH][t2]);
      Eulan2[2]=MCAngles[CHI][t2];

      istop=0;
      if(RotDenType == 0)
      {
         rotden_(Eulan1,Eulan2,Eulrel,&rho,&erot,&esq,rhoprp,erotpr,erotsq,&istop);
         if(istop == 1)
         {
          cerr<<"large matrix test error"<<endl;
          exit(0);
         }
      }
      else if(RotDenType == 1)
      {
         rsrot_(Eulan1,Eulan2,&X_Rot,&Y_Rot,&Z_Rot,&MCRotTau,&RotOdEvn,&RotEoff,&rho,&erot);
      }

#ifdef PIGSTYPE
	if (t1 == 0 || t1 == (NumbRotTimes-1))
   	if (t1 == 0)
	{
       	dens_new = rho;
        rhonew   = rho;
	}
    else 
	{
       	dens_new = dens_new; 
        rhonew   =  rhonew;
	}
   	else 
	{
#endif
   	dens_new     =dens_new*rho;
   	rhonew       = rhonew + rho;
#ifdef PIGSTYPE
   	}
#endif

      if (fabs(dens_new)<RZERO) dens_new = 0.0;
//    if (dens_new<0.0) nrerror("Rotational Moves: ","Negative rot density");
//    toby's temporary treatment for negative rho of paraH2O
      if (dens_new<0.0) dens_new=fabs(dens_new);


      double pot_new  = 0.0;

#ifndef FREE_ROTORS
      for (int it=itr0;it<itr1;it++)  // average over tr time slices
      pot_new  += (PotRotE3D(gatom,Eulan1,it));
#endif
//    Toby: pot_new can be calculated with newcoords
     
      double rd;

//    distinginsh between Noya and RS
	if(RotDenType == 0)
    {
    	if (dens_old>RZERO)
        rd = dens_new/dens_old;
        else rd = 1.0;
#ifndef FREE_ROTORS
        rd *= exp(- MCTau*(pot_new-pot_old));
#endif
    }
    else if(RotDenType == 1)
   	{
		//rd = dens_new/dens_old;
		//use logarithmic for RS
		rd = (rhonew - rhoold)/(4.0*(MCRotTau/WNO2K));
		//cout<<"in cc:"<<rhonew<<" "<<rhoold<<" "<<rd<<" "<<4.0*(MCRotTau/WNO2K)<<endl;
		//rd = rhonew - rhoold;
		//rd = exp(rd);
		rd -= MCTau*(pot_new-pot_old);
	}

	//rd *= exp(- MCTau*(pot_new-pot_old));   

	bool Accepted = false;
	if(RotDenType == 0)
	{
		if (rd>1.0)         Accepted = true;
		//else if (rd>rnd7()) Accepted = true;
		else if (rd>rand4) Accepted = true;
    }
    else
    {
    	if (rd > 0.0)       Accepted = true;
		//else if (rd > log(rnd7())) Accepted = true;
		else if (rd > log(rand4)) Accepted = true;
	}

	//MCTotal[type][MCROTAT] += 1.0;  
	MCRotChunkTot +=1.0;
      
	if (Accepted)
    {
		//MCAccep[type][MCROTAT] += 1.0;
		MCRotChunkAcp +=1.0;

    	MCAngles[CTH][t1] = cost;
    	MCAngles[PHI][t1] = phi;
    	MCAngles[CHI][t1] = chi; //toby adds
  
    	sint=sqrt(1.0-cost*cost);
    	MCCosine [AXIS_X][t1] = sint*cos(phi);
    	MCCosine [AXIS_Y][t1] = sint*sin(phi);
    	MCCosine [AXIS_Z][t1] = cost;
		//This MCCosine will be used in estimating correlation function of the orientation of one molecule-fixed axis in GetRCF
		//and Ieff about and perpendicular to one molecule-ixed axis.

#ifdef MOLECULEINCAGE
		MCCosinex[AXIS_X][t1] = cost*cos(phi)*cos(chi)-sin(phi)*sin(chi);
        MCCosinex[AXIS_Y][t1] = cost*sin(phi)*cos(chi)+cos(phi)*sin(chi);
        MCCosinex[AXIS_Z][t1] = -sint*cos(chi);

        MCCosiney[AXIS_X][t1] = -cost*cos(phi)*sin(chi)-sin(phi)*cos(chi);
        MCCosiney[AXIS_Y][t1] = -cost*sin(phi)*sin(chi)+cos(phi)*cos(chi);
        MCCosiney[AXIS_Z][t1] = sint*sin(chi);
#endif
	}	      
}

double PotEnergy(int atom0, double **pos)   
// 
//  interaction of atom0 with other atoms/molecules
//  only two atom types so far, with the number of second particles 0 or 1 
//
{
	int type0   = MCType[atom0];
   	int offset0 = NumbTimes*atom0;

   	double dr[NDIM];
   	double spot =  0.0;

   	for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
   	if (atom1 != atom0)                      // skip "self-interaction"
   	{	    
       	int type1   = MCType[atom1];
       	int offset1 = NumbTimes*atom1; 

       	double spot_pair=0.0;

       	#pragma omp parallel for reduction(+: spot_pair)
       	for (int it = 0; it < NumbTimes; it++) 	    
       	{ 
      		int t0 = offset0 + it;
       		int t1 = offset1 + it;

	        string stype = MCAtom[type0].type;
			if (stype == H2)
	        {
                double s1 = 0.0;
                double s2 = 0.0;
                double dr2 = 0.0;
				double dr[NDIM];

                for (int id = 0; id < NDIM; id++)
                {
                    dr[id]  = (pos[id][t0] - MCCoords[id][t1]);
                    dr2    += (dr[id]*dr[id]);
                    double cst1 = (MCCoords[id][t1] - pos[id][t0])*MCCosine[id][t0];
                    double cst2 = (MCCoords[id][t1] - pos[id][t0])*MCCosine[id][t1];
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
                    b2[id] = (MCCoords[id][t1] - pos[id][t0])/r;
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
                double r2 = r1;
                double rd = r/BOHRRADIUS;
                double potl;
                vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
                spot_pair += potl*CMRECIP2KL;
			}  //stype

#ifdef IOWRITE
           	bool wline = true;                  // skip if the time slice between ira and masha

          	if (WORM && Worm.exists && (Worm.type == type1))  
          	wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
          
          	if (wline)
          	{
          		int t0 = offset0 + it;
          		int t1 = offset1 + it;

          		double dr2 = 0.0;  		 
          		for (int id=0;id<NDIM;id++)
	  			{
             		dr[id]  = (pos[id][t0] - MCCoords[id][t1]);
            
             		if (MINIMAGE)
             		dr[id] -= (BoxSize*rint(dr[id]/BoxSize));

             		dr2    += (dr[id]*dr[id]);
          		}
	       	 
//#ifdef _CUTOFF_	     
//       		if (dr2<dljcutoff2)
//#endif
          		double r = sqrt(dr2);

//-------------- MOLECULES ----------------------

          		int tm;

          		if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
          		{ 
              		int sgn = 1;                // set to -1 to correct the orientaion of dr

              		tm = offset1 + it/RotRatio;

              		int typep = type1;           // define the model of interaction

              		if (MCAtom[type0].molecule == 1)  // does not work for two molecules
              		{
                  		sgn = -1;   

                  		tm  = offset0 + it/RotRatio;

                  		typep = type0; 
              		}

              		double cost = 0.0;
              		for (int id=0;id<NDIM;id++) // n*dr = r*cos(theta) 
              		cost += (MCCosine[id][tm]*dr[id]);   	 
	 
              		cost /= r;                  // cos(theta)
              		cost *= sgn;                // correct the orientation 

              		spot_pair += LPot2D(r,cost,typep);   
          		}
//----------------ATOM-NON-LINEAR MOLECULES----------------
          		else if ((((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && ISPHER == 0) && (MCAtom[type0].molecule != MCAtom[type1].molecule)) // 3D interaction
          		{
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
              		int    ivcord=0;
              		if(MCAtom[type0].molecule == 2)
              		{
                 		tm  = offset0 + it/RotRatio;
                 		for (int id=0;id<NDIM;id++)
                 		{
                    		RCOM[id] = pos[id][t0];
                    		Rpt[id]  = MCCoords[id][t1];
                 		}
              		}
              		else
              		{
                 		tm  = offset1 + it/RotRatio;
                 		for (int id=0;id<NDIM;id++)
                 		{
                    		Rpt[id]  = pos[id][t0];
                    		RCOM[id] = MCCoords[id][t1];
                 		}
              		}
              		Eulang[PHI]=MCAngles[PHI][tm];
              		Eulang[CTH]=acos(MCAngles[CTH][tm]);
              		Eulang[CHI]=MCAngles[CHI][tm];

              		vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

              		spot_pair += vpot3d;

          		}
//----------------ATOM-NON-LINEAR MOLECULES spherical----------------
          		else if ((((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && ISPHER == 1) && (MCAtom[type0].molecule != MCAtom[type1].molecule)) // spherical treatment for non-linear rotor
          		{

              		double radret,vpot3d;
              		radret = r;
              		vspher_(&radret,&vpot3d);

              		spot_pair += vpot3d;

          		}
//----------------- NonLinear ---- NonLinear------------------
          		else if ( ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
          		{
        		// GG:
				//           cout<<"PotEnergy: ((MCAtom[type0].molecule == 2) && (MCAtom[type1].molecule == 2))"<<endl;
      			//    if (MCType[atom1] == IMTYPE)
        		//     {
        		//      int t0 = offset0 + it;
        //      	int t1 = offset1 + it;
              		double com_1[3];
              		double com_2[3];
              		double Eulang_1[3];
              		double Eulang_2[3];
              		double E_2H2O;
              		int t0 = offset0 + it;
              		int t1 = offset1 + it;
              		for (int id=0;id<NDIM;id++)
              		{
             //  		cout<<"id it pos[id][t0] "<<id<<" "<<it<<" "<<pos[id][t0]<<endl;
             //  		cout<<"id it MCCoords[id][t1] "<<id<<" "<<it<<" "<<MCCoords[id][t1]<<endl;
                   		com_1[id] = pos[id][t0];
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
              		caleng_(com_1, com_2, &E_2H2O, Eulang_1, Eulang_2);
//            		cout<<t0<<" "<<t1<<" "<<com_1[0]<<" "<<com_1[1]<<" "<<com_1[2]<<" "<<com_2[0]<<" "<<com_2[1]<<" "<<com_2[2]<<endl;
//            		cout<<Eulang_1[0]<<" "<<Eulang_1[1]<<" "<<Eulang_1[2]<<" "<<Eulang_2[0]<<" "<<Eulang_2[1]<<" "<<Eulang_2[2]<<" "<<E_2H2O<<endl;
              		spot_pair += E_2H2O;
           //   }
          		}
//----------------------------------------------- 
          		else 
          		spot_pair += SPot1D(r,type1);    // 1D interaction

// it shoud be SPot1D(r,type0,type1) or  SPot1D(r,ind) with ind =type0*NumbTypes+type1
       		} // wline  
#endif
       	} 
      	spot += spot_pair;
    } 
    return (spot);
}

void Reflect_MF_XZ(void)
{
// reflect the coordinates of point like particles with respect to the xz plane in the MFF of HCOOCH3
// cout<<"in Reflect_MF_XZ"<<" "<<PrintYrfl<<endl;

// print out coordinates before reflection

   if(PrintYrfl)
   {
      string fname="before_refY";
      IOxyzAng(IOWrite,fname.c_str());
   }


///*
// calculate kinetic and potential energy
   double kin_before= GetKinEnergy();
   double pot_before =GetPotEnergy();
   double rot_before = GetRotE3D();
// cout<<"before:"<<kin<<" "<<spot<<" "<<srot<<endl;
//*/

// reflect Euler angles for all rotors
   for (int molec=0;molec<MCAtom[IMTYPE].numb;molec++)
   {

      int offset = MCAtom[IMTYPE].offset + molec*NumbTimes;
      #pragma omp parallel for
      for (int it_rot=0;it_rot<NumbTimes/RotRatio;it_rot++)
      {

         int pMF = it_rot+offset;

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
//          Rpt[id]  = MCCoords[id][it_rot];
            RCOM[id] = MCCoords[id][pMF];
         }

         Eulang[PHI]=MCAngles[PHI][pMF];
         Eulang[CTH]=acos(MCAngles[CTH][pMF]);
         Eulang[CHI]=MCAngles[CHI][pMF];

         vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

/*
         cout<<"in C"<<endl;
         for (int dim=0;dim<NDIM;dim++)
         {
            cout<<RCOM[dim]<<" "<<hatx[dim]<<" "<<haty[dim]<<" "<<hatz[dim]<<" "<<Eulang[dim]<<endl;
         }
*/

         rflmfy_(RCOM,hatx,haty,hatz,Eulang);

         MCAngles[PHI][pMF]=Eulang[PHI];
         MCAngles[CTH][pMF]=cos(Eulang[CTH]);
         MCAngles[CHI][pMF]=Eulang[CHI];
//       cout<<"in C:"<<MCAngles[PHI][pMF]<<" "<<MCAngles[CTH][pMF]<<" "<<MCAngles[CHI][pMF]<<endl;

      }
   }

// reflect positions
   for (int atype=0;atype<NumbTypes;atype++)
   for (int atom=0;atom<MCAtom[atype].numb;atom++)
   {
   #pragma omp parallel for
   for (int it=0;it<NumbTimes;it++)
   {

      int offset=MCAtom[atype].offset+NumbTimes*atom;
      MCCoords[1][offset+it] *=-1.0;

   }
   }

// print out coordinates after reflection
   if(PrintYrfl)
   {
     string fname="after_refY";
     IOxyzAng(IOWrite,fname.c_str());
     PrintYrfl = 0;
   }

// calculate kinetic and potential energy
///*
   double kin_after= GetKinEnergy();
   double pot_after =GetPotEnergy();
   double rot_after = GetRotE3D();
// cout<<"after:"<<skin<<" "<<spot<<" "<<srot<<endl;
   if((abs(kin_after-kin_before) > 0.01) || (abs(pot_after-pot_before) > 0.01) || (abs(rot_after-rot_before) > 0.05))
   cout<<"Warning in Reflect_MF_XZ"<<kin_before<<" "<<kin_after<<" "<<pot_before<<" "<<pot_after<<" "<<rot_before<<" "<<rot_after<<endl;
//*/

}

void Reflect_MF_YZ(void)
{
// reflect the coordinates of point like particles with respect to the yz plane in the MFF
// cout<<"in Reflect_MF_YZ"<<" "<<PrintYrfl<<endl;

// print out coordinates before reflection

   if(PrintXrfl)
   {
      string fname="before_refX";
      IOxyzAng(IOWrite,fname.c_str());
   }


///*
// calculate kinetic and potential energy
   double kin_before= GetKinEnergy();
   double pot_before =GetPotEnergy();
   double rot_before = GetRotE3D();
// cout<<"before:"<<kin<<" "<<spot<<" "<<srot<<endl;
//*/

// reflect Euler angles for all rotors
   for (int molec=0;molec<MCAtom[IMTYPE].numb;molec++)
   {
      int offset = MCAtom[IMTYPE].offset + molec*NumbTimes;
      #pragma omp parallel for
      for (int it_rot=0;it_rot<NumbTimes/RotRatio;it_rot++)
      {

         int pMF = it_rot+offset;

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
//          Rpt[id]  = MCCoords[id][it_rot];
            RCOM[id] = MCCoords[id][pMF];
         }

         Eulang[PHI]=MCAngles[PHI][pMF];
         Eulang[CTH]=acos(MCAngles[CTH][pMF]);
         Eulang[CHI]=MCAngles[CHI][pMF];

         vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

/*
         cout<<"in C"<<endl;
         for (int dim=0;dim<NDIM;dim++)
         {
            cout<<RCOM[dim]<<" "<<hatx[dim]<<" "<<haty[dim]<<" "<<hatz[dim]<<" "<<Eulang[dim]<<endl;
         }
*/

         rflmfx_(RCOM,hatx,haty,hatz,Eulang);

         MCAngles[PHI][pMF]=Eulang[PHI];
         MCAngles[CTH][pMF]=cos(Eulang[CTH]);
         MCAngles[CHI][pMF]=Eulang[CHI];
//       cout<<"in C:"<<MCAngles[PHI][pMF]<<" "<<MCAngles[CTH][pMF]<<" "<<MCAngles[CHI][pMF]<<endl;

      }
   }

// reflect positions
   for (int atype=0;atype<NumbTypes;atype++)
   for (int atom=0;atom<MCAtom[atype].numb;atom++)
   {
   #pragma omp parallel for
   for (int it=0;it<NumbTimes;it++)
   {

      int offset=MCAtom[atype].offset+NumbTimes*atom;
      MCCoords[1][offset+it] *=-1.0;

   }
   }

// print out coordinates after reflection
   if(PrintXrfl)
   {
     string fname="after_refX";
     IOxyzAng(IOWrite,fname.c_str());
     PrintXrfl = 0;
   }

// calculate kinetic and potential energy
///*
   double kin_after= GetKinEnergy();
   double pot_after =GetPotEnergy();
   double rot_after = GetRotE3D();
// cout<<"after:"<<skin<<" "<<spot<<" "<<srot<<endl;
   if((abs(kin_after-kin_before) > 0.01) || (abs(pot_after-pot_before) > 0.01) || (abs(rot_after-rot_before) > 0.05))
   cout<<"Warning in Reflect_MF_YZ"<<kin_before<<" "<<kin_after<<" "<<pot_before<<" "<<pot_after<<" "<<rot_before<<" "<<rot_after<<endl;
//*/

}

void Reflect_MF_XY(void)
{
// reflect the coordinates of point like particles with respect to the xz plane in the MFF of HCOOCH3
// cout<<"in Reflect_MF_XY"<<" "<<PrintZrfl<<endl;

// print out coordinates before reflection

   if(PrintZrfl)
   {
      string fname="before_refZ";
      IOxyzAng(IOWrite,fname.c_str());
   }


///*
// calculate kinetic and potential energy
   double kin_before= GetKinEnergy();
   double pot_before =GetPotEnergy();
   double rot_before = GetRotE3D();
// cout<<"before:"<<kin<<" "<<spot<<" "<<srot<<endl;
//*/

// reflect Euler angles for all rotors
   for (int molec=0;molec<MCAtom[IMTYPE].numb;molec++)
   {
      int offset = MCAtom[IMTYPE].offset + molec*NumbTimes;
      #pragma omp parallel for
      for (int it_rot=0;it_rot<NumbTimes/RotRatio;it_rot++)
      {

         int pMF = it_rot+offset;

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
//          Rpt[id]  = MCCoords[id][it_rot];
            RCOM[id] = MCCoords[id][pMF];
         }

         Eulang[PHI]=MCAngles[PHI][pMF];
         Eulang[CTH]=acos(MCAngles[CTH][pMF]);
         Eulang[CHI]=MCAngles[CHI][pMF];

         vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

/*
         cout<<"in C"<<endl;
         for (int dim=0;dim<NDIM;dim++)
         {
            cout<<RCOM[dim]<<" "<<hatx[dim]<<" "<<haty[dim]<<" "<<hatz[dim]<<" "<<Eulang[dim]<<endl;
         }
*/

         rflmfz_(RCOM,hatx,haty,hatz,Eulang);

         MCAngles[PHI][pMF]=Eulang[PHI];
         MCAngles[CTH][pMF]=cos(Eulang[CTH]);
         MCAngles[CHI][pMF]=Eulang[CHI];
//       cout<<"in C:"<<MCAngles[PHI][pMF]<<" "<<MCAngles[CTH][pMF]<<" "<<MCAngles[CHI][pMF]<<endl;

      }
   }

// reflect positions
   for (int atype=0;atype<NumbTypes;atype++)
   for (int atom=0;atom<MCAtom[atype].numb;atom++)
   {
   #pragma omp parallel for
   for (int it=0;it<NumbTimes;it++)
   {

      int offset=MCAtom[atype].offset+NumbTimes*atom;
      MCCoords[1][offset+it] *=-1.0;

   }
   }

// print out coordinates after reflection
   if(PrintZrfl)
   {
     string fname="after_refZ";
     IOxyzAng(IOWrite,fname.c_str());
     PrintZrfl = 0;
   }

// calculate kinetic and potential energy
///*
   double kin_after= GetKinEnergy();
   double pot_after =GetPotEnergy();
   double rot_after = GetRotE3D();
// cout<<"after:"<<skin<<" "<<spot<<" "<<srot<<endl;
   if((abs(kin_after-kin_before) > 0.01) || (abs(pot_after-pot_before) > 0.01) || (abs(rot_after-rot_before) > 0.05))
   cout<<"Warning in Reflect_MF_XY"<<kin_before<<" "<<kin_after<<" "<<pot_before<<" "<<pot_after<<" "<<rot_before<<" "<<rot_after<<endl;
//*/

}

void RotSymConfig(void)
{
// increase all the body fixed angle of the rotors by a symmetry angle in property evaluation
// cout<<"in RotSymConfig"<<endl;

///*
// calculate kinetic and potential energy
   double kin_before= GetKinEnergy();
   double pot_before =GetPotEnergy();
   double rot_before;
   if(MCAtom[IMTYPE].molecule == 2)
   {
      rot_before = GetRotE3D();
   }
   else if(MCAtom[IMTYPE].molecule == 1)
   {
      rot_before = GetRotEnergy(); 
   }
// cout<<"before:"<<kin<<" "<<spot<<" "<<srot<<endl;
//*/
   double rand=rnd1();

// reflect Euler angles for one arbitrarily chosen rotor
   for (int molec=0;molec<MCAtom[IMTYPE].numb;molec++)
   {
      if(rand > (double)molec/(double)MCAtom[IMTYPE].numb && rand <= (double)(molec+1)/(double)MCAtom[IMTYPE].numb)
      {
         int offset = MCAtom[IMTYPE].offset + molec*NumbTimes;
         #pragma omp parallel for
         for (int it_rot=0;it_rot<NumbTimes/RotRatio;it_rot++)
         {

            int pMF = it_rot+offset;

            if(MCAtom[IMTYPE].molecule == 2) // non-linear rotor
            {
               double chi=MCAngles[CHI][pMF] + 2.0*M_PI/(double)NFOLD_ROT;

               chi = fmod(chi,2.0*M_PI);
               if(chi<0.0) chi = 2.0*M_PI + chi;

               MCAngles[CHI][pMF] = chi;
            }
///*
            else if(MCAtom[IMTYPE].molecule == 1) // linear rotor, regardless of NFOLD_ROT
            {
               double phi  =  MCAngles[PHI][pMF] + M_PI;
               MCAngles[CTH][pMF] *= -1.0;

               phi = fmod(phi,2.0*M_PI);
               if(phi<0.0) phi = 2.0*M_PI + phi;
               MCAngles[PHI][pMF] = phi;
               double cost=MCAngles[CTH][pMF];
               double sint=sqrt(1.0 - cost*cost);
              MCCosine[AXIS_X][pMF] = sint*cos(phi);
              MCCosine[AXIS_Y][pMF] = sint*sin(phi);
              MCCosine[AXIS_Z][pMF] = cost;

            }
//*/

         }
      }
   }


// calculate kinetic and potential energy
///*
   double kin_after= GetKinEnergy();
   double pot_after =GetPotEnergy();
   double rot_after;
   if(MCAtom[IMTYPE].molecule == 2)
   {
      rot_after = GetRotE3D();
   }
   else if(MCAtom[IMTYPE].molecule == 1)
   {
      rot_after = GetRotEnergy();
   }
// cout<<"after:"<<skin<<" "<<spot<<" "<<srot<<endl;
   if((abs(kin_after-kin_before) > 0.01) || (abs(pot_after-pot_before) > 0.01) || (abs(rot_after-rot_before) > 0.05))
   cout<<"Warning in RotSymConfig"<<kin_before<<" "<<kin_after<<" "<<pot_before<<" "<<pot_after<<" "<<rot_before<<" "<<rot_after<<endl;
//*/

}

double PotEnergy(int atom0, double **pos, int it)   
// 
//  interaction of atom0 with other atoms/molecules
//  only two atom types so far, with number of second particles 0 or 1 
//
{
   	int type0   = MCType[atom0];
   	int offset0 = NumbTimes*atom0;

   	double dr[NDIM];
   	double spot = 0.0;

   	for (int atom1=0;atom1<NumbAtoms;atom1++)
   	if (atom1 != atom0)                    // skip "self-interaction"
   	{	
     	int type1   = MCType[atom1];
     	int offset1 = NumbTimes*atom1; 
        int t0 = offset0 + it;
        int t1 = offset1 + it;

	    string stype = MCAtom[type0].type;
		if (stype == H2)
	    {
   			double s1 = 0.0;
            double s2 = 0.0;
            double dr2 = 0.0;
			double dr[NDIM];

            for (int id = 0; id < NDIM; id++)
            {
                dr[id]  = (pos[id][t0] - MCCoords[id][t1]);
                dr2    += (dr[id]*dr[id]);
                double cst1 = (MCCoords[id][t1] - pos[id][t0])*MCCosine[id][t0];
                double cst2 = (MCCoords[id][t1] - pos[id][t0])*MCCosine[id][t1];
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
               	b2[id] = (MCCoords[id][t1] - pos[id][t0])/r;
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
           	double r2 = r1;
           	double rd = r/BOHRRADIUS;
           	double potl;
           	vh2h2_(&rd, &r1, &r2, &th1, &th2, &phi, &potl);
           	spot += potl*CMRECIP2KL;
		}  //stype

#ifdef IOWRITE
     bool wline = true;                  // skip if the time slice between ira and masha

     if (WORM && Worm.exists && (Worm.type == type1))  
     wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
    
// PotEnergy (int, double) vs PotEnergy (int,int, double, int):
// the difference only one line below (commented out)
// this from PotEnergy (int, double) for (int it=0;it<NumbTimes;it++)
     if (wline)
     {  
        int t0 = offset0 + it;
        int t1 = offset1 + it;

        double dr2 = 0.0;  		 
        for (int id=0;id<NDIM;id++)
        {
           dr[id]  = (pos[id][t0] - MCCoords[id][t1]);

           if (MINIMAGE)
           dr[id] -= (BoxSize*rint(dr[id]/BoxSize));
 
           dr2    += (dr[id]*dr[id]);
        }
	       	 
//#ifdef _CUTOFF_	     
//     if (dr2<dljcutoff2)
//#endif
       double r = sqrt(dr2);

//-------------- MOLECULES ----------------------

       int tm;

       if ((MCAtom[type0].molecule == 1)||(MCAtom[type1].molecule == 1))  // 2D interaction 
       { 
          int sgn = 1;             // set to -1 to correct the orientaion of dr

          tm = offset1 + it/RotRatio;

          int typep = type1;         // define the model of interaction

          if (MCAtom[type0].molecule == 1)  // does not work for two molecules
          {
             sgn = -1;   
             tm  = offset0 + it/RotRatio;
             typep = type0; 
          }

          double cost = 0.0;
          for (int id=0;id<NDIM;id++) // n*dr = r*cos(theta) 
          cost += (MCCosine[id][tm]*dr[id]);   	 
	 
          cost /= r;                  // cos(theta)
          cost *= sgn;                // correct the orientation 

          spot += LPot2D(r,cost,typep);   
       }
//----------------ATOM-NON-LINEAR MOLECULES----------------
       else if ((((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && ISPHER == 0) &&(MCAtom[type0].molecule != MCAtom[type1].molecule) ) // 3D interaction
       {
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
           int    ivcord=0;
           if(MCAtom[type0].molecule == 2)
           {
              tm  = offset0 + it/RotRatio;
              for (int id=0;id<NDIM;id++)
              {
                 RCOM[id] = pos[id][t0];
                 Rpt[id]  = MCCoords[id][t1];
              }
           }
           else
           {
              tm  = offset1 + it/RotRatio;
              for (int id=0;id<NDIM;id++)
              {
                 Rpt[id]  = pos[id][t0];
                 RCOM[id] = MCCoords[id][t1];
              }
           }
           Eulang[PHI]=MCAngles[PHI][tm];
           Eulang[CTH]=acos(MCAngles[CTH][tm]);
           Eulang[CHI]=MCAngles[CHI][tm];

           vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

           spot += vpot3d;

       }
//----------------ATOM-NON-LINEAR MOLECULES spherical ----------------
       else if ((((MCAtom[type0].molecule == 2)||(MCAtom[type1].molecule == 2)) && ISPHER == 1)&& (MCAtom[type0].molecule != MCAtom[type1].molecule)) // spherical treatment for non-linear rotor
       {

           double radret,vpot3d;
           radret = r;
           vspher_(&radret,&vpot3d);

           spot += vpot3d;

       }
//----------- NonLinear ---- NonLinear----------------
       else if  ( ((MCAtom[type0].molecule == 2)&&(MCAtom[type1].molecule == 2)) && (MCAtom[IMTYPE].numb > 1) )
       {
        // GG:
      //    if (MCType[atom1] == IMTYPE)
        //     {
        //      int t0 = offset0 + it;
        //      int t1 = offset1 + it;
           double com_1[3];
           double com_2[3];
           double Eulang_1[3];
           double Eulang_2[3];
           double E_2H2O;
           int t0 = offset0 + it;
           int t1 = offset1 + it;
           for (int id=0;id<NDIM;id++)
           {
                com_1[id] = pos[id][t0];
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
           spot += E_2H2O;
           //   }
       }
       else
       spot += SPot1D(r,type1);       // 1D interaction

// it shoud be SPot1D(r,type0,type1) or  SPot1D(r,ind) with ind =type0*NumbTypes+type1
     	} // END sum over time slices 	   
#endif
	}   // END sum over atoms
   	return (spot);
}

#ifdef IOWRITE
//double PotRotEnergy(int atom0, double **cosine, int it)   
double PotRotEnergy(int atom0, double *Eulang0, int it)   
//  Orientational energy 
{
	int type0   =  MCType[atom0];
#ifdef IOWRITE
#ifdef DEBUG_PIMC
	const char *_proc_=__func__;         //  PotRotEnergy()

	if ((type0 != IMTYPE) || (MCAtom[type0].molecule == 0))
	nrerror(_proc_,"Use PotEnergy(int atom0, double **pos, int it)");

	if (MCAtom[type0].numb != 1)
	nrerror(_proc_,"Only one molecular impurity");
#endif
#endif

	double spot = 0.0;

#ifdef IOWRITE
	double dr[NDIM];

	int offset0 =  atom0*NumbTimes;

	for (int atom1=0;atom1<NumbAtoms;atom1++)
	if (atom1 != atom0)                    // skip "self-interaction"
	{	
		int offset1 = atom1*NumbTimes;
		int type1   = MCType[atom1];

#ifdef DEBUG_PIMC
		if ((MCAtom[type1].molecule == 1) || (MCAtom[type1].molecule == 2) )
		nrerror(_proc_,"More then one molecular impurity type");
#endif

		bool wline = true;                  // skip if the time slice between ira and masha

		if (WORM && Worm.exists && (Worm.type == type1))  
		wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
    
		if (wline)
		{  
			int t0 = offset0 + it;
			int t1 = offset1 + it;

			double dr2 = 0.0;  		 
			for (int id=0;id<NDIM;id++)
			{
				dr[id]  = (MCCoords[id][t0] - MCCoords[id][t1]);

				if (MINIMAGE)
				{
					cout << "MIN IMAGE for orient pot" << endl; exit(0);
					dr[id] -= (BoxSize*rint(dr[id]/BoxSize));
				}
 
				dr2    += (dr[id]*dr[id]);
			}
	       	 
//#ifdef _CUTOFF_	     
//     		if (dr2<dljcutoff2)
//#endif
			double r = sqrt(dr2);

			int sgn = -1;   
			int tm  = offset0 + it/RotRatio;
//   		int tm  = offset0 + (int)floor((double)it/(double)RotRatio);

			double cost = 0.0;
			for (int id=0;id<NDIM;id++) // n*dr = r*cos(theta) 
			cost += (cosine[id][tm]*dr[id]);   	 
	 
			cost /= r;                  // cos(theta)
			cost *= sgn;                // correct the orientation 
          	spot += (LPot2D(r,cost,type0));  
     	} 
	}   // END sum over atoms
#endif

	double spotSwap = 0.0;
    double weight, weight1;
	weight = 1.0;
#ifdef PIGSTYPE
    if (it == 0 || it == (NumbRotTimes - 1)) weight = 0.5;
#endif

	if ( (MCAtom[type0].molecule == 4) && (MCAtom[type0].numb > 1) )
	{
	    int offset0 =  atom0*NumbRotTimes;
        int t0  = offset0 + it;
		double cosine[NDIM][NumbAtoms*NumbRotTimes];
		cosine[0][t0] = sin(Eulang0[CTH])*cos(Eulang0[PHI]);
		cosine[1][t0] = sin(Eulang0[CTH])*sin(Eulang0[PHI]);
		cosine[2][t0] = cos(Eulang0[CTH]);

		int atom1Init, NumbAtoms1;
#ifdef ENTANGLEMENT
   	    int particleA1Min = (NumbAtoms/2) - NumbParticle;
   	    int particleA1Max = particleA1Min + NumbParticle - 1;
   	    int particleA2Min = particleA1Max + 1;
   	    int particleA2Max = particleA2Min + NumbParticle - 1;

        if (atom0 <= particleA1Max)
        {
            atom1Init  = 0;
            NumbAtoms1 = (particleA1Max+1);
        }
        else
        {
            atom1Init  = particleA2Min;
            NumbAtoms1 = NumbAtoms;
        }
#else
        atom1Init  = 0;
        NumbAtoms1 = NumbAtoms;
#endif
        spot = 0.0;
        for (int atom1 = atom1Init; atom1 < NumbAtoms1; atom1++)
        if (atom1 != atom0)                    
        {
            int offset1 = atom1*NumbTimes;
            int t1  = offset1 + it;

	        string stype = MCAtom[type0].type;
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
                    double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*cosine[id][t0];
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
                    b1[id] = cosine[id][t0];
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
			}  //stype

		    if (stype == HF )
            {
				weight1 = 1.0;
#ifndef REGULARPATH
#ifdef ENTANGLEMENT
				if (((atom0 < particleA1Min) || (atom0 > particleA2Max)) && ((atom1 >= particleA1Min) && (atom1 <= particleA2Max)))
            	{
           			if (it == ((NumbRotTimes - 1)/2))
              		{
                   		weight1 = 0.5;
					}
            	} 
            	if (((atom0 >= particleA1Min) && (atom0 <= particleA2Max)) && ((atom1 < particleA1Min) || (atom1 > particleA2Max)))
            	{
            		if (it == ((NumbRotTimes - 1)/2))
               		{
                   		weight1 = 0.5;
					}
            	} 
#endif
#endif
				double Eulang1[NDIM];
				Eulang1[PHI] = MCAngles[PHI][t1];
        		Eulang1[CTH] = acos(MCAngles[CTH][t1]);
        		Eulang1[CHI] = 0.0;
        		spot += weight*weight1*PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            }  //stype
        } //loop over atom1 (molecules)
#ifdef SWAP
		if (it == ((NumbRotTimes - 1)/2))
		{
			int atomSwapMin, atomSwapMax;
			if (atom0 < particleA1Min)
    	   	{
				atomSwapMin = particleA2Min;
                atomSwapMax = (particleA2Max + 1);
			}
			if (atom0 > particleA2Max)
    	   	{
				atomSwapMin = particleA1Min;
                atomSwapMax = (particleA1Max + 1);
			}
	    	if ((atom0 >= particleA1Min) && (atom0 <= particleA1Max))
       	    {
			    atomSwapMin = (particleA2Max + 1);
                atomSwapMax = NumbAtoms;
		    }
     		if ((atom0 >= particleA2Min) && (atom0 <= particleA2Max))
         	{
	     		atomSwapMin = 0;
                atomSwapMax = particleA1Min;
		    }
			spotSwap = 0.0;
     		for (int atomSwap = atomSwapMin; atomSwap < atomSwapMax; atomSwap++)
	    	{
                int offsetSwap = NumbRotTimes*atomSwap;
                int tSwap  = offsetSwap + it;
	        	string stype = MCAtom[type0].type;
                if (stype == HF )
                {
					double EulangSwap[NDIM];
					EulangSwap[PHI] = MCAngles[PHI][tSwap];
        			EulangSwap[CTH] = acos(MCAngles[CTH][tSwap]);
        			EulangSwap[CHI] = 0.0;
        			spotSwap += 0.5*PotFunc(atom0, atomSwap, Eulang0, EulangSwap, it);
                }  //stype
		    }
        }
#endif
    }

	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
	{
        double E12 = -2.0*DipoleMomentAU2*cos(Eulang0[CTH])/(RR*RR*RR);
        spot        = E12*AuToKelvin;
    }
    double spot_cage;
#ifdef CAGEPOT
    double cost = cos(Eulang0[CTH]);
    double phi = Eulang0[PHI];
    if (phi < 0.0) phi = 2.0*M_PI + phi;
    phi = fmod(phi,2.0*M_PI);
    spot_cage = weight*LPot2DRotDOF(cost,phi,type0);
#else
    spot_cage = 0.0;
#endif
	double spotReturn = (spot + spotSwap + spot_cage);
    return spotReturn;
}
#endif

double PotRotE3D(int atom0,double * Eulang,int it)   
//  Orientational energy 
{
   int type0   =  MCType[atom0];

#ifdef DEBUG_PIMC
   const char *_proc_=__func__;         //  PotRotEnergy()

   if ((type0 != IMTYPE) || (MCAtom[type0].molecule == 0))
   nrerror(_proc_,"Use PotEnergy(int atom0, double **pos, int it)");

   if (MCAtom[type0].numb > NumbRotLim)
   nrerror(_proc_,"Too many non-linear rotors");
#endif

   double spot = 0.0;

   int offset0 =  atom0*NumbTimes;

   for (int atom1=0;atom1<NumbAtoms;atom1++)
   if (atom1 != atom0)                    // skip "self-interaction"
   {	
      int offset1 = atom1*NumbTimes;
      int type1   = MCType[atom1];

#ifdef DEBUG_PIMC
//    if ((MCAtom[type1].molecule == 1) || (MCAtom[type1].molecule == 2) )
//    nrerror(_proc_,"More then one molecular impurity type");
      if(MCAtom[type1].molecule == 1)
      nrerror(_proc_,"No support of non-linear-linear interaction yet");
#endif

      if (type1 != IMTYPE) // atom-rotor interaction
      {
         bool wline = true;                  // skip if the time slice between ira and masha

         if (WORM && Worm.exists && (Worm.type == type1))  
         wline = WorldLine((atom1-MCAtom[type1].offset/NumbTimes), it);
    
         if (wline)
         {  
            int t0 = offset0 + it;
            int t1 = offset1 + it;

            double RCOM[3];
            double Rpt[3];
            double vpot3d;
            double radret;
            double theret;
            double chiret;
            double hatx[3];
            double haty[3];
            double hatz[3];
            int    ivcord = 0;
            for (int id=0;id<NDIM;id++)
            {
               RCOM[id] = MCCoords[id][t0];
               Rpt[id]  = MCCoords[id][t1];
            }

            vcord_(Eulang,RCOM,Rpt,vtable,&Rgrd,&THgrd,&CHgrd,&Rvmax,&Rvmin,&Rvstep,&vpot3d,&radret,&theret,&chiret,hatx,haty,hatz,&ivcord);

//       for(int id=0;id<NDIM;id++)
/*  Toby's printing
         cout<<Eulang[id]<<" "<<RCOM[id]<<" "<<Rpt[id]<<endl;
         cout<<vpot3d<<endl;
*/

            spot += vpot3d;
 
         } // END sum over time slices 	   
      }
      else if (MCType[atom1] == IMTYPE)
      {
         int t0 = offset0 + it;
         int t1 = offset1 + it;
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
                   Eulang, Eulang_2);
         spot += E_2H2O;
//       cout<<"in PotRotE3D "<<it<<" "<<t0<<" "<<t1<<" "<<tm0<<" "<<tm1<<" "<<" "<<offset0<<" "<<offset1<<" "<<E_2H2O<<endl;
//       cout<<Eulang[PHI]<<" "<<Eulang[CTH]<<" "<<Eulang[CHI]<<" "<<Eulang_2[PHI]<<" "<<Eulang_2[CTH]<<" "<<Eulang_2[CHI]<<endl;
//       cout<<com_1[0]<<" "<<com_1[1]<<" "<<com_1[2]<<" "<<com_2[0]<<" "<<com_2[1]<<" "<<com_2[2]<<endl;
      }

   }   // END sum over atoms
   return (spot);
}

void ResetMCCounts(void)
{
   for (int type=0;type<NumbTypes;type++)
   for (int im=0;im<MCMAXMOVES;im++)
   {
      MCTotal[type][im] = 0;
      MCAccep[type][im] = 0;
   }
}

void MemAllocMCCounts(void)
{
   MCTotal = doubleMatrix(NumbTypes,MCMAXMOVES);
   MCAccep = doubleMatrix(NumbTypes,MCMAXMOVES);
}

void MFreeMCCounts(void)
{
   free_doubleMatrix(MCTotal);
   free_doubleMatrix(MCAccep);
}

#ifdef IOWRITE
double MCQuaternions(double Aa, double Bb, double Cc, double Dd)
{
	double RMat[NDIM][NDIM];

	RMat[0][0] = Aa*Aa+Bb*Bb-Cc*Cc-Dd*Dd;
	RMat[0][1] = 2.0*Bb*Cc-2.0*Aa*Dd;
	RMat[0][2] = 2.0*Bb*Dd+2.0*Aa*Cc;

	RMat[1][0] = 2.0*Bb*Cc+2.0*Aa*Dd;
	RMat[1][1] = Aa*Aa-Bb*Bb+Cc*Cc-Dd*Dd;
	RMat[1][2] = 2.0*Cc*Dd-2.0*Aa*Bb;

	RMat[2][0] = 2.0*Bb*Dd-2.0*Aa*Cc;
	RMat[2][1] = 2.0*Cc*Dd+2.0*Aa*Bb;
	RMat[2][2] = Aa*Aa-Bb*Bb-Cc*Cc+Dd*Dd;
	
	return MCCosines;
}
#endif
#ifdef PROPOSED
int myRand(double *freq, double rand1)
{
    // Create and fill prefix array
    int nn = NCOST*NPHI;
    double prefix[nn];
    prefix[0] = freq[0];
    for (int i = 1; i < nn; ++i)
    {
        prefix[i] = prefix[i - 1] + freq[i];
    }

    // prefix[n-1] is sum of all frequencies. Generate a random number
    // with value from 1 to this sum

    // Find index of ceiling of r in prefix arrat
    int indexc = findCeil(prefix, rand1);
    return indexc;
}

// Utility function to find ceiling of r in arr[l..h]
int findCeil(double *arr, double rand1)
{
    int l = 0;
    int h = NCOST*NPHI-1;
    int mid;
    while (l < h)
    {
        mid = l + ((h - l) >> 1); // Same as mid = (l+h)/2
        (rand1 > arr[mid]) ? (l = mid + 1) : (h = mid);
    }
    return (arr[l] >= rand1) ? l : -1;
}
#endif

//double PotRotEnergy(int atom0, double **cosine, int it)   
double PotRotEnergyPIMC(int atom0, double *Eulang0, int it)   
{
	int type0   =  MCType[atom0];

	double spot;

	if ( (MCAtom[type0].molecule == 4) && (MCAtom[type0].numb > 1) )
	{
	    int offset0 =  atom0*NumbRotTimes;
        int t0  = offset0 + it;
		double cosine[NDIM][NumbAtoms*NumbRotTimes];
		cosine[0][t0] = sin(Eulang0[CTH])*cos(Eulang0[PHI]);
		cosine[1][t0] = sin(Eulang0[CTH])*sin(Eulang0[PHI]);
		cosine[2][t0] = cos(Eulang0[CTH]);

        spot = 0.0;
        for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
        if (atom1 != atom0)                    
        {
            int offset1 = atom1*NumbRotTimes;
            int t1  = offset1 + it;

	        string stype = MCAtom[type0].type;
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
                    double cst1 = (MCCoords[id][t1] - MCCoords[id][t0])*cosine[id][t0];
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
                    b1[id] = cosine[id][t0];
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
			}  //stype

		    if (stype == HF )
            {
				double Eulang1[NDIM];
				Eulang1[PHI] = MCAngles[PHI][t1];
        		Eulang1[CTH] = acos(MCAngles[CTH][t1]);
        		Eulang1[CHI] = 0.0;
        		spot += PotFunc(atom0, atom1, Eulang0, Eulang1, it);
            }  //stype
        } //loop over atom1 (molecules)
    }

	if ( (MCAtom[IMTYPE].molecule == 4) && (MCAtom[IMTYPE].numb == 1) )
	{
        double E12 = -2.0*DipoleMomentAU2*cos(Eulang0[CTH])/(RR*RR*RR);
        spot       = E12*AuToKelvin;
#ifdef POTZERO
		spot       = 0.0;
#endif
    }
    double spot_cage;
#ifdef CAGEPOT
    double cost = cos(Eulang0[CTH]);
    double phi = Eulang0[PHI];
    if (phi < 0.0) phi = 2.0*M_PI + phi;
    phi = fmod(phi,2.0*M_PI);
    spot_cage = weight*LPot2DRotDOF(cost,phi,type0);
#else
    spot_cage = 0.0;
#endif
	double spotReturn = (spot + spot_cage);
    return spotReturn;
}


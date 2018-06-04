//
//         main()
//

#include <string>
#include <stdlib.h>
//#include <stdlib>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <vector>

#include "mc_confg.h"
#include "mc_setup.h"
#include "mc_input.h"
#include "mc_randg.h"
#include "mc_utils.h"
#include "mc_poten.h"
#include "mc_piqmc.h"
#include "mc_qworm.h"
#include "mc_estim.h"
#include "mc_const.h"

#include <stdio.h>
//#include <mpi.h>
#include <omp.h>
#include "rngstream.h"
#include "omprng.h"


void MCWormAverage(void);
void MCWormAverageReset(void);

double avergCount;   // # of calls of get_estim inside a block
#ifdef PIGSENTBOTH
double avergCountENT;   // # of calls of get_estim inside a block
double totalCountENT;
#endif
double totalCount;   // sum avergCount   
double sumsCount;
double totalStep;

void PIMCPass(int,int);

void MCGetAverage(void);
void MCResetBlockAverage(void);
void MCGetAveragePIMC(void);
void MCResetBlockAveragePIMC(void);
void MCGetAveragePIGS(void);
void MCResetBlockAveragePIGS(void);
void MCGetAveragePIGSENT(void);
void MCResetBlockAveragePIGSENT(void);

void MCSaveBlockAverages(long int);
void MCSaveBlockAveragesPIGSENT(long int);

void MCSaveAcceptRatio(long int,long int,long int);

//--------------- BLOCK AVERAGE ------------

double _dbpot;       // potential energy differencies, block average  added by Hui Li
double _bpot;       // kinetic   energy, block average
double _btotal;
double _bkin;       // potential energy, block average
double _bCv;        // heat capacity, block average
double _bCv_trans;  // translational heat capacity, block average
double _bCv_rot;   //  rotational heat capacity, block average

double _dpot_total;  // potential energy differences, global average  added by Hui Li
double _pot_total;  // kinetic   energy, global average
double _kin_total;  // potential energy, global average 
double _total;  // potential energy, global average 
double _bcostheta;
double _ucompx;
double _ucompy;
double _ucompz;
double _costheta_total;
double _ucompx_total;
double _ucompy_total;
double _ucompz_total;
vector<double> _cdipoleXYZ;
vector<double> _cdipoleX;
vector<double> _cdipoleY;
vector<double> _cdipoleZ;
vector<double> _cdipoleXY;
vector<double> _cdipoleXYZ_total;
vector<double> _cdipoleX_total;
vector<double> _cdipoleY_total;
vector<double> _cdipoleZ_total;
vector<double> _cdipoleXY_total;
//
double _bnm;
double _bdm;
double _nm_total;
double _dm_total;
double _trOfDensitySq;
double _trOfDensitySq_total;

double _brot;       // rotational kin energy, block average
double _brot1;       // rotational kin energy, block average
double _rot_total;  // rotational kin energy, global average
double _rot_total1;  // rotational kin energy, global average
double _brotsq;     // rotational energy square, block average
double _rotsq_total; // rotational energy square, global average
double _Cv_total;    // heat capacity, global average
double _Cv_trans_total;    // translational heat capacity, global average
double _Cv_trans_1_total;    // translational heat capacity, global average
double _Cv_trans_2_total;    // translational heat capacity, global average
double _Cv_rot_total;    // rotational heat capacity, global average

fstream _feng;      // save accumulated energy
fstream _fang;      // save accumulated energy
fstream _fdc;      // save accumulated energy
fstream _fangins;      // save accumulated energy
fstream _fanginsc;      // save accumulated energy
fstream _fengins;      // save accumulated energy
fstream _fentropy;
fstream _fswap;

//---------------- ESTIMATORS ---------------

void SaveEnergy    (const char [],double,long int); // block average
void SaveSumEnergy (double,double);                 // accumulated average

void SaveInstantAngularDOF(long int); 
void SaveInstantEnergy ();                 // accumulated average

#ifdef DDCORR
void SaveDipoleCorr(const char [],double,long int);
void SaveSumDipoleCorr(double,double);              // accumulated average
#endif
void SaveAngularDOF(const char [],double,long int);
void SaveSumAngularDOF(double,double);              // accumulated average
void SaveTrReducedDens(const char [], double , long int );
void SaveSumTrReducedDens(double , double);

void InitTotalAverage(void);
void DoneTotalAverage(void);

//---------------- INITIAL MCCOORDS -----------

extern "C" void initconf_(double *MCCooInit, double *MCAngInit, int *Ntotal, int *NBoson, int *PIndex,int *Rindex);

//---------------- PRINT PERMUTATION TABLE -----

extern "C" void prtper_(int *PIndex,int *NBoson,long int *blockCount);

//----------- subroutine that performs some initialization work for vh2h2.f-----

extern "C" void vinit_();
string Distribution = "unSwap";
double MCAccepSwap;
double MCAccepUnSwap;
int iSwap;
int iUnSwap;
void SaveInstantSwap(double, double, long int);                 // accumulated average
#ifdef PROPOSED
int iChooseOld = 0;
int iChooseNew;
int iChoose;
#endif

//-------------------------------------------

int main(int argc, char *argv[])
{
#ifdef PAIRDENSITY
	readPairDensity();
	exit(11);
#endif
   	randomseed(); //set seed according to clock
//Tapas added
#ifdef SWAPTOUNSWAP
srand( time(NULL) );
#endif
#ifdef POTH2
   vinit_();
#endif

//int numprocs, rank, namelen;
//char processor_name[MPI_MAX_PROCESSOR_NAME];


//MPI_Init(&argc, &argv);
//MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//MPI_Get_processor_name(processor_name, &namelen);
//MPI_Status mpistatus;

 
//-------------------------
   MPIsize = 1;            // this is for InitRandom()          
   MPIrank = MPI_MASTER;   // this is for InitRandom()          
//-----------------------------

   int restart = 0;
   IOReadParams(FINPUT,restart); // read input parameters 
// get the number of rotational steps treated by worker CPUs
// chunksize = NumbRotTimes / numprocs;

// pass numprocs to NProcs, which is a global variable
// NProcs = numprocs;

   MCInitParams();               // set system dependent parameters
   MCSetUnits(); 
   MCMemAlloc();

   MemAllocMCCounts();
   MemAllocQWCounts();

//head CPU comes in after memory declaration.
//if ( rank == MPI_MASTER ){// in the master cpu section

// only internal system of units after this point 
 
	MCInit();
#ifdef PROPOSED
proposedGrid();
#endif
ParamsPotential();
#ifdef GAUSSIANMOVE
	ProposedMCCoords();
#endif

	if (WORM)
	MCWormInit();

// clean the permutation sampling from the last run
   	if (FileExist(FPERMU))
   	_io_error("QMC ->",IO_ERR_FEXST,FPERMU);

   	if (!restart) // new run, generate new status, rnd() streams and config files     
   	{
// check the existence of the old files first
  
      	if (FileExist(FSTATUS))    
     	_io_error("QMC ->",IO_ERR_FEXST,FSTATUS);
       
      	if (FileExist(FCONFIG)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FCONFIG);

      	if (FileExist(FRANDOM)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FRANDOM);

      	if (FileExist(FSEED)) 
     	_io_error("QMC ->",IO_ERR_FEXST,FSEED);

//--------------------------------------------------------
//    	generate tables - potentals, configurations etc, status

//    	MCStartBlock = 0; 
//    	SEED for head CPU
      	SEED = 985456376;
      	RandomInit(MPIrank,MPIsize);

    	MCConfigInit();                // generate initial configurations
    	for (int it=0;it<NumbAtoms*NumbTimes;it++)
    	{
        	for (int id=0;id<NDIM;id++)
        	{
            	cout<<"it " << it<<" id "<< id<< " "<< MCCoords[id][it]<<endl;
        	}
    	} 
    	cout<<"  "<<endl;
    	cout<<"  "<<endl;
		for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
		{
        	int offset0 = NumbTimes*atom0;
			for (int atom1 = 0; atom1 < NumbAtoms; atom1++)
			{
				if (atom0 != atom1)
				{
            		int offset1 = NumbTimes*atom1;

            		int it = ((NumbRotTimes - 1)/2);
            		int t0 = offset0 + it;
            		int t1 = offset1 + it;

            		double dr;
            		double dr2 = 0.0;
            		for (int id=0;id<NDIM;id++)
            		{
                		dr   = (MCCoords[id][t0] - MCCoords[id][t1]);
                		dr2 += (dr*dr);
            		}
            		double r = sqrt(dr2);
					cout<<"atom0 "<<atom0<<" atom1 "<<atom1<<" RCOM "<<r<<endl;
				}
			}
			cout<<"  "<<endl;
		}
    	cout<<"  "<<endl;
    	cout<<"  "<<endl;
    	for (int it=0;it<NumbAtoms*NumbTimes;it++)
    	{
        	for (int id=0;id<NDIM;id++)
        	{
            	cout<<"it " << it<<" id "<< id<< " "<< MCCosine[id][it]<<endl;
        	}
    	}	 
//    	read in initial MCCoords and MCAngles
      	if(InitMCCoords)
      	{
         	cout<<"read in MCCoords and MCAngles from xyz.init file"<<endl;
         	int ntotal=NumbAtoms*NumbTimes;
         	int iperm,nboson;
         	nboson=0;
         	if(BOSONS)
         	{
            	nboson=MCAtom[BSTYPE].numb;
         	}
//       	initconf_(MCCooInit,MCAngInit,&ntotal,&MCAtom[BSTYPE].numb,PIndex,RIndex);
         	initconf_(MCCooInit,MCAngInit,&ntotal,&nboson,PIndex,RIndex);
         	if(BOSONS)
         	{
            	cout<<"read PIndex: "<<" ";
            	for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            	cout<<PIndex[id]<<" ";

            	cout<<endl;

            	cout<<"read RIndex: "<<" ";
            	for(int id=0;id<MCAtom[BSTYPE].numb;id++)
            	cout<<RIndex[id]<<" ";

            	cout<<endl;
         	}

         	for (int it=0;it<NumbAtoms*NumbTimes;it++)
         	{
            	for(int id=0;id<NDIM;id++)
            	{
               		MCCoords[id][it]=MCCooInit[it*NDIM+id];
               		MCAngles[id][it]=MCAngInit[it*NDIM+id];
            	}
         	}
      	}
		delete [] MCCooInit;
		delete [] MCAngInit;

    	//string fname = MCFileName + IO_EXT_XYZ;
    	//IOxyz(IOWrite,fname.c_str());  // test output of initial config
      	string fname = MCFileName;
      	IOxyzAng(IOWrite,fname.c_str()); // test output of initial config

//--------------------------------------------------------

      	StatusIO(IOWrite,FSTATUS);
      	ConfigIO(IOWrite,FCONFIG);
      	TablesIO(IOWrite,FTABLES);
      	RandomIO(IOWrite,FRANDOM);
		SeedIO(IOWrite,FSEED);

      	if (WORM)
      	QWormsIO(IOWrite,FQWORMS);
   	}

   	MCWormAverageReset();      // debug worm

// --- RESTART/START NEW RUN ----------------------------

   	if (restart) // new run, generate new status, rnd() streams and config files     
   	{
		StatusIO(IORead,FSTATUS);  // load MCStartBlock
   		ConfigIO(IORead,FCONFIG);  // load atoms/molecules positions
   		TablesIO(IORead,FTABLES);  // load permutation tables
   		RandomIO(IORead,FRANDOM);  // load rnd streams 
		SeedIO(IORead,FSEED);
      	string fname = MCFileName;
      	IOxyzAng(IOWrite,fname.c_str()); // test output of initial config
	}

    for (int it=0;it<NumbAtoms*NumbTimes;it++)
    {

        double phi  = MCAngles[PHI][it];
        double cost = MCAngles[CTH][it];
        double sint = sqrt(1.0 - cost*cost);
        double chi  = MCAngles[CHI][it];

        MCCosine[AXIS_X][it] = sint*cos(phi);
        MCCosine[AXIS_Y][it] = sint*sin(phi);
        MCCosine[AXIS_Z][it] = cost;
    }

   	if (WORM)
   	QWormsIO(IORead,FQWORMS);

// BEGIN loop over blocks ------------------------

   	InitPotentials();

   	if (ROTATION)
    	InitRotDensity();

   	if((IREFLX == 1 || IREFLY == 1 || IREFLZ == 1) && MCAtom[IMTYPE].molecule != 2)
   	{
      	cout<<IREFLX<<" "<<IREFLY<<" "<<IREFLZ<<" "<<MCAtom[IMTYPE].molecule<<endl;
      	nrerror("main","SYMMETRIZATION IS ONLY IMPLEMENTED FOR TOPS");
   	}
// MPI BLOCK 1 in MASTER start
// send run information to all the worker CPUs
/*
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&NumbTypes,1,MPI_INT,iproc,tagNumbTypes,MPI_COMM_WORLD);
      MPI_Send(&NumbAtoms,1,MPI_INT,iproc,tagNumbAtoms,MPI_COMM_WORLD);
      MPI_Send(&rhoprp[0],SizeRotDen,MPI_DOUBLE,iproc,tagrho,MPI_COMM_WORLD);
      MPI_Send(&erotpr[0],SizeRotDen,MPI_DOUBLE,iproc,tagero,MPI_COMM_WORLD);
      MPI_Send(&vtable[0],SizePotTab,MPI_DOUBLE,iproc,tagpot,MPI_COMM_WORLD);
      MPI_Send(&WORM,1,MPI_INT,iproc,tagWORM,MPI_COMM_WORLD);
      MPI_Send(&MCType[0],NumbAtoms,MPI_INT,iproc,tagMCType,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].molecule,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_molecule,MPI_COMM_WORLD);
      MPI_Send(&RotRatio,1,MPI_INT,iproc,tagRotRatio,MPI_COMM_WORLD);
      MPI_Send(&Worm.type,1,MPI_INT,iproc,tagWormtype,MPI_COMM_WORLD);
      MPI_Send(&NumbRotTimes,1,MPI_INT,iproc,tagNumbRotTimes,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].offset,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_offset,MPI_COMM_WORLD);
      MPI_Send(&IMTYPE,1,MPI_INT,iproc,tagIMTYPE,MPI_COMM_WORLD);
      MPI_Send(&MCAtom[0].numb,MAX_NUMBER_TYPES,MPI_INT,iproc,tagMCAtom_numb,MPI_COMM_WORLD);
      MPI_Send(&MCTau,1,MPI_DOUBLE,iproc,tagMCTau,MPI_COMM_WORLD);
      MPI_Send(&ROTATION,1,MPI_INT,iproc,tagROTATION,MPI_COMM_WORLD);
      MPI_Send(&MCStartBlock,1,MPI_LONG,iproc,tagMCStartBlock,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCBlocks,1,MPI_LONG,iproc,tagNumberOfMCBlocks,MPI_COMM_WORLD);
      MPI_Send(&NumberOfMCPasses,1,MPI_LONG,iproc,tagNumberOfMCPasses,MPI_COMM_WORLD);
//    MPI_Send(&TZMAT[1][1],NDIM,MPI_DOUBLE,iproc,tagTZMAT,MPI_COMM_WORLD);
      MPI_Send(&RotDenType,1,MPI_INT,iproc,tagRotDenType,MPI_COMM_WORLD);
      MPI_Send(&RotOdEvn,1,MPI_INT,iproc,tagRotOdEvn,MPI_COMM_WORLD);
      MPI_Send(&RotEoff,1,MPI_DOUBLE,iproc,tagRotEoff,MPI_COMM_WORLD);
      MPI_Send(&X_Rot,1,MPI_DOUBLE,iproc,tagARot,MPI_COMM_WORLD);
      MPI_Send(&Y_Rot,1,MPI_DOUBLE,iproc,tagBRot,MPI_COMM_WORLD);
      MPI_Send(&Z_Rot,1,MPI_DOUBLE,iproc,tagCRot,MPI_COMM_WORLD);
      MPI_Send(&MCRotTau,1,MPI_DOUBLE,iproc,tagMCRotTau,MPI_COMM_WORLD);
   }
// tagrunning = (tagNumberOfMCPasses + 2) % tagUpper;
// tagrunning = (tagTZMAT + 2) % tagUpper;
   tagrunning = (tagMCRotTau + 2) % tagUpper;

// the following loop passes the same tagrunning to all workers, in order to guarantee all CPUs have the same running tag
   for(int iproc=1;iproc<numprocs;iproc++)
   {
      MPI_Send(&tagrunning,1,MPI_INT,iproc,tagTAGRUN,MPI_COMM_WORLD);
   }
// MPI BLOCK 1 in MASTER done
*/


    InitMCEstims();
   	InitTotalAverage();      // DUMP 
   
   	ResetMCCounts();
   	ResetQWCounts();

/*
// 	try openmp for loop
   	int chunksize,nthrds;
   	#pragma omp parallel
   	{
   	int tid=omp_get_thread_num();
   	nthrds=omp_get_num_threads();
   	chunksize=NumbRotTimes/nthrds;
   	cout<<"chunksize="<<chunksize<<endl;
   int itini=chunksize*tid;
   int itfnl=itini+chunksize;
// #pragma omp for
   for (int itrot=itini;itrot<itfnl-1;itrot++)
   {
      cout<<"itrot="<<itrot<<" tid="<<tid<<endl;
   }
   }

   for (int itrot=chunksize-1;itrot<NumbRotTimes;itrot=itrot+chunksize)
   {
      cout<<"itrot="<<itrot<<endl;
   }
   for (int itrot=nthrds*chunksize;itrot<NumbRotTimes;itrot++)
   {
      cout<<"itrot="<<itrot<<endl;
   }
*/

// get the chunksize and # of threads for the parallel rotation
	#pragma omp parallel
	{
    	NThreads=omp_get_num_threads();
      	chunksize=NumbRotTimes/NThreads;
   	} // end omp

   	cout<<"NThreads="<<NThreads<<" chunksize="<<chunksize<<endl;

   	printf("OpenMP version: %d\n", _OPENMP);
   	//randomseed(); //set seed according to clock
	//RngStream Rng[omp_get_num_procs()];     // initialize a parallel RNG named "Rng"

   	long int blockCount = MCStartBlock;  
   	while (blockCount<NumberOfMCBlocks) // START NEW BLOCK      
   	{      
    	blockCount++; 
#ifdef PIMCTYPE
       	MCResetBlockAveragePIMC();
#endif
//
#ifdef PIGSTYPE
       	MCResetBlockAveragePIGS();
#endif
//
#ifdef PIGSENTTYPE
       	MCResetBlockAveragePIGSENT();
#endif
     
       	long int passCount = 0;        // BEGIN NEW MC PASS
       	long int passTotal = 0;        // total number of Moves = passCount*time 

       	//for (int time=0; time<1; time++)
       	while (passCount++ < NumberOfMCPasses) 
       	for (int time=0; time<NumbTimes; time++)
       	{
        	passTotal++;  

            //toby's printing
            //cout<<blockCount<<" "<<passCount<<" "<<time<<endl;
           	for (int type=0;type<NumbTypes;type++)
           	if (WORM && (type == Worm.type)) 
           	{

           		MCWormMove();

           		if (!Worm.exists)
           		{
					/* reactive */
               		int rt = nrnd2(NumbTimes);   // select a time slice

               		if ((type == BSTYPE) || (type == FERMTYPE))
               			MCBisectionMoveExchange (type,rt);
                 	else
                 		MCBisectionMove (type,rt);
						/* reactive */
#ifdef MOVECENTROIDTEST
               		if  (time == 0)
           			if ((type == BSTYPE) || (type == FERMTYPE))
           				MCMolecularMoveExchange(type);        
               		else
               			MCMolecularMove(type);        
#endif
           		}
        	} 
       		else
			{
		    	PIMCPass(type,time);
			}
#ifdef INSTANT
			totalStep++;
#ifdef IOWRITE
            SaveInstantEnergy (); 
#endif
			if (blockCount > (NumberOfMCBlocks - 10))
			{
		    	SaveInstantAngularDOF(totalStep);
#ifdef SWAPTOUNSWAP
	    		SaveInstantSwap(iUnSwap, iSwap, totalStep);
#endif
			}
#endif

   			if (blockCount>NumberOfEQBlocks)        // skip equilibration steps
   			{
				// evaluate averages 
           		if (passTotal % MCSKIP_AVERG == 0)   // skip correlated configurations
           		{
               		if (WORM)
               		{
                		if (!Worm.exists)
                		{
                			//MCGetAverage();
	                        //prtper_(PIndex,&MCAtom[BSTYPE].numb,&blockCount);

                   			// print the instantaneous xyz and prl info for the closed path
                   			if(PrintXYZprl)
                   			{
                       			stringstream bc;                // convert block # to string
                       			bc.width(IO_BLOCKNUMB_WIDTH);
                       			bc.fill('0');
                       			bc<<blockCount;
                       			string fname = MCFileName + bc.str();  // file name prefix including block #
                       			IOxyzAng(IOWrite,fname.c_str());
                       			PrintXYZprl = 0;
                   			}

                   		}
                   		else  
                		MCWormAverage();  /// TEST ONLY
                	}	   
                	else
                 	{
#ifdef IOWRITE
                    	MCGetAverage();
#endif
//
#ifdef PIMCTYPE
                    	MCGetAveragePIMC();
#endif
//
#ifdef PIGSTYPE
                    	MCGetAveragePIGS();
#endif
//
#ifdef PIGSENTTYPE
                    	MCGetAveragePIGSENT();
#endif
//
					    //omp_set_num_threads(1);
					    //MCGetAverage();
					    //omp_set_num_threads(NThreads);
                        //print the instantaneous xyz and prl info for the closed path
                    	if(PrintXYZprl)
                    	{
#ifdef IOWRITE
                    		stringstream bc;                // convert block # to string
                    		bc.width(IO_BLOCKNUMB_WIDTH);
                    		bc.fill('0');
                    		bc<<blockCount;
                    		string fname = MCFileName + bc.str();  // file name prefix including block #
                    		IOxyzAng(IOWrite,fname.c_str());
#endif
                    		PrintXYZprl = 0;
                    	}
                 	}
              	}

// DUMP, save global average. if avergCount = 0, then MCGetAverage is never called in this block.  All Saving steps are skipped.

              	if (passTotal % MCSKIP_TOTAL == 0 && avergCount)  
              	{
               		sumsCount += 1.0;                 
#ifdef PIGSENTTYPE
               		SaveSumEnergy (totalCountENT,sumsCount);
                    SaveSumTrReducedDens(totalCount, sumsCount);
#else
               		SaveSumEnergy (totalCount,sumsCount);
#endif
            	}
			}  
          
        	if (passTotal % MCSKIP_RATIO == 0)
        	MCSaveAcceptRatio(passTotal,passCount,blockCount);
		}                               
       
		// END  loop over MC passes (time slices) -----------------

		if (blockCount>NumberOfEQBlocks && avergCount)   // skip equilibration steps
		{
#ifdef PIGSENTTYPE
   			MCSaveBlockAveragesPIGSENT(blockCount);
#else
   			MCSaveBlockAverages(blockCount);
#endif
			// save accumulated interatomic distribution

#ifdef HISTOGRAM
			if (blockCount > 1500)
			SaveGxyzSum(MCFileName.c_str(),totalCount);
#endif
#ifdef NEWDENSITY
			if (blockCount > 600)
			SaveGraSum(MCFileName.c_str(),totalCount);
#endif

#ifdef IOWRITE
			if(IMPURITY && MCAtom[IMTYPE].molecule == 1)
	    	SaveDensities2D(MCFileName.c_str(),totalCount,MC_TOTAL);

	   		if(IMPURITY && MCAtom[IMTYPE].molecule == 3)
	     	SaveDensities2D(MCFileName.c_str(),totalCount,MC_TOTAL);

           	if(IMPURITY && MCAtom[IMTYPE].molecule == 2)
           	{
            	SaveDensities3D(MCFileName.c_str(),totalCount,MC_TOTAL); // this step takes lots of space.  temporarily turnned off
              	SaveRho1D(MCFileName.c_str(),totalCount,MC_TOTAL);
				// SaveRhoThetaChi(MCFileName.c_str(),totalCount,MC_TOTAL); // we don't need 2d angular distribution comparison now
           	}

#endif
#ifdef PIMCTYPE
			if (ROTATION)                  // DUMP  accumulated average
				SaveRCF(MCFileName.c_str(),totalCount,MC_TOTAL);
#endif
		}

		//  CHECKPOINT: save status, rnd streams and configs ------

		MCStartBlock = blockCount; 

		if (blockCount%1000 == 0)
		{
			IOFileBackUp(FSTATUS); StatusIO(IOWrite,FSTATUS);
			IOFileBackUp(FCONFIG); ConfigIO(IOWrite,FCONFIG);
			IOFileBackUp(FTABLES); TablesIO(IOWrite,FTABLES);
			IOFileBackUp(FRANDOM); RandomIO(IOWrite,FRANDOM);      
			IOFileBackUp(FSEED); SeedIO(IOWrite,FSEED);
		}
      
		string fname = MCFileName;// + IO_EXT_XYZ;
		IOxyzAng(IOWrite,fname.c_str());  // test output of initial config

       	if (WORM)
       	{
        	IOFileBackUp(FQWORMS);
        	QWormsIO(IOWrite,FQWORMS);
    	}
	}  
	// END  loop over blocks ------------------------
	// } // master cpu if block ends

	//MPI_Finalize();

	DoneTotalAverage();  // DUMP

	if (WORM)
		MCWormDone();

	DoneMCEstims();
	DonePotentials();

	if (ROTATION)
		DoneRotDensity();

	MCMemFree();
	RandomFree();

	MFreeMCCounts();
	MFreeQWCounts();

	return 1; 
}

void PIMCPass(int type,int time)
{
  // skip solvent and translation moves for rotations only
#ifdef MOVECOM
#ifdef GAUSSIANMOVE
	MCMolecularMoveGauss(type);
#else
   if (time == 0)
   MCMolecularMove(type);        
// move the solvent particles
   MCBisectionMove(type,time);
#endif
#endif

	if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 1)  // rotational degrees of freedom
    	MCRotationsMove(type);
   	if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 3)  // rotational degrees of freedom
    	MCRotationsMove(type);
	if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 2)  // non-linear rotor rotation added by Toby
    	MCRotations3D(type);
	if ((type == IMTYPE) && ROTATION && MCAtom[type].molecule == 4)  // linear rotor rotation added by Tapas Sahoo
    	MCRotationsMove(type);
}

#ifdef PIMCTYPE
void MCResetBlockAveragePIMC(void) 
{
	avergCount = 0.0;

	ResetMCEstims();

	ResetMCCounts();
	ResetQWCounts();

	_dbpot     = 0.0;  //added by Hui Li
	_bpot      = 0.0;
	_btotal    = 0.0;
	_bkin        = 0.0;

	_brot        = 0.0;
	_brotsq      = 0.0;
	_bCv         = 0.0;
	_bCv_trans   = 0.0;
	_bCv_rot     = 0.0;

	PrintXYZprl = 1;

	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;
#ifdef DDCORR
	int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
    _cdipoleXYZ.resize(NDIMDP);
    _cdipoleX.resize(NDIMDP);
    _cdipoleY.resize(NDIMDP);
    _cdipoleZ.resize(NDIMDP);
  	_cdipoleXY.resize(NDIMDP);
	
	for (int idp = 0; idp < NDIMDP; idp++)
	{
        _cdipoleXYZ[idp] = 0.0;
        _cdipoleX[idp] = 0.0;
        _cdipoleY[idp] = 0.0;
        _cdipoleZ[idp] = 0.0;
    	_cdipoleXY[idp] = 0.0;
	}
#endif
	_bcostheta = 0.0;
	_ucompx     = 0.0;
	_ucompy     = 0.0;
	_ucompz     = 0.0;
}
#endif

#ifdef PIGSTYPE
void MCResetBlockAveragePIGS(void) 
{
	avergCount = 0.0;

	ResetMCEstims();

	ResetMCCounts();
	ResetQWCounts();

	_dbpot     = 0.0;  //added by Hui Li
	_bpot      = 0.0;
	_btotal    = 0.0;
	_bcostheta = 0.0;
	_ucompx     = 0.0;
	_ucompy     = 0.0;
	_ucompz     = 0.0;
#ifdef DDCORR
	int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
    _cdipoleXYZ.resize(NDIMDP);
    _cdipoleX.resize(NDIMDP);
    _cdipoleY.resize(NDIMDP);
    _cdipoleZ.resize(NDIMDP);
  	_cdipoleXY.resize(NDIMDP);
	
	for (int idp = 0; idp < NDIMDP; idp++)
	{
        _cdipoleXYZ[idp] = 0.0;
        _cdipoleX[idp] = 0.0;
        _cdipoleY[idp] = 0.0;
        _cdipoleZ[idp] = 0.0;
    	_cdipoleXY[idp] = 0.0;
	}
#endif
	_bkin        = 0.0;

	_brot        = 0.0;
	_brot1        = 0.0;
	_brotsq      = 0.0;
	_bCv         = 0.0;
	_bCv_trans   = 0.0;
	_bCv_rot     = 0.0;

	PrintXYZprl = 1;

	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;

}
#endif

#ifdef PIGSENTTYPE
void MCResetBlockAveragePIGSENT(void) 
{
#ifdef SWAPTOUNSWAP
    MCAccepSwap = 0.0;
    MCAccepUnSwap = 0.0;
#endif
	avergCount = 0.0;

	ResetMCEstims();
	ResetMCCounts();
	ResetQWCounts();

    _bnm       = 0.0;
    _bdm       = 0.0;
	_trOfDensitySq = 0.0;

#ifdef PIGSENTBOTH
	avergCountENT = 0.0;
    _dbpot     = 0.0;  //added by Hui Li
    _bpot      = 0.0;
    _btotal    = 0.0;
    _bcostheta = 0.0;
    _ucompx     = 0.0;
    _ucompy     = 0.0;
    _ucompz     = 0.0;
#ifdef DDCORR
	int NumbAtoms1 = NumbAtoms/2;
	int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;
    _cdipoleXYZ.resize(NDIMDP);
    _cdipoleX.resize(NDIMDP);
    _cdipoleY.resize(NDIMDP);
    _cdipoleZ.resize(NDIMDP);
    _cdipoleXY.resize(NDIMDP);

    for (int idp = 0; idp < NDIMDP; idp++)
    {
        _cdipoleXYZ[idp] = 0.0;
        _cdipoleX[idp] = 0.0;
        _cdipoleY[idp] = 0.0;
        _cdipoleZ[idp] = 0.0;
        _cdipoleXY[idp] = 0.0;
    }
#endif
    _bkin        = 0.0;

    _brot        = 0.0;
    _brot1        = 0.0;
#endif
	PrintXYZprl = 1;
	PrintYrfl   = 1;
	PrintXrfl   = 1;
	PrintZrfl   = 1;
}
#endif

#ifdef PIGSENTTYPE
void MCGetAveragePIGSENT(void) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  

#ifdef SWAPTOUNSWAP
    double snm        = MCAccepSwap/(MCAccepSwap+MCAccepUnSwap);
    double sdm        = MCAccepUnSwap/(MCAccepSwap+MCAccepUnSwap);
#else
    double snm        = GetEstimNM();
    double sdm        = GetEstimDM();
#endif
    _bnm             += snm;
    _bdm             += sdm;
    _nm_total        += snm;
    _dm_total        += sdm;
#ifdef SWAP
	_trOfDensitySq   += sdm/snm;
	_trOfDensitySq_total += sdm/snm;
#endif
#ifdef UNSWAP
	_trOfDensitySq   += snm/sdm;
	_trOfDensitySq_total += snm/sdm;
#endif
#ifdef PIGSENTBOTH
	if (Distribution == "unSwap")
	{
		avergCountENT    += 1.0;
		totalCountENT    += 1.0;  

		double spot       = GetPotEnergyPIGSENT();//GetPotEnergyPIGS(); // pot energy and density distributions
		_bpot            += spot;                     // block average for pot energy
		_pot_total       += spot;
		double stotal     = GetTotalEnergyPIGSENT();//GetTotalEnergy();         // Total energy
		_btotal          += stotal;                   // kin+pot
		_total           += stotal;
		double srot1      = (stotal - spot);
		_brot1           += srot1;
		_rot_total1      += srot1;

		double cosTheta   = 0.0;
		double compxyz[NDIM];
		compxyz[0] = 0.0;
		compxyz[1] = 0.0;
		compxyz[2] = 0.0;

		GetCosThetaPIGSENT(cosTheta, compxyz);
		double scostheta  = cosTheta;
		double scompx     = compxyz[0];
		double scompy     = compxyz[1];
		double scompz     = compxyz[2];

		_bcostheta       += scostheta; 
		_costheta_total  += scostheta;

		_ucompx          += scompx;
		_ucompy          += scompy;
		_ucompz          += scompz;
		_ucompx_total    += scompx;
		_ucompy_total    += scompy;
		_ucompz_total    += scompz;
#ifdef DDCORR
		int NumbAtoms1 = NumbAtoms/2;
		int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;
    	double DipoleCorrXYZ[NDIMDP];
    	double DipoleCorrX[NDIMDP];
    	double DipoleCorrY[NDIMDP];
    	double DipoleCorrZ[NDIMDP];
    	double DipoleCorrXY[NDIMDP];

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			DipoleCorrXYZ[idp] = 0.0;
			DipoleCorrX[idp]   = 0.0;
			DipoleCorrY[idp]   = 0.0;
			DipoleCorrZ[idp]   = 0.0;
			DipoleCorrXY[idp]  = 0.0;
		}

		GetDipoleCorrelationPIGSENT(DipoleCorrXYZ, DipoleCorrX, DipoleCorrY, DipoleCorrZ, DipoleCorrXY);

		for (int idp = 0; idp < NDIMDP; idp++)
		{
			_cdipoleXYZ[idp] += DipoleCorrXYZ[idp];
			_cdipoleX[idp]   += DipoleCorrX[idp];
			_cdipoleY[idp]   += DipoleCorrY[idp];
			_cdipoleZ[idp]   += DipoleCorrZ[idp];
			_cdipoleXY[idp]  += DipoleCorrXY[idp];
			_cdipoleXYZ_total[idp] += DipoleCorrXYZ[idp];
			_cdipoleX_total[idp]   += DipoleCorrX[idp];
			_cdipoleY_total[idp]   += DipoleCorrY[idp];
			_cdipoleZ_total[idp]   += DipoleCorrZ[idp];
			_cdipoleXY_total[idp]  += DipoleCorrXY[idp];
		}
#endif
		double srot;
		if (ROTATION)
		{
			if(MCAtom[IMTYPE].molecule == 1)
			{
				srot      = GetRotEnergy();           // kin energy
			}

			if(MCAtom[IMTYPE].molecule == 3)
			{
				srot      = GetRotPlanarEnergy();     // kin energy
			}

			if(MCAtom[IMTYPE].molecule == 3)
			{
				srot          = GetRotE3D();
			}
        
			if(MCAtom[IMTYPE].molecule == 4)
			{
				srot      = GetRotEnergyPIGS();           // kin energy
			}
			_brot        += srot;
			_rot_total   += srot;
		}
	}
#endif
//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}
#endif

#ifdef PIGSTYPE
void MCGetAveragePIGS(void) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  

	double spot       = GetPotEnergyPIGS(); // pot energy and density distributions
	_bpot            += spot;                     // block average for pot energy
	_pot_total       += spot;
	double stotal     = GetTotalEnergy();         // Total energy
	_btotal          += stotal;                   // kin+pot
	_total           += stotal;
	double srot1      = (stotal - spot);
	_brot1           += srot1;
	_rot_total1      += srot1;

	double cosTheta   = 0.0;
	double compxyz[NDIM];
	compxyz[0] = 0.0;
	compxyz[1] = 0.0;
	compxyz[2] = 0.0;

	GetCosThetaPIGS(cosTheta, compxyz);
	double scostheta  = cosTheta;
	double scompx     = compxyz[0];
	double scompy     = compxyz[1];
	double scompz     = compxyz[2];

	_bcostheta       += scostheta; 
	_costheta_total  += scostheta;

	_ucompx          += scompx;
	_ucompy          += scompy;
	_ucompz          += scompz;
	_ucompx_total    += scompx;
	_ucompy_total    += scompy;
	_ucompz_total    += scompz;

#ifdef DDCORR
	int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
    double DipoleCorrXYZ[NDIMDP];
    double DipoleCorrX[NDIMDP];
    double DipoleCorrY[NDIMDP];
    double DipoleCorrZ[NDIMDP];
    double DipoleCorrXY[NDIMDP];

	for (int idp = 0; idp < NDIMDP; idp++)
	{
		DipoleCorrXYZ[idp] = 0.0;
		DipoleCorrX[idp]   = 0.0;
		DipoleCorrY[idp]   = 0.0;
		DipoleCorrZ[idp]   = 0.0;
		DipoleCorrXY[idp]  = 0.0;
	}

	GetDipoleCorrelationPIGS(DipoleCorrXYZ, DipoleCorrX, DipoleCorrY, DipoleCorrZ, DipoleCorrXY);

	for (int idp = 0; idp < NDIMDP; idp++)
	{
		_cdipoleXYZ[idp] += DipoleCorrXYZ[idp];
		_cdipoleX[idp]   += DipoleCorrX[idp];
		_cdipoleY[idp]   += DipoleCorrY[idp];
		_cdipoleZ[idp]   += DipoleCorrZ[idp];
		_cdipoleXY[idp]  += DipoleCorrXY[idp];
		_cdipoleXYZ_total[idp] += DipoleCorrXYZ[idp];
		_cdipoleX_total[idp]   += DipoleCorrX[idp];
		_cdipoleY_total[idp]   += DipoleCorrY[idp];
		_cdipoleZ_total[idp]   += DipoleCorrZ[idp];
		_cdipoleXY_total[idp]  += DipoleCorrXY[idp];
	}
#endif
	double srot;
	if (ROTATION)
	{
		if(MCAtom[IMTYPE].molecule == 1)
		{
			srot      = GetRotEnergy();           // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot          = GetRotE3D();
		}
        
		if(MCAtom[IMTYPE].molecule == 4)
		{
			srot      = GetRotEnergyPIGS();           // kin energy
		}
		_brot        += srot;
		_rot_total   += srot;
	}

//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}
#endif

#ifdef PIMCTYPE
void MCGetAveragePIMC(void) 
{
	avergCount       += 1.0;
	totalCount       += 1.0;  
#ifdef HISTOGRAM
	GetDensities();
#endif

    double skin       = 0.;
#ifdef MOVECOM
	skin              = GetKinEnergy();           // kin energy
#endif
	_bkin            += skin;                     // block average for kin energy
	_kin_total       += skin;                     // accumulated average 

	double spot       = GetPotEnergy_Densities(); // pot energy and density distributions
	_bpot            += spot;                     // block average for pot energy
	_pot_total       += spot;


    double cosTheta   = 0.0;
    double compxyz[NDIM];
    compxyz[0] = 0.0;
    compxyz[1] = 0.0;
    compxyz[2] = 0.0;

    GetCosThetaPIMC(cosTheta, compxyz);
    double scostheta  = cosTheta;
    double scompx     = compxyz[0];
    double scompy     = compxyz[1];
    double scompz     = compxyz[2];

    _bcostheta       += scostheta;
    _costheta_total  += scostheta;

    _ucompx          += scompx;
    _ucompy          += scompy;
    _ucompz          += scompz;
    _ucompx_total    += scompx;
    _ucompy_total    += scompy;
    _ucompz_total    += scompz;

#ifdef DDCORR
	int NDIMDP = NumbAtoms*(NumbAtoms+1)/2;
    double DipoleCorrXYZ[NDIMDP];
    double DipoleCorrX[NDIMDP];
    double DipoleCorrY[NDIMDP];
    double DipoleCorrZ[NDIMDP];
    double DipoleCorrXY[NDIMDP];

    for (int idp = 0; idp < NDIMDP; idp++)
    {
        DipoleCorrXYZ[idp] = 0.0;
        DipoleCorrX[idp]   = 0.0;
        DipoleCorrY[idp]   = 0.0;
        DipoleCorrZ[idp]   = 0.0;
        DipoleCorrXY[idp]  = 0.0;
    }

    GetDipoleCorrelationPIMC(DipoleCorrXYZ, DipoleCorrX, DipoleCorrY, DipoleCorrZ, DipoleCorrXY);

    for (int idp = 0; idp < NDIMDP; idp++)
    {
        _cdipoleXYZ[idp] += DipoleCorrXYZ[idp];
        _cdipoleX[idp]   += DipoleCorrX[idp];
        _cdipoleY[idp]   += DipoleCorrY[idp];
        _cdipoleZ[idp]   += DipoleCorrZ[idp];
        _cdipoleXY[idp]  += DipoleCorrXY[idp];
        _cdipoleXYZ_total[idp] += DipoleCorrXYZ[idp];
        _cdipoleX_total[idp]   += DipoleCorrX[idp];
        _cdipoleY_total[idp]   += DipoleCorrY[idp];
        _cdipoleZ_total[idp]   += DipoleCorrZ[idp];
        _cdipoleXY_total[idp]  += DipoleCorrXY[idp];
    }
#endif

	double srot;
	if (ROTATION)
	{
		if(MCAtom[IMTYPE].molecule == 1)
		{
			srot      = GetRotEnergy();           // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}

		if(MCAtom[IMTYPE].molecule == 3)
		{
			srot          = GetRotE3D();
		}
        
		if(MCAtom[IMTYPE].molecule == 4)
		{
			srot      = GetRotPlanarEnergy();     // kin energy
		}
		_brot        += srot;
		_rot_total   += srot;

		_brotsq      += ErotSQ;
		_rotsq_total += ErotSQ;

//Tapas commented out
	     GetRCF(); 

	}

	double stotal     = skin + srot + spot;
	_btotal          += stotal;                   // kin+pot
	_total           += stotal;

#ifdef IOWRITE
// accumulate terms for Cv
   double sCv;
   sCv = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin - spot - srot)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin) + Erot_termSQ - ErotSQ;

   _bCv += sCv;
   _Cv_total += sCv;

// accumulate terms for translational Cv
   double sCv_trans;
   double sCv_trans_1=-((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin);
   double sCv_trans_2=(2.0/MCBeta)*0.5*(double)(NDIM*NumbAtoms)/MCTau - skin;
   sCv_trans = -0.5*(double)(NDIM*NumbAtoms)/(MCBeta*MCTau)
         -((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)*((double)(NDIM*NumbAtoms)*0.5/MCTau - skin)
         + (2.0/MCBeta)*(0.5*(double)(NDIM*NumbAtoms)/MCTau - skin);

   _bCv_trans += sCv_trans;
   _Cv_trans_total += sCv_trans;
   _Cv_trans_1_total += sCv_trans_1;
   _Cv_trans_2_total += sCv_trans_2;

// accumulate terms for rotational Cv
   double sCv_rot;
   sCv_rot = -srot*srot + Erot_termSQ - ErotSQ;

   _bCv_rot += sCv_rot;
   _Cv_rot_total += sCv_rot;

// check whether there're bosons in the system

	if (BOSONS) 
	{
		GetExchangeLength();
		GetPermutation();
	}

#ifdef DEBUG_PIMC
	if (NDIM !=3) 
		nrerror("MCGetAverage()","Area estimators for 3D case only"); 
#endif
	if (BOSONS and MCAtom[IMTYPE].molecule == 1) 
	GetAreaEstimators();

	if (BOSONS)
	{
		int iframe = 0;
		GetAreaEstim3D(iframe);
	}

	if ((BOSONS and MCAtom[IMTYPE].molecule == 2) && ISPHER == 0) 
	{
		int iframe = 1;
		GetAreaEstim3D(iframe);
	}
#endif

//  reflect for MF molecule
	if(IREFLY == 1)
	{
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_XZ();
        }
//      else
//      {
//      	cout<<"Not in Reflect_MF_XZ"<<endl;
//      }
	}

    if(IREFLX == 1)
    {
    	double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
        	Reflect_MF_YZ();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_YZ"<<endl;
//      }
    }

    if(IREFLZ == 1)
    {
        double rndrfl=rnd7();
        if(rndrfl < 0.5)
        {
            Reflect_MF_XY();
        }
//      else
//      {
//          cout<<"Not in Reflect_MF_XY"<<endl;
//      }
    }

    if(IROTSYM == 1)
    {
        double rndrot=rnd7();

        if(rndrot < 0.5)
        RotSymConfig();
    }
}
#endif

void MCWormAverageReset(void)
{
}

void MCWormAverage(void)
{
}

void MCSaveBlockAverages(long int blocknumb) 
{
  	const char *_proc_=__func__;    //   MCSaveBlockAverages()

  	stringstream bc;                // convert block # to string
  	bc.width(IO_BLOCKNUMB_WIDTH); 
  	bc.fill('0');
  	bc<<blocknumb;

  	string fname = MCFileName + bc.str();  // file name prefix including block #

  	//-----------------------------------------------------------------
  	// densities

  	// toby temporarily by-pass the density save for the non-linear rotor
#ifdef IOWRITE 
  	if( IMPURITY && MCAtom[IMTYPE].molecule == 1)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);       
      	SaveDensities2D    (fname.c_str(),avergCount,MC_BLOCK);
    }
    if( IMPURITY && MCAtom[IMTYPE].molecule == 3)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);       
      	SaveDensities2D    (fname.c_str(),avergCount,MC_BLOCK);
    }
  	if( IMPURITY && MCAtom[IMTYPE].molecule == 2)
    {
      	SaveDensities1D    (fname.c_str(),avergCount);

      	SaveRho1D          (fname.c_str(),avergCount,MC_BLOCK);

      	// SaveDensities3D    (fname.c_str(),avergCount,MC_BLOCK); // this step takes lots of space. temporarily turnned off

      	SaveRhoThetaChi    (fname.c_str(),avergCount,MC_BLOCK);

      	// IOxyzAng(IOWrite,fname.c_str());
    }

	if (ROTATION) 
	{
    	SaveRCF            (fname.c_str(),avergCount,MC_BLOCK); 
	}
#endif
	SaveEnergy(MCFileName.c_str(),avergCount,blocknumb);
	SaveAngularDOF(MCFileName.c_str(),avergCount,blocknumb);
#ifdef DDCORR
	SaveDipoleCorr(MCFileName.c_str(),avergCount,blocknumb);
#endif

#ifdef IOWRITE
	if (BOSONS) 
    SaveExchangeLength (MCFileName.c_str(),avergCount,blocknumb);

	if (BOSONS && MCAtom[IMTYPE].molecule ==1)  
    SaveAreaEstimators (MCFileName.c_str(),avergCount,blocknumb);

	if (BOSONS)
    {
    	int iframe = 0;
      	SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
    }

 	if ((BOSONS && MCAtom[IMTYPE].molecule ==2) && ISPHER == 0)
    {
    	int iframe = 1;
      	SaveAreaEstim3D (MCFileName.c_str(),avergCount,blocknumb,iframe);
    }
#endif

}

#ifdef PIGSENTTYPE
void MCSaveBlockAveragesPIGSENT(long int blocknumb) 
{
  	const char *_proc_=__func__;    //   MCSaveBlockAverages()

  	stringstream bc;                // convert block # to string
  	bc.width(IO_BLOCKNUMB_WIDTH); 
  	bc.fill('0');
  	bc<<blocknumb;

  	string fname = MCFileName + bc.str();  // file name prefix including block #

	SaveEnergy(MCFileName.c_str(),avergCountENT,blocknumb);
	SaveAngularDOF(MCFileName.c_str(),avergCountENT,blocknumb);
#ifdef DDCORR
	SaveDipoleCorr(MCFileName.c_str(),avergCountENT,blocknumb);
#endif
    SaveTrReducedDens(MCFileName.c_str(),avergCount,blocknumb);
}
#endif

void SaveEnergy (const char fname [], double acount, long int blocknumb)
{
	const char *_proc_=__func__;    //  SaveEnergy()
 
	fstream fid;
	string fenergy;

	fenergy  = fname; 
	fenergy += IO_EXT_ENG; 

	fid.open(fenergy.c_str(),ios::app | ios::out);
	io_setout(fid);

	if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

#ifdef PIMCTYPE
	fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
	fid << setw(IO_WIDTH) << _bkin*Units.energy/avergCount << BLANK;    // potential anergy 2
	fid << setw(IO_WIDTH) << _brot*Units.energy/avergCount << BLANK;    // rot energy 5  
	fid << setw(IO_WIDTH) << _bpot*Units.energy/avergCount << BLANK;    // potential anergy 2
	fid << setw(IO_WIDTH) << _btotal*Units.energy/avergCount << BLANK;  //total energy including rot energy 
	fid << endl;
#endif
//
#ifdef PIGSTYPE
	fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
	fid << setw(IO_WIDTH) << _brot*Units.energy/avergCount << BLANK;    // rot energy 5  
	fid << setw(IO_WIDTH) << _brot1*Units.energy/avergCount << BLANK;    // rot energy 5  
	fid << setw(IO_WIDTH) << _bpot*Units.energy/avergCount << BLANK;    // potential anergy 2
	fid << setw(IO_WIDTH) << _btotal*Units.energy/avergCount << BLANK;  //total energy including rot energy 
	fid << endl;
#endif
//
#ifdef PIGSENTBOTH
	if (acount != 0.0)
	{
		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
		fid << setw(IO_WIDTH) << _brot*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _brot1*Units.energy/acount << BLANK;    // rot energy 5  
		fid << setw(IO_WIDTH) << _bpot*Units.energy/acount << BLANK;    // potential anergy 2
		fid << setw(IO_WIDTH) << _btotal*Units.energy/acount << BLANK;  //total energy including rot energy 
	    fid << endl;
	}
#endif
//
#ifdef IOWRITE
    fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;                 // block number 1 
    fid << setw(IO_WIDTH) << _bkin*Units.energy/avergCount << BLANK;    // kinetic energy 2
    fid << setw(IO_WIDTH) << _bpot*Units.energy/avergCount << BLANK;    // potential anergy 3
    fid << setw(IO_WIDTH) <<(_bkin+_bpot)*Units.energy/avergCount << BLANK;  // 4
   
    // if (ROTATION)
    fid << setw(IO_WIDTH) << _brot*Units.energy/avergCount << BLANK;    // rot energy 5
    fid << setw(IO_WIDTH) << _brotsq*(Units.energy*Units.energy)/avergCount << BLANK;    // rot energy square
    // fid << setw(IO_WIDTH) << _dbpot*Units.energy/avergCount << BLANK;    // potential differencies added by Hui Li 
    fid << setw(IO_WIDTH) <<(_bkin+_bpot+_brot)*Units.energy/avergCount << BLANK;  //total energy including rot energy 
    fid << setw(IO_WIDTH) << _bCv/avergCount << BLANK; // heat capacity
    fid << setw(IO_WIDTH) << _bCv_trans/avergCount << BLANK; // translational heat capacity
    fid << setw(IO_WIDTH) << _bCv_rot/avergCount << BLANK; // rotational heat capacity
    fid << endl;
#endif
    fid.close();
}

#ifdef DDCORR
void SaveDipoleCorr(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;    

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += "Dipole.corr";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
	
	if (acount != 0)
	{
   		fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;   
#ifdef PIGSENTBOTH
		int NumbAtoms1 = NumbAtoms/2;
#else
		int NumbAtoms1 = NumbAtoms;
#endif
		int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;

		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleXYZ[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleX[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleY[idp]/acount<<BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleZ[idp]/acount<< BLANK;
		}
		for (int idp = 0; idp < NDIMDP; idp++)
		{ 
    		fid << setw(IO_WIDTH) << _cdipoleXY[idp]/acount<< BLANK;
		}
    	fid << endl;
	}
    fid.close();
}
#endif

void SaveAngularDOF(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;    

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += ".dof";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

	if (acount != 0)
	{
        fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;   
        fid << setw(IO_WIDTH) << _bcostheta/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompx/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompy/acount << BLANK;
        fid << setw(IO_WIDTH) << _ucompz/acount << BLANK;
        fid << endl;
	}
    fid.close();
}

void SaveTrReducedDens(const char fname [], double acount, long int blocknumb)
{
    const char *_proc_=__func__;

    fstream fid;
    string fenergy;

    fenergy  = fname;
    fenergy += ".rden";

    fid.open(fenergy.c_str(),ios::app | ios::out);
    io_setout(fid);

    if (!fid.is_open()) _io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());

    fid << setw(IO_WIDTH_BLOCK) << blocknumb  << BLANK;
    fid << setw(IO_WIDTH) << _bnm/avergCount << BLANK;
    fid << setw(IO_WIDTH) << _bdm/avergCount << BLANK;
    fid << setw(IO_WIDTH) << _trOfDensitySq/avergCount << BLANK;
    fid << endl;
    fid.close();
}

void SaveSumEnergy (double acount, double numb)  // global average
{
	const char *_proc_=__func__;    //  SaveSumEnergy()
 
#ifdef PIMCTYPE
	_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
	_feng << setw(IO_WIDTH) << _kin_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _rot_total*Units.energy/acount << BLANK;   
	_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _total*Units.energy/acount << BLANK;   
#endif
//
#ifdef PIGSTYPE
	_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
	_feng << setw(IO_WIDTH) << _rot_total*Units.energy/acount << BLANK;   
	_feng << setw(IO_WIDTH) << _rot_total1*Units.energy/acount << BLANK;   
	_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _total*Units.energy/acount << BLANK;   
#endif
//
#ifdef PIGSENTBOTH
	_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
	_feng << setw(IO_WIDTH) << _rot_total*Units.energy/acount << BLANK;   
	_feng << setw(IO_WIDTH) << _rot_total1*Units.energy/acount << BLANK;   
	_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _total*Units.energy/acount << BLANK;   
#endif
//
#ifdef IOWRITE
	_feng << setw(IO_WIDTH_BLOCK) << numb <<BLANK;    
	_feng << setw(IO_WIDTH) << _kin_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) << _pot_total*Units.energy/acount << BLANK;    
	_feng << setw(IO_WIDTH) <<(_kin_total+_pot_total)*Units.energy/acount << BLANK;   

// if (ROTATION)
	{
		_feng << setw(IO_WIDTH) <<_rot_total*Units.energy/acount << BLANK;   
		_feng << setw(IO_WIDTH) <<_rotsq_total*(Units.energy*Units.energy)/acount << BLANK;   
	}
//_feng << setw(IO_WIDTH) << _dpot_total*Units.energy/acount << BLANK;   //added by Hui Li 
	_feng << setw(IO_WIDTH) <<(_kin_total+_pot_total+_rot_total)*Units.energy/acount << BLANK;  //total energy including rot  
// Cv
	double Cv = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-(_kin_total+_pot_total+_rot_total)*Units.energy/acount;
	Cv = Cv*Cv + _Cv_total*Units.energy*Units.energy/acount;
	Cv = -Cv*MCBeta/Temperature;
	_feng << setw(IO_WIDTH) <<Cv << BLANK;

	double Cv_trans = 0.5*(double)(NDIM*NumbAtoms*NumbTimes*Temperature)-_kin_total*Units.energy/acount;
	Cv_trans = Cv_trans*Cv_trans+ _Cv_trans_total*Units.energy*Units.energy/acount;
	Cv_trans = -Cv_trans*MCBeta/Temperature;
	_feng << setw(IO_WIDTH) <<Cv_trans << BLANK;

// if (ROTATION)
	{
		double Cv_rot = -_rot_total*Units.energy/acount;
		Cv_rot = Cv_rot*Cv_rot+ _Cv_rot_total*Units.energy*Units.energy/acount;
		Cv_rot = -Cv_rot*MCBeta/Temperature;
		_feng << setw(IO_WIDTH) <<Cv_rot << BLANK;
	}
#endif

//_feng << setw(IO_WIDTH) << _Cv_trans_1_total*Units.energy*Units.energy/acount << BLANK;
//_feng << setw(IO_WIDTH) << _Cv_trans_2_total*Units.energy*Units.energy/acount << BLANK;

	_feng << endl;
}

#ifdef DDCORR
void SaveSumDipoleCorr(double acount, double numb)
{
    const char *_proc_=__func__;

    _fdc << setw(IO_WIDTH_BLOCK) << numb <<BLANK;
#ifdef PIGSENTBOTH
	int NumbAtoms1 = NumbAtoms/2;
#else
	int NumbAtoms1 = NumbAtoms;
#endif
	int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;

	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleXYZ_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleX_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleY_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleZ_total[idp]/acount << BLANK;
	}
	for (int idp = 0; idp < NDIMDP; idp++)
	{ 
    	_fdc << setw(IO_WIDTH) << _cdipoleXY_total[idp]/acount << BLANK;
	}
    _fdc << endl;
}
#endif

void SaveSumAngularDOF(double acount, double numb)
{
    const char *_proc_=__func__;

    _fang << setw(IO_WIDTH_BLOCK) << numb <<BLANK;
    _fang << setw(IO_WIDTH) << _costheta_total/acount << BLANK;
    _fang << setw(IO_WIDTH) << _ucompx_total/acount << BLANK;
    _fang << setw(IO_WIDTH) << _ucompy_total/acount << BLANK;
    _fang << setw(IO_WIDTH) << _ucompz_total/acount << BLANK;
    _fang << endl;
}

void SaveInstantAngularDOF(long int numb)
{
    const char *_proc_=__func__;

#ifdef PIGSENTTYPE
	for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
   	{
   		for (int it = ((NumbRotTimes - 1)/2-1); it <= ((NumbRotTimes-1)/2); it++) 
		{
        	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
            int t0      = offset0 + it;

			_fangins << setw(IO_WIDTH) << MCCosine[2][t0] << BLANK;
        }
    }
	_fangins << endl;
#endif
//
#ifdef PIMCTYPE
   	_fanginsc << setw(IO_WIDTH) << numb << BLANK;
	double cosTheta   = 0.0;
	double compxyz[NDIM];
	compxyz[0] = 0.0;
	compxyz[1] = 0.0;
	compxyz[2] = 0.0;

	GetCosThetaPIMC(cosTheta, compxyz); //Centroid distribution
   	_fanginsc << setw(IO_WIDTH) << cosTheta << BLANK;
   	_fanginsc << setw(IO_WIDTH) << compxyz[0] << BLANK;
   	_fanginsc << setw(IO_WIDTH) << compxyz[1] << BLANK;
   	_fanginsc << setw(IO_WIDTH) << compxyz[2] << BLANK;
    _fanginsc << endl;
//
   	//_fangins << setw(IO_WIDTH) << numb << BLANK;
	for (int atom0 = 0; atom0 < NumbAtoms; atom0++)
   	{
   		for (int it = 0; it < NumbRotTimes; it++) // Rotational Time slices, P
		{
        	int offset0 = MCAtom[IMTYPE].offset + NumbRotTimes*atom0;
            int t0      = offset0 + it;

			_fangins << setw(IO_WIDTH) << MCCosine[0][t0] << BLANK;
			_fangins << setw(IO_WIDTH) << MCCosine[1][t0] << BLANK;
			_fangins << setw(IO_WIDTH) << MCCosine[2][t0] << BLANK;
			_fangins << endl;
        }
    }
#endif

#ifdef PIGSTYPE
#ifdef BINARY
	double instArray[5];
	instArray[0] = numb;
	double cosTheta   = 0.0;
	double compxyz[NDIM];
	GetCosThetaPIGS(cosTheta, compxyz);
	instArray[1] = cosTheta;
	instArray[2] = compxyz[0];
	instArray[3] = compxyz[1];
	instArray[4] = compxyz[2];
	_fangins.write((char*)instArray, 5*sizeof(double));
#else
#ifndef GAUSSIANMOVE
   	_fangins << setw(IO_WIDTH) << numb << BLANK;
	double cosTheta   = 0.0;
	double compxyz[NDIM];
	compxyz[0] = 0.0;
	compxyz[1] = 0.0;
	compxyz[2] = 0.0;

	GetCosThetaPIGS(cosTheta, compxyz);
   	_fangins << setw(IO_WIDTH) << cosTheta << BLANK;
   	_fangins << setw(IO_WIDTH) << compxyz[0] << BLANK;
   	_fangins << setw(IO_WIDTH) << compxyz[1] << BLANK;
   	_fangins << setw(IO_WIDTH) << compxyz[2] << BLANK;
    _fangins << endl;
#else
   	_fangins << setw(IO_WIDTH) << numb << BLANK;
	int it = (NumbTimes - 1)/2;
	int atom0 = 0;
	int t0 = it + atom0*NumbTimes;
	for (int id = 0; id < NDIM; id++)
	{
   		_fangins << setw(IO_WIDTH) << MCCoords[id][t0]<< BLANK;
	}
	for (int id = 0; id < NDIM; id++)
	{
   		_fangins << setw(IO_WIDTH) << MCCoords[id][0]<< BLANK;
	}
	for (int id = 0; id < NDIM; id++)
	{
   		_fangins << setw(IO_WIDTH) << MCCoords[id][NumbTimes-1]<< BLANK;
	}
   	_fangins << endl;
#endif
#endif 
#endif 

}

void SaveInstantEnergy()
{
    const char *_proc_=__func__;

    double srotinst, spotinst, stotalinst;
    if (MCAtom[IMTYPE].molecule == 4)
    {
#ifdef PIGSTYPE
    srotinst   = GetRotEnergyPIGS();
	spotinst   = GetPotEnergyPIGS(); 
	stotalinst = GetTotalEnergy();
#endif
//
#ifdef PIMCTYPE
    srotinst   = GetRotPlanarEnergy(); 
	spotinst   = GetPotEnergy_Densities(); 
    stotalinst = 0.0;
#endif
    }
    _fengins << setw(IO_WIDTH) << spotinst << BLANK;
    _fengins << setw(IO_WIDTH) << stotalinst << BLANK;
    _fengins << setw(IO_WIDTH) << srotinst << BLANK;
    _fengins << endl;
}

#ifdef SWAPTOUNSWAP
void SaveInstantSwap(double numbSwap, double numbUnSwap, long int numb)
{
    const char *_proc_=__func__;

    _fswap << setw(IO_WIDTH) << numb << BLANK;
    _fswap << setw(IO_WIDTH) << numbUnSwap << BLANK;
    _fswap << setw(IO_WIDTH) << numbSwap << BLANK;
    _fswap << endl;
}
#endif

void SaveSumTrReducedDens(double acount, double numb)
{
    const char *_proc_=__func__;

    _fentropy << setw(IO_WIDTH_BLOCK) << numb <<BLANK;
    _fentropy << setw(IO_WIDTH) << _nm_total/acount << BLANK;
    _fentropy << setw(IO_WIDTH) << _dm_total/acount << BLANK;
    _fentropy << setw(IO_WIDTH) << _trOfDensitySq_total/acount << BLANK;
    _fentropy << endl;
}

void InitTotalAverage(void)  // DUMP
{
	const char *_proc_=__func__;    

	totalCount = 0.0;  // need to save in the status file 
#ifdef PIGSENTBOTH
	totalCountENT = 0.0;  // need to save in the status file 
#endif
   	sumsCount = 0.0;  // counter for accum sums
   	totalStep = 0.0;

	_kin_total = 0.0;   // need a function to reset all global average
	_pot_total = 0.0;
	_total = 0.0;
	_costheta_total= 0.0;
	_ucompx_total  = 0.0;
	_ucompy_total  = 0.0;
	_ucompz_total  = 0.0;
#ifdef PIGSENTBOTH
	int NumbAtoms1 = NumbAtoms/2;
#else
	int NumbAtoms1 = NumbAtoms;
#endif
	int NDIMDP = NumbAtoms1*(NumbAtoms1+1)/2;

    _cdipoleXYZ_total.resize(NDIMDP);
    _cdipoleX_total.resize(NDIMDP);
    _cdipoleY_total.resize(NDIMDP);
    _cdipoleZ_total.resize(NDIMDP);
  	_cdipoleXY_total.resize(NDIMDP);
	
	for (int idp = 0; idp < NDIMDP; idp++)
	{
		 _cdipoleXYZ_total[idp] = 0.0;
		 _cdipoleX_total[idp] = 0.0;
		 _cdipoleY_total[idp] = 0.0;
		 _cdipoleZ_total[idp] = 0.0;
		 _cdipoleXY_total[idp] = 0.0;
	}
    _nm_total  = 0.0;
    _dm_total  = 0.0;
    _trOfDensitySq_total = 0.0;
	_dpot_total = 0.0;  //added by Hui Li

	_rot_total = 0.0;
	_rot_total1 = 0.0;
	_rotsq_total = 0.0;
	_Cv_total = 0.0;
	_Cv_trans_total = 0.0;
	_Cv_trans_1_total = 0.0;
	_Cv_trans_2_total = 0.0;
	_Cv_rot_total = 0.0;

	// open files for output
	// ENERGY
//
#ifndef PIGSENTTYPE
	string fenergy;

	fenergy  = MCFileName + IO_SUM; 
	fenergy += IO_EXT_ENG; 
 
	if (FileExist(fenergy.c_str()))   // backup the output of previous simulations 
	IOFileBackUp(fenergy.c_str());

	_feng.open(fenergy.c_str(), ios::out);
	io_setout(_feng);

	if (!_feng.is_open())
	_io_error(_proc_,IO_ERR_FOPEN,fenergy.c_str());
//

    string fangular;

    fangular  = MCFileName + IO_SUM;
    fangular += ".dof";

    if (FileExist(fangular.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fangular.c_str());

    _fang.open(fangular.c_str(), ios::out);
    io_setout(_fang);

    if (!_fang.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fangular.c_str());

//
#ifdef DDCORR
    string fdipolecorr;

    fdipolecorr  = MCFileName + IO_SUM;
    fdipolecorr += "Dipole.corr";

    if (FileExist(fdipolecorr.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fdipolecorr.c_str());

    _fdc.open(fdipolecorr.c_str(), ios::out);
    io_setout(_fdc);

    if (!_fdc.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fdipolecorr.c_str());
#endif
#endif

//
#ifdef INSTANT
#ifdef PIMCTYPE
    string fangularinsc;

    fangularinsc  = MCFileName + "_instant_centroid";
    fangularinsc += ".dof";

    if (FileExist(fangularinsc.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fangularinsc.c_str());

#ifdef BINARY
    _fanginsc.open(fangularinsc.c_str(), ios::out | ios::binary);
#else
    _fanginsc.open(fangularinsc.c_str(), ios::out);
#endif
    io_setout(_fanginsc);

    if (!_fanginsc.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fangularinsc.c_str());
#endif
//
    string fangularins;

    fangularins  = MCFileName + "_instant";
    fangularins += ".dof";

    if (FileExist(fangularins.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fangularins.c_str());

#ifdef BINARY
    _fangins.open(fangularins.c_str(), ios::out | ios::binary);
#else
    _fangins.open(fangularins.c_str(), ios::out);
#endif
    io_setout(_fangins);

    if (!_fangins.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fangularins.c_str());
#endif

#ifdef IOFILES
    string fenergyins;

    fenergyins  = MCFileName + "_instant";
    fenergyins += ".eng";

    if (FileExist(fenergyins.c_str()))   // backup the output of previous simulations
    IOFileBackUp(fenergyins.c_str());

    _fengins.open(fenergyins.c_str(), ios::out);
    io_setout(_fengins);

    if (!_fengins.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fenergyins.c_str());
#endif

#ifdef PIGSENTTYPE
    string fentropy;

    fentropy  = MCFileName + IO_SUM;
    fentropy += ".rden";

    if (FileExist(fentropy.c_str()))   // backup the output of previous simulations 
    IOFileBackUp(fentropy.c_str());

    _fentropy.open(fentropy.c_str(), ios::out);
    io_setout(_fentropy);

    if (!_fentropy.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fentropy.c_str());


#ifdef INSTANT
    string fswapins;

    fswapins  = MCFileName + "SwapUnSwapRatio";
    fswapins += ".inst";

    if (FileExist(fswapins.c_str()))   // backup the output of previous simulations
    IOFileBackUp(fswapins.c_str());

    _fswap.open(fswapins.c_str(), ios::out);
    io_setout(_fswap);

    if (!_fswap.is_open())
    _io_error(_proc_,IO_ERR_FOPEN,fswapins.c_str());
#endif
#endif
}

void DoneTotalAverage(void)
{
  _feng.close();
  _fang.close();
}

void MCSaveAcceptRatio(long int step,long int pass,long int block)
{
   int w = 8; 
 
   cout << "BLOCK:" << setw(w) << block << BLANK;
   cout << "PASS:"  << setw(w) << pass  << BLANK;
   cout << "STEP:"  << setw(w) << step  << BLANK;

#ifdef MOVECOM
   for (int type=0;type<NumbTypes;type++)
   if (WORM && (type == Worm.type))
   {
      double ratio_open    = QWAccep[0][QW_OPEN] /QWTotal[0][QW_OPEN];
      double ratio_close   = QWAccep[0][QW_CLOSE]/QWTotal[0][QW_CLOSE];

      double ratio_advance = QWAccep[0][QW_ADVANCE]/QWTotal[0][QW_ADVANCE];
      double ratio_recede  = QWAccep[0][QW_RECEDE] /QWTotal[0][QW_RECEDE];

      double ratio_swap    = QWAccep[0][QW_SWAP]/QWTotal[0][QW_SWAP];

      cout<<setw(w)<<"open/close";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_OPEN]/countQW<<"-"<<QWTotal[0][QW_CLOSE]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
 
      cout<<setw(w)<<ratio_open   << BLANK; // accept ratio for "open"   move
      cout<<setw(w)<<ratio_close  << BLANK; // accept ratio for "close"  move

      cout<<setw(w)<<"advance/recede";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_ADVANCE]/countQW<<"-"<<QWTotal[0][QW_RECEDE]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
 
      cout<<setw(w)<<ratio_advance<< BLANK;  // accept ratio for "advance" move 
      cout<<setw(w)<<ratio_recede << BLANK;  // accept ratio for "recede"  move

      cout<<setw(w)<<"swap";
      cout<<" [ "; 
      cout<<QWTotal[0][QW_SWAP]/countQW; 
      cout<<" ] "; 
      cout<<BLANK;                     
      cout<<setw(w)<<ratio_swap << BLANK;  // accept ratio for "recede"  move
   }
   else
   {
#ifndef GAUSSIANMOVE
      double ratio_molec = MCAccep[type][MCMOLEC]/MCTotal[type][MCMOLEC];
      double ratio_multi = MCAccep[type][MCMULTI]/MCTotal[type][MCMULTI];
 
      cout<<setw(w)<<MCAtom[type].type<<BLANK; // atom type

      cout<<setw(w)<<ratio_molec<<BLANK;       // accept ratio for "molecular" move
      cout<<setw(w)<<ratio_multi<<BLANK;       // accept ratio for multilevel move 
#endif
   }
#endif

   if (ROTATION)
   {
      double ratio_rotat = MCAccep[IMTYPE][MCROTAT]/MCTotal[IMTYPE][MCROTAT];
      cout<<BLANK<< "Rot: ";
      cout<<setw(w)<<ratio_rotat<<BLANK;       // accept ratio for "molecular" move
   }

   cout << endl;

}

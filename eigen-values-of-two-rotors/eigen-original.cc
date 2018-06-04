#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include "eigen.h"
#include "setup.h"
#include "input.h"

using namespace std;

int main()
{
	vinit_();

	IOReadParams(FINPUT);

    double cosTheta[nSize];
	double weightsTheta[nSize];
	thetagrid(nSize, cosTheta, weightsTheta);

	normalizationCheck(nSize, nSizeRot);

    double phi[2*nSize+1];
	double weightsPhi;
	phigrid(nSize, phi, weightsPhi);
	
	int nSizeTotal = sizeloop(nSizeRot, nSkip);

	double hMatrixRe[nSizeTotal][nSizeTotal];
	double hMatrixIm[nSizeTotal][nSizeTotal];

//Bra loop started here//

	int indexBra = 0;

	for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
	{
		if (i1Bra%nSkip) continue;

		for (int i2Bra = 0; i2Bra<nSizeRot; i2Bra++)
		{
			if (i2Bra%nSkip) continue;

			for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
			{
				int phase1Bra;
				if (m1Bra>0) phase1Bra = pow((-1.),(m1Bra));
        		else phase1Bra=1.;

				for (int m2Bra = -i2Bra; m2Bra<=i2Bra; m2Bra++)
				{
					int phase2Bra;
					if (m2Bra>0) phase2Bra = pow((-1.),(m2Bra));
        			else phase2Bra=1.;

//Ket loop started here//

					int indexKet = 0;

					for (int i1Ket = 0; i1Ket<nSizeRot; i1Ket++)
					{		
						if (i1Ket%nSkip) continue;

						for (int i2Ket = 0; i2Ket<nSizeRot; i2Ket++)
						{
							if (i2Ket%nSkip) continue;

							for (int m1Ket = -i1Ket; m1Ket<=i1Ket; m1Ket++)
							{
								int phase1Ket;
								if (m1Ket>0) phase1Ket = pow((-1.),(m1Ket));
        						else phase1Ket=1.;

								for (int m2Ket = -i2Ket; m2Ket<=i2Ket; m2Ket++)
								{
									int phase2Ket;
									if (m2Ket>0) phase2Ket = pow((-1.),(m2Ket));
        							else phase2Ket=1.;

//loop over coordinates start here//

                                 	double sumRotor1Re 		 = 0.0;
                                   	double sumRotor1Im 		 = 0.0;

									for (int iXRotor1 = 0; iXRotor1<nSize; iXRotor1++)
				                    {
                				        double xRotor1        = cosTheta[iXRotor1];
										double theta1         = acos(xRotor1);

										double associatedLegBra1 = phase1Bra*normalizedPlm(i1Bra, m1Bra, xRotor1)*sqrt(weightsTheta[iXRotor1]);
										double associatedLegKet1 = phase1Ket*normalizedPlm(i1Ket, m1Ket, xRotor1)*sqrt(weightsTheta[iXRotor1]);

                                   		double sumPhiRotor1Re = 0.0;
                                   		double sumPhiRotor1Im = 0.0;

										for (int iPhiRotor1 = 0; iPhiRotor1<2*nSize+1; iPhiRotor1++)
										{
											double spherHarmonicsReBra1  = associatedLegBra1*cos((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);
											double spherHarmonicsImBra1  = associatedLegBra1*sin((double)m1Bra*phi[iPhiRotor1])*sqrt(weightsPhi);
											double spherHarmonicsReKet1  = associatedLegKet1*cos((double)m1Ket*phi[iPhiRotor1])*sqrt(weightsPhi);
											double spherHarmonicsImKet1  = associatedLegKet1*sin((double)m1Ket*phi[iPhiRotor1])*sqrt(weightsPhi);
                                       		double sumRotor2Re 		 = 0.0;
                                       		double sumRotor2Im 		 = 0.0;

											for (int iXRotor2 = 0; iXRotor2<nSize; iXRotor2++)
                                    		{
                                        		double xRotor2        = cosTheta[iXRotor2];
												double theta2            = acos(xRotor2);

                                        		double associatedLegBra2 = phase2Bra*normalizedPlm(i2Bra, m2Bra, xRotor2)*sqrt(weightsTheta[iXRotor2]);
                                        		double associatedLegKet2 = phase2Ket*normalizedPlm(i2Ket, m2Ket, xRotor2)*sqrt(weightsTheta[iXRotor2]);

                                        		double sumPhiRotor2Re = 0.0;
                                        		double sumPhiRotor2Im = 0.0;

                                        		for (int iPhiRotor2 = 0; iPhiRotor2<2*nSize+1; iPhiRotor2++)
                                        		{
                                            		double spherHarmonicsReBra2 = associatedLegBra2*cos((double)m2Bra*phi[iPhiRotor2])*sqrt(weightsPhi);
                                            		double spherHarmonicsImBra2 = associatedLegBra2*sin((double)m2Bra*phi[iPhiRotor2])*sqrt(weightsPhi);
                                            		double spherHarmonicsReKet2 = associatedLegKet2*cos((double)m2Ket*phi[iPhiRotor2])*sqrt(weightsPhi);
                                            		double spherHarmonicsImKet2 = associatedLegKet2*sin((double)m2Ket*phi[iPhiRotor2])*sqrt(weightsPhi);

//potential is called here
													double potl;
													double dihedralAngle		= phi[iPhiRotor1]-phi[iPhiRotor2];
													vh2h2_(&rCOM, &bondLength1, &bondLength2, &theta1, &theta2, &dihedralAngle, &potl);
#ifndef TESTNORM 
													double vpot 				= potl/wavenumber;
#else
													double vpot					=1.0;
#endif
//potential call is end here

													double psiReBra				= spherHarmonicsReBra1*spherHarmonicsReBra2-spherHarmonicsImBra1*spherHarmonicsImBra2;
													double psiImBra				= spherHarmonicsReBra1*spherHarmonicsImBra2+spherHarmonicsImBra1*spherHarmonicsReBra2;
													double psiReKet				= spherHarmonicsReKet1*spherHarmonicsReKet2-spherHarmonicsImKet1*spherHarmonicsImKet2;
													double psiImKet				= spherHarmonicsReKet1*spherHarmonicsImKet2+spherHarmonicsImKet1*spherHarmonicsReKet2;

                                            		sumPhiRotor2Re          += (psiReBra*psiReKet+psiImBra*psiImKet)*vpot;
                                            		sumPhiRotor2Im          += (psiReBra*psiImKet-psiReKet*psiImBra)*vpot;
                                        		}
												sumRotor2Re          += sumPhiRotor2Re;
												sumRotor2Im          += sumPhiRotor2Im;
                                    		}
                                           	sumPhiRotor1Re   += sumRotor2Re;
                                           	sumPhiRotor1Im   += sumRotor2Im;

										}
										sumRotor1Re          += sumPhiRotor1Re;
										sumRotor1Im          += sumPhiRotor1Im;
									}

//loop over coordinates end here//


									double energyRotor;
									if (i1Bra == i1Ket && m1Bra == m1Ket && i2Bra == i2Ket && m2Bra == m2Ket) 
									{
										energyRotor = rotEnergy(i1Bra) + rotEnergy(i2Bra);
									}
									else
									{
										energyRotor = 0.0;
									}

#ifndef TESTNORM 
#ifdef EXCLUDEROTATION
									hMatrixRe[indexBra][indexKet] = sumRotor1Re;
									hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#else
									hMatrixRe[indexBra][indexKet] = sumRotor1Re + energyRotor;
									hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#endif
#else
									hMatrixRe[indexBra][indexKet] = sumRotor1Re;
									hMatrixIm[indexBra][indexKet] = sumRotor1Im;
#endif

#ifdef IOREAD
									if (sumRotor1Re > 1e-15)
									{
									cout<<" i1Bra, i1Ket, m1Bra, m1Ket, i2Bra, i2Ket, m2Bra, m2Ket "<<i1Bra<<"  "<<i1Ket<<"  "<<m1Bra<<"  "<<m1Ket<<"  "<<i2Bra<<"  "<<i2Ket<<"  "<<m2Bra<<"  "<<m2Ket<<"  "<< setw(10)<<sumRotor1Re<<"     "<<setw(10)<<sumRotor1Im<<endl;
									}
#endif

									indexKet += 1;
								}
							}
						}
					}

//Ket loop ended here//

					indexBra += 1;
				}
			}
		}
	}

//Bra loop ended here//
	
	double c[nSizeTotal*nSizeTotal];
	double d[nSizeTotal*nSizeTotal];
	int ii = 0;
	for (int i=0; i<nSizeTotal; i++)
	{
		for (int j=0; j<nSizeTotal; j++)
		{
			c[ii] = hMatrixRe[i][j];
			d[ii] = hMatrixIm[i][j];
			if(hMatrixRe[i][j] != hMatrixRe[j][i])
			{
				cout<<i<<"   "<<j<<"   "<<setw(20)<<hMatrixRe[i][j]<<"      "<<setw(20)<<hMatrixRe[j][i]<<endl;
			}
			ii++;
		}
	}
    matdiaghermitian(nSizeTotal, c, d);
	return 0;
}

void matdiagsymmetric(int nSizeTotal, double *a)
{
//Lapack starts here
    /* Locals */

	int N = nSizeTotal;
	int LDA = N;
	a[LDA*N];
    int n = N, lda = LDA, info, lwork;
    double wkopt;
    double* work;

    /* Local arrays */
    double w[N];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    char msg1[] = "Vectors";
    char msg2[] = "Upper";
    dsyev_( msg1, msg2, &n, a, &lda, w, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );

    /* Solve eigenproblem */
    dsyev_( msg1, msg2, &n, a, &lda, w, work, &lwork, &info );

    /* Check for convergence */
    if(info>0)
    {
        cout<<"The algorithm failed to compute eigenvalues.\n";
        exit(1);
    }

    /* Print eigenvalues */
    char msg3[] = "Eigenvalues";
    print_matrix(msg3, 1, n, w, 1);

    /* Print eigenvectors */
/*    char msg4[] = "Eigenvectors (stored columnwise)";
    print_matrix(msg4, n, n, a, lda );*/

    /* Free workspace */
    free((void*)work);
    exit(0);

// End of DSYEV Example 
}

void matdiaghermitian(int nSizeTotal, double *c, double *d)
{
//Lapack starts here
    /* Locals */

	c[nSizeTotal*nSizeTotal];
	d[nSizeTotal*nSizeTotal];

	int N = nSizeTotal;
	int LDA = N;
    int n = N, lda = LDA, info, lwork;
	dcomplex a[LDA*N];

	for (int i = 0; i < nSizeTotal*nSizeTotal; i++)
	{
		a[i] = {c[i],d[i]};
	}
	dcomplex wkopt;
	dcomplex* work;

    /* Local arrays */
    double w[N], rwork[3*N-2];

    /* Query and allocate the optimal workspace */
    lwork = -1;
    char msg1[] = "Vectors";
    char msg2[] = "Lower";
	zheev_( msg1, msg2, &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
	lwork = (int)wkopt.re;
	work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );

    /* Solve eigenproblem */
    zheev_( msg1, msg2, &n, a, &lda, w, work, &lwork, rwork, &info );

    /* Check for convergence */
    if( info > 0 ) 
	{
		cout<< "The algorithm failed to compute eigenvalues."<<endl;
        exit( 1 );
    }

    /* Print eigenvalues */
    char msg3[] = "Eigenvalues";
    char msg4[] = "Eigenvectors (stored columnwise)";

	print_rmatrix( msg3, 1, n, w, 1 );
	/* Print eigenvectors */
//    print_matrix( msg4, n, n, a, lda );

    /* Free workspace */
    free( (void*)work );
    exit( 0 );
} /* End of ZHEEV Example */

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double* a, int lda ) 
{
	int i, j;

	stringstream bc;
    bc.width(1);
    bc.fill('0');
    bc<<nSizeRot;
    string fname = "Eigen-values-nSizeRot" + bc.str()+".txt";

    ofstream myfile;
    myfile.open(fname.c_str());

    for( i = 0; i < m; i++ ) 
	{
    	for( j = 0; j < n; j++ )
		{
#ifndef TESTNORM 
			myfile<<setw(14)<< a[i+j*lda]*wavenumber <<endl;
#else
			myfile<<setw(14)<< a[i+j*lda] <<endl;
#endif
        }
	}
	myfile.close();
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, double* a, int lda )
{
    int i, j;

    stringstream bc;
	bc.width(1);
	bc.fill('0');
    bc<<nSizeRot;
	string fname = "Eigen-values-nSizeRot" + bc.str()+".txt";

    ofstream myfile;
    myfile.open(fname.c_str());

    for(i=0; i<m; i++)
    {
        for(j = 0; j<n; j++)
        {
            myfile<<setw(14)<<a[i+j*lda]*wavenumber<<endl;
        }
        cout<<endl;
    }
	myfile.close();
}

void normalizationCheck(int nSize, int nSizeRot)
{
    double cosTheta[nSize];
	double weightsTheta[nSize];
	thetagrid(nSize, cosTheta, weightsTheta);

    double phi[2*nSize+1];
	double weightsPhi;
	phigrid(nSize, phi, weightsPhi);

    ofstream myfile;
    myfile.open("Norm-Legendre-Polynomial.txt");
	for (int i = 0; i<nSize; i++)
	{
		for (int j = 0; j<nSize; j++)
		{
			double sum=0.0;
			for (int k = 0; k<nSize; k++)
			{
				double x=cosTheta[k];
				sum += normalizedPl0(i, x)*normalizedPl0(j, x)*weightsTheta[k];
			}
			myfile<<" i "<<i<<" j "<<j<<"   "<<sum<<endl;
		}
	}
	myfile.close();

    myfile.open("Norm-Spherical-Harmonics.txt");

	for (int i = 0; i<nSizeRot; i++)
	{
		for (int j = 0; j<nSizeRot; j++)
		{
			for (int m = 0; m<i+1; m++)
			{
				for (int n = 0; n<j+1; n++)
				{
	            	double sumThetaPhi=0.0;
	        	    for (int k = 0; k<nSize; k++)
    	        	{
       	        		double x=cosTheta[k];
						double associatedLegBra = normalizedPlm(i, m, x)*sqrt(weightsTheta[k]);
						double associatedLegKet = normalizedPlm(j, n, x)*sqrt(weightsTheta[k]);
						
						double sumPhi = 0.0;
	        	    	for (int l = 0; l<2*nSize+1; l++)
    	        		{
							double spherHarmonicsReBra = associatedLegBra*cos((double)m*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsImBra = associatedLegBra*sin((double)m*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsReKet = associatedLegKet*cos((double)n*phi[l])*sqrt(weightsPhi);
							double spherHarmonicsImKet = associatedLegKet*sin((double)n*phi[l])*sqrt(weightsPhi);
            	    		sumPhi += spherHarmonicsReBra*spherHarmonicsReKet+spherHarmonicsImBra*spherHarmonicsImKet;
            			}
						sumThetaPhi += sumPhi;
            		}
            		myfile<<" i "<<i<<" j "<<j<<" m "<<m<<" n "<<n<<"   "<<sumThetaPhi<<endl;
				}
			}
		}
	}
    myfile.close();

}	

double normalizedPlm(int j,int m, double x)
{
	if (j < m) cerr<<"j < m"<<endl;
	double jv=(double)j;
	double mv=fabs((double)m);
	return sqrt((jv+.5)*exp(lgamma(jv-mv+1.)-lgamma(jv+mv+1.)))*plgndr(j,abs(m),x); 

	// for the normalization constant sqrt(...) refer to 6.55 in Quantum Chemistry. 
	//lgamma is the natural logarithim of the gamma function: gamma(n) = (n-1)! 
}

double normalizedPl0(int j,double x)
{
	return sqrt((double)j+.5)*plgndr(j,0,x);
}

void thetagrid(int nSize, double *cosTheta, double *weightsTheta)
{
	//Gauss-Legendre quadrature points and wiights are generated 
	//by gauleg.f

    double x1=-1.;
    double x2=1.;
    double x[nSize], w[nSize];
    cosTheta[nSize];
	weightsTheta[nSize];

    gauleg(x1,x2,x,w,nSize);

    ofstream myfile;
    myfile.open("DVR-theta.txt");

    for (int i=0; i<nSize; i++)
    {
		cosTheta[i]     = x[i];
		weightsTheta[i] = w[i];

        myfile<<"DVR theta = "<<cosTheta[i]<<"  weights =  "<<weightsTheta[i]<<endl;
    }
    myfile.close();
}

void phigrid(int nSize, double *phi, double &weightsPhi)
{

	const double PI  =3.141592653589793238463;
    phi[2*nSize+1];

    ofstream myfile;
    myfile.open("DVR-phi.txt");

    for (int i=0; i<2*nSize+1; i++)
    {
        phi[i]          = 2.0*PI*((double)i+1.0)/(2.0*(double)nSize+1.0);
        weightsPhi      = 1.0/(2.0*(double)nSize+1.0);

        myfile<<"DVR phi = "<<phi[i]<<"  weights =  "<<weightsPhi<<endl;
    }
    myfile.close();
}

int sizeloop(int nSizeRot, int nSkip)
{
    int indexBra = 0;

    for (int i1Bra = 0; i1Bra<nSizeRot; i1Bra++)
    {
        if (i1Bra%nSkip) continue;

        for (int i2Bra = 0; i2Bra<nSizeRot; i2Bra++)
        {
            if (i2Bra%nSkip) continue;

            for (int m1Bra = -i1Bra; m1Bra<=i1Bra; m1Bra++)
            {
                for (int m2Bra = -i2Bra; m2Bra<=i2Bra; m2Bra++)
                {
					indexBra++;
				}
			}
		}
	}
	return indexBra;
}


double rotEnergy(int l)
{
	double energyj0_v0    = -36117.5942855;  // unit cm^-1
	double energyj1_v0    = -35999.1009407;  // unit cm^-1
//	double rotConstant    = 0.5*(energyj1_v0 - energyj0_v0)/wavenumber;
	double rotConstant    = 59.0622/wavenumber;
	double energy         = rotConstant*(double)l*((double)l + 1.0);

	return energy;
}

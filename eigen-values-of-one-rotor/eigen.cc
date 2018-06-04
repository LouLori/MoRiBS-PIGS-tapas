#include <stdio.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>
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

	double hMatrix[nSizeRot][nSizeRot];

//Bra loop started here//

	for (int iBra = 0; iBra<nSizeRot; iBra++)
	{
		for (int iKet = 0; iKet<nSizeRot; iKet++)
		{		
           	double sum 		 = 0.0;
			for (int iXRotor = 0; iXRotor<nSize; iXRotor++)
            {
                double xRotor        = cosTheta[iXRotor];

				double associatedLegBra = normalizedPlm(iBra, 0, xRotor)*sqrt(weightsTheta[iXRotor]);
				double associatedLegKet = normalizedPlm(iKet, 0, xRotor)*sqrt(weightsTheta[iXRotor]);

               	sum          += associatedLegBra*associatedLegKet;
			}
            hMatrix[iBra][iKet] = sum;
		}
	}

	for (int i=0; i<nSizeRot; i++)
	{
		for (int j=0; j<nSizeRot; j++)
		{
			if(hMatrix[i][j] != hMatrix[j][i])
			{
				cout<<i<<"   "<<j<<"   "<<setw(20)<<hMatrix[i][j]<<"      "<<setw(20)<<hMatrix[j][i]<<endl;
			}
		}
	}
    exit(0);
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
//    print_matrix(msg3, 1, n, w, 1);

    /* Print eigenvectors */
/*    char msg4[] = "Eigenvectors (stored columnwise)";
    print_matrix(msg4, n, n, a, lda );*/

    /* Free workspace */
    free((void*)work);
    exit(0);

// End of DSYEV Example 
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

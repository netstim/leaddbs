
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <map>
#include <string.h>
#include <time.h>

using namespace std;

#define REAL float
#define PI 3.1415926536



#include "MersenneTwister.h"
MTRand mtrand;
#include "auxilary_classes.cpp"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 2 ) {
	printf("\nUsage: Df = STderivative(f)\n\n");
	printf(" Computes xxx\n");
	printf(" Parameters:\n");
	printf("   f - 2D input image of type REAL \n");
	printf(" Return value Df contains xxx.\n\n");
    return;
	} else if(nlhs>2) {
	printf("Too many output arguments\n");
    return;
	}
	
	int ppcnt = 0;
	const mxArray *Points;
	Points = prhs[ppcnt++];       
	int numPoints = mxGetN(Points);
    int numDiffP = mxGetM(Points)-3;
    REAL *points = (REAL*) mxGetData(Points);

    const mxArray *Param;
	Param = prhs[ppcnt++];       
    double *param = mxGetPr(Param);

    
	Particle* particles = (Particle*) malloc(sizeof(Particle)*numPoints);
	Particle particles_proposal[10000];
    REAL *diff = (REAL*) malloc(sizeof(REAL)*numPoints*numDiffP);
    REAL *diff_prop = (REAL*) malloc(sizeof(REAL)*numDiffP*10000);
    
    
    
    int i;
	for (i = 0; i < numPoints; i++)
	{
		particles[i].R.setXYZ(points[(numDiffP+3)*i],points[(numDiffP+3)*i+1],points[(numDiffP+3)*i+2]);
        for (int k = 0; k < numDiffP;k++)
            diff[i*(numDiffP)+k] = points[(numDiffP+3)*i+k+3];
	}
	pVector source = particles[0].R;
	pVector sink = particles[numPoints-1].R;
    REAL *diff_source = (REAL*) malloc(sizeof(REAL)*numDiffP);
    REAL *diff_sink = (REAL*) malloc(sizeof(REAL)*numDiffP);
    for (int k = 0; k < numDiffP;k++)
    {
        diff_source[k] = diff[k];
        diff_prop[k] = diff[k];
        diff_sink[k] = diff[k+numDiffP*(numPoints-1)];
    }


	int pcnt = numPoints;
	REAL len = param[0];
    REAL Leng = 0;

    particles_proposal[0].R = source;
    REAL dtau = 0;
    int cur_p = 1;
    int cur_i = 1;
    pVector dR;
    REAL *dD = (REAL*) malloc(sizeof(REAL)*numDiffP);
    REAL normdR;
    for (;;)
    {
        while (dtau <= len && cur_p < pcnt)
        {
            dR  = particles[cur_p].R - particles[cur_p-1].R;
            normdR = sqrt(dR.norm_square());
            dtau += normdR;
            Leng += normdR;
            for (int k = 0; k < numDiffP; k++)
               dD[k] = diff[cur_p*numDiffP + k] - diff[(cur_p-1)*numDiffP + k];
            
            cur_p++;
        }

        if (dtau >= len)
        {
            particles_proposal[cur_i].R = particles[cur_p-1].R - dR*( (dtau-len)/normdR );
            for (int k = 0; k < numDiffP; k++)
                diff_prop[cur_i*numDiffP+k] = diff[(cur_p-1)*numDiffP+k]; // -dD[k]*(dtau-len)/len;
        }
        else
        {
            for (int k = 0; k < numDiffP; k++)
                diff_prop[cur_i*numDiffP+k] = diff_sink[k];
            particles_proposal[cur_i].R = sink;
            break;
        }

        dtau = dtau-len;

        cur_i++;
        if (cur_i >= 10000)
        {
            mexPrintf("bugy");
            break;
        }			

    }




	const mwSize sz[] = {3+numDiffP, cur_i+1};

 	plhs[0] = mxCreateNumericArray(2,sz,mxGetClassID(Points),mxREAL);
 	REAL *pk = (REAL*) mxGetData(plhs[0]);	
	for (i = 0; i < cur_i+1 ; i++)
	{
		pk[(3+numDiffP)*i] = 	particles_proposal[i].R.x;
		pk[(3+numDiffP)*i+1] = particles_proposal[i].R.y;
		pk[(3+numDiffP)*i+2] = particles_proposal[i].R.z;
        for (int k = 0; k < numDiffP;k++)
        {
            pk[(3+numDiffP)*i+3+k] = diff_prop[i*numDiffP+k];
        }
	}


	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
 	double *plen =  mxGetPr(plhs[1]);	
	plen[0] = Leng;
    
    free(dD);
    free(particles);
    free(diff);
    free(diff_prop);
    free(diff_source);
    free(diff_sink);
}


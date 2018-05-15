/***************************************************************************************************************
    Copyright (c) 2011, Marco Reisert, Valerij G. Kiselev, Medical Physics, University Medical Center Freiburg
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
	* Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	* Neither the name of the <organization> nor the
	  names of its contributors may be used to endorse or promote products
	  derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****************************************************************************************************************/

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
    REAL *points = (REAL*) mxGetData(Points);

    const mxArray *Param;
	Param = prhs[ppcnt++];       
    double *param = mxGetPr(Param);

    
	Particle* particles = (Particle*) malloc(sizeof(Particle)*numPoints);
	Particle particles_proposal[10000];

    int i;
	for (i = 0; i < numPoints; i++)
	{
		particles[i].R.setXYZ(points[3*i],points[3*i+1],points[3*i+2]);
	}
	pVector source = particles[0].R;
	pVector sink = particles[numPoints-1].R;


	int pcnt = numPoints;
	REAL len = param[0];
    REAL Leng = 0;

		particles_proposal[0].R = source;
		REAL dtau = 0;
		int cur_p = 1;
		int cur_i = 1;
		pVector dR;
		REAL normdR;
		for (;;)
		{
			while (dtau <= len && cur_p < pcnt)
			{
				dR  = particles[cur_p].R - particles[cur_p-1].R;
				normdR = sqrt(dR.norm_square());
				dtau += normdR;
                Leng += normdR;
				cur_p++;
			}

			if (dtau >= len)
			{
				particles_proposal[cur_i].R = particles[cur_p-1].R - dR*( (dtau-len)/normdR );
			}
			else
			{
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




	const mwSize sz[] = {3, cur_i+1};

 	plhs[0] = mxCreateNumericArray(2,sz,mxGetClassID(Points),mxREAL);
 	REAL *pk = (REAL*) mxGetData(plhs[0]);	
	for (i = 0; i < cur_i+1 ; i++)
	{
		particles_proposal[i].R.storeXYZ();
		pk[3*i] = pVector::store[0];
		pk[3*i+1] = pVector::store[1];
		pk[3*i+2] = pVector::store[2];
	}


	plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
 	double *plen =  mxGetPr(plhs[1]);	
	plen[0] = Leng;
    
    free(particles);
}


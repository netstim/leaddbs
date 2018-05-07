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


#include <math.h>
#include "mex.h"
#include "matrix.h"
#define REAL float








void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs!=2) {
	mexPrintf("wrong usage!!\n",nrhs);
    return;
	} else if(nlhs>4) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;
    const mxArray *Dim;
    Dim = prhs[pcnt++];       
    REAL *dim = (REAL*) mxGetData(Dim);
    int w = (int) dim[0];
    int h = (int) dim[1];
    int d = (int) dim[2];    

    const mxArray *Pts;
    Pts = prhs[pcnt++];       
    const int numdim = mxGetNumberOfDimensions(Pts);
    const mwSize *pdims = mxGetDimensions(Pts);
    int numPts = pdims[1];
    REAL *pts = (REAL*) mxGetData(Pts);
    

    mwSize dims[3];
    dims[0] = w;
    dims[1] = h;
    dims[2] = d;
    plhs[0] = mxCreateNumericArray(3,dims,mxGetClassID(Dim),mxREAL);
    REAL *accu = (REAL*) mxGetData(plhs[0]);


   for (int i = 0; i < numPts; i++)
	{
		int idx = 3*i;

		int px = (int) (pts[idx]);
		if (px < 0 || px >= w-1)
			continue;
		int py = (int) (pts[idx+1]);
		if (py < 0 || py >= h-1)
			continue;
		int pz = (int) (pts[idx+2]);
		if (pz < 0 || pz >= d-1)
			continue;
		float frac_x = pts[idx  ] - px;
		float frac_y = pts[idx+1] - py;
		float frac_z = pts[idx+2] - pz;


		accu[px + w*(py+h*pz)] += (1-frac_x)*(1-frac_y)*(1-frac_z);

		accu[px+1 + w*(py+h*pz)] += (frac_x)*(1-frac_y)*(1-frac_z);
		accu[px + w*(py+1+h*pz)] += (1-frac_x)*(frac_y)*(1-frac_z);
		accu[px + w*(py+h*pz+h)] += (1-frac_x)*(1-frac_y)*(frac_z);

		accu[px + w*(py+1+h*pz+h)] += (1-frac_x)*(frac_y)*(frac_z);
		accu[px+1 + w*(py+h*pz+h)] += (frac_x)*(1-frac_y)*(frac_z);
		accu[px+1 + w*(py+1+h*pz)] += (frac_x)*(frac_y)*(1-frac_z);

		accu[px+1 + w*(py+1+h*pz+h)] += (frac_x)*(frac_y)*(frac_z);
	
	}


}







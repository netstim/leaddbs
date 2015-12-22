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
#include <string.h>
#define REAL float

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 3) {
	printf("\nUsage: res = aniosDiffusion(img,Dcoeff,[alpha it])\n\n",nrhs);
	printf(" Propagates Heat equation with position dependedd conductivity Dcoeff\n\n");
	printf(" img - inital state\n");
	printf(" Dcoeff - contains Dtensor of size [ 6 size(img)] with compnents xx yy zz xz yz xz \n");
	printf(" alpha - timestep \n");
	printf(" it - number of iterations \n");

    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}

    int pcnt = 0;
    const mxArray *Img;
    Img = prhs[pcnt++];       
    const int numdim = mxGetNumberOfDimensions(Img);
    const int *dims = mxGetDimensions(Img);
    REAL *img = (REAL*) mxGetData(Img);
    if (numdim != 3)
    {
	mexPrintf("input has to be 3dim\n");
	return;
    }
    int w = dims[0];
    int h = dims[1];
    int d = dims[2];

    const mxArray *BG;
    BG = prhs[pcnt++];       
    const int *dimsbg = mxGetDimensions(BG);
    if (dimsbg[0] != 6)
    {
	mexPrintf("blurrguide has to be rank 2\n");
	return;
    }
    REAL *blurrg = (REAL*) mxGetData(BG);


    const mxArray *Factors;
    Factors = prhs[pcnt++];          
    REAL *factors = (REAL*) mxGetData(Factors);
    REAL alpha = factors[0];
    int itmax = (int) factors[1];

    plhs[0] = mxCreateNumericArray(3,dims,mxGetClassID(Img),mxREAL);
    REAL *result = (REAL*) mxGetData(plhs[0]);
    memcpy(result,img,sizeof(REAL)*w*h*d);



    int ndims[4];
    ndims[0] = 3; ndims[1] = dims[0]; ndims[2] = dims[1]; ndims[3] = dims[2];
    const mxArray *Tmp1 = mxCreateNumericArray(4,ndims,mxGetClassID(Img),mxREAL);
    REAL *tmp1 = (REAL*) mxGetData(Tmp1);


	
    int wh = w*h;


    for (int k = 0; k < itmax ; k++)
    {
	
	for (int z = 1; z < d-1; z++)
		for (int y = 1; y < h-1; y++)
			for (int x = 1; x < w-1; x++)
			{
				int idx = z*wh + y*w + x;
				REAL gx = result[idx-1] - result[idx];
				REAL gy = result[idx-w] - result[idx];
				REAL gz = result[idx-wh] - result[idx];
				REAL *D = &(blurrg[6*idx]);
				REAL *tmp = &(tmp1[3*idx]);
				tmp[0] = gx*D[0] + gy*D[3] + gz*D[4];
				tmp[1] = gx*D[3] + gy*D[1] + gz*D[5];
				tmp[2] = gx*D[4] + gy*D[5] + gz*D[2];
			}
	
	for (int z = 1; z < d-1; z++)
		for (int y = 1; y < h-1; y++)
			for (int x = 1; x < w-1; x++)
			{
				int idx = (z*wh + y*w + x);
				REAL *tmp = &(tmp1[3*idx]);
				result[idx] = result[idx] - alpha*(tmp[3]-tmp[0] + tmp[1+3*w]-tmp[1] + tmp[2+3*wh]-tmp[2]);
			}
	/*
    
    
	for (int z = 1; z < d-1; z++)
		for (int y = 1; y < h-1; y++)
			for (int x = 1; x < w-1; x++)
			{
				int idx = z*wh + y*w + x;
				REAL gxx = result[idx+1] + result[idx-1] - 2*result[idx];
				REAL gyy = result[idx+w] + result[idx-w] - 2*result[idx];
				REAL gzz = result[idx+wh] + result[idx-wh] - 2*result[idx];
                REAL gxy = (result[idx+1+w] - result[idx-1+w] - (result[idx+1-w] - result[idx-1-w]))/4;
                REAL gxz = (result[idx+1+wh] - result[idx-1+wh] - (result[idx+1-wh] - result[idx-1-wh]))/4;
                REAL gyz = (result[idx+w+wh] - result[idx-w+wh] - (result[idx+w-wh] - result[idx-w-wh]))/4;
 
                REAL *D = &(blurrg[6*idx]);
				//tmp[0] = gxx*D[0] + gxy*D[3] + gxz*D[4];
				//tmp[1] = gxy*D[3] + gyy*D[1] + gyz*D[5];
				//tmp[2] = gxz*D[4] + gyz*D[5] + gzz*D[2];
                                                                
                result[idx] = result[idx] - alpha*( gxx*D[0] + gxy*D[3] + gxz*D[4]  +  gxy*D[3] + gyy*D[1] + gyz*D[5] + gxz*D[4] + gyz*D[5] + gzz*D[2]);
            }
     */
    }






}






/***************************************************************************************************************
    Copyright (c) 2011, Marco Reisert, Valerij G. Kiselev, Medical Physics, University Medical Center Freiburg
    All rights reserved.

    Reference: Reisert M, Mader I, Anastasopoulos C, Weigel M, Schnell S, Kiselev V. 
    Global fiber reconstruction becomes practical.
    Neuroimage.  2011 Jan 15;54(2):955-62. Epub 2010 Sep 18.


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
#include <vector>
#include <string.h>

using namespace std;

#define REAL float
#define PI 3.1415926536



#include "MersenneTwister.h"
MTRand mtrand;




#include "ParticleGrid.cpp"







void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 4) {

    return;
	} else if(nlhs>2) {
	printf("Too many output Arguments %i\n",nlhs);
    return;
	}
    
	fprintf(stderr,"computing CM \n"); fflush(stderr);
    
    
    ParticleGrid<Particle> pcontainer;

  
    
	int pcnt = 0;
	const mxArray *Labeling;
	Labeling = prhs[pcnt++];       
	REAL *labeling = (REAL*) mxGetData(Labeling);
    const mwSize *dims = mxGetDimensions(Labeling);
    

    REAL Nsize = REAL(*mxGetPr(prhs[pcnt++]));
      
  	const mxArray *Conns;
	Conns = prhs[pcnt++];       
	int numConns = mxGetN(Conns);
	REAL *conns = (REAL*) mxGetData(Conns);  
    
  	const mxArray *Lens;
	Lens = prhs[pcnt++];       
	REAL *lens = (REAL*) mxGetData(Lens);  
    
    
    
    
    
    double cellsize = Nsize*2;

    
    REAL lam = 1/(cellsize*cellsize);
    
    int cellcnt_x = (int)((REAL)dims[0]/cellsize) +1;
    int cellcnt_y = (int)((REAL)dims[1]/cellsize) +1;
    int cellcnt_z = (int)((REAL)dims[2]/cellsize) +1;
    int cell_capacity = 1024;
//   mexPrintf("%f %f %f \n",sz[0],sz[1],sz[2]);

   mexPrintf("%i %i %i \n",cellcnt_x,cellcnt_y,cellcnt_z);
    int err = pcontainer.allocate(1000000, cellcnt_x, cellcnt_y, cellcnt_z, cellsize,cell_capacity);

    if (err == -1)
    {
        fprintf(stderr,"RJMCMCBase: out of Memory!\n");
        return;
    }

    
    
    
    mexPrintf("total number fibers: %i \n",numConns/2);
    

    int attrcnt = 3;
    int maxlabel = 0;
    for (int z =0 ; z < dims[2]; z++)
    for (int y =0 ; y < dims[1]; y++)
    {
        for (int x =0 ; x < dims[0]; x++)
        {
                int idx = x+y*dims[0]+z*dims[0]*dims[1];
                if (labeling[idx] > 0)
                {
                    if (labeling[idx] > maxlabel)
                        maxlabel = labeling[idx];
                    REAL xf = REAL(x);
                    REAL yf = REAL(y);
                    REAL zf = REAL(z);
                    pVector R(xf,yf,zf);
                    Particle *p = pcontainer.newParticle(R);
                    if (p!=0) 
                    {
                        p->label = labeling[idx]-1;
                    }
                    else
                    {
                        fprintf(stderr,"error: cannot allocate particle (%i,%i,%i)\n",x,y,z);
                        return;
                    }
                }
        }
      
    }
    mexPrintf("maxlabel %i\n",maxlabel);
   
    const mwSize dims2[] = {maxlabel,maxlabel};
    plhs[0] = mxCreateNumericArray(2,dims2,mxGetClassID(Labeling),mxREAL);
	REAL *cc = (REAL*) mxGetData(plhs[0]);	
    plhs[1] = mxCreateNumericArray(2,dims2,mxGetClassID(Labeling),mxREAL);
	REAL *lenmat = (REAL*) mxGetData(plhs[1]);	
    
    vector<Particle*> n1;
    vector<Particle*> n2;
    n1.resize(2056);
    n2.resize(2056);
    
   
    for (int k = 0; k < numConns/2;k++)
    {
            
        if (k%(numConns/15) == 0)
        {
            mexPrintf("%.0f%%, ",100.0*(float)k/numConns*2);
            mexEvalString("drawnow;");
        }
        
        n1.clear();
        n2.clear();
            
        pVector R1(conns[k*6],conns[k*6+1],conns[k*6+2]);
        pcontainer.computeNeighbors(R1);		
		for (;;)
		{
			Particle *p =  pcontainer.getNextNeighbor();
			if (p == 0) break;
            n1.push_back(p);
        }      
        pVector R2(conns[k*6+3],conns[k*6+4],conns[k*6+5]);
        pcontainer.computeNeighbors(R2);		
		for (;;)
		{
			Particle *p =  pcontainer.getNextNeighbor();
			if (p == 0) break;
            n2.push_back(p);
        }
        
   
        for (int i = 0; i < n1.size();i++)
        {
             
            for(int j = 0; j < n2.size();j++)
            {
                if (n1[i]->label >= 0 && n1[i]->label < maxlabel && n2[j]->label >= 0 && n2[j]->label < maxlabel)
                {
    
                    REAL dif1 = 0;                
                    n1[i]->R.storeXYZ();
                    for (int q = 0; q < 3 ; q++)                    
                        dif1 += (conns[k*6+q]-pVector::store[q])*(conns[k*6+q]-pVector::store[q]);
                    if (dif1 > Nsize*Nsize)
                        continue;
                     REAL dif2 = 0;
                    n2[j]->R.storeXYZ();
                    for (int q = 0; q < 3 ; q++)
                        dif2 += (conns[k*6+3+q]-pVector::store[q])*(conns[k*6+3+q]-pVector::store[q]);
                    if (dif2 > Nsize*Nsize)
                        continue;


                    REAL fac = exp(-(dif1+dif2)*lam);

                    cc[n1[i]->label + maxlabel*n2[j]->label]+=fac;
                    cc[n2[j]->label + maxlabel*n1[i]->label]+=fac;

                    lenmat[n1[i]->label + maxlabel*n2[j]->label] += fac*lens[k];
                    lenmat[n2[j]->label + maxlabel*n1[i]->label] += fac*lens[k];
                }
            }
        }
    }
    
    
}



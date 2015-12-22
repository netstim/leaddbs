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

#include "auxilary_classes.cpp"


class SphereInterpolator
{
public:
    REAL *barycoords;
    int *indices;
    int size;  // size of LUT
    int sN;   // (sizeofLUT-1)/2
    int nverts; // number of data vertices


    REAL beta;
    REAL inva;
    REAL b;
    
    
    int *idx;
    REAL *interpw;
    
    SphereInterpolator(REAL *barycoords, int *indices, int numverts, int sizeLUT, REAL beta)
    {
        this->barycoords = barycoords;
        this->indices = indices;
        this->size = sizeLUT;
        this->sN = (sizeLUT-1)/2;
	this->nverts = numverts;
        this->beta = beta;

        inva = (sqrt(1+beta)-sqrt(beta));
        b = 1/(1-sqrt(1/beta + 1));    
       
    }

    
    inline void getInterpolation(pVector &N)
    {
	    N.storeXYZ();
	    REAL nx = pVector::store[0];
	    REAL ny = pVector::store[1];
	    REAL nz = pVector::store[2];

            if (nz > 0.5)
            {
                int x = real2int(nx);
                int y = real2int(ny);
                int i = 3*6*(x+y*size);  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nz < -0.5)
            {
                int x = real2int(nx);
                int y = real2int(ny);
                int i = 3*(1+6*(x+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nx > 0.5)
            {
                int z = real2int(nz);
                int y = real2int(ny);
                int i = 3*(2+6*(z+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nx < -0.5)
            {
                int z = real2int(nz);
                int y = real2int(ny);
                int i = 3*(3+6*(z+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (ny > 0)
            {
                int x = real2int(nx);
                int z = real2int(nz);
                int i = 3*(4+6*(x+z*size));  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }     
            else
            {
                int x = real2int(nx);
                int z = real2int(nz);
                int i = 3*(5+6*(x+z*size));  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }            
                        
    }
    
    
    inline REAL invrescale(REAL f)
    {
        REAL x = (fabs(f)-b)*inva;
        if (f>0)            
            return (x*x-beta);
        else
            return beta - x*x;
    }
    
    inline int real2int(REAL x)
    {
       return int((invrescale(x)+1)*sN-0.5);
    
    }

    
};




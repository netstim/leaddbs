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


#ifndef _ENCOMPINTERFACE
#define _ENCOMPINTERFACE


#include "SphereInterpolator.cpp"

inline REAL myATAN2(REAL y,REAL x)
{
        REAL phi = acos(x);
     //   REAL phi = ((x>=1.0) ? ((0.0000*x+-0.0000)) : ((x>=1.0) ? ((-10.0167*x+10.0167)) : ((x>=0.9) ? ((-3.1336*x+3.2713)) : ((x>=0.8) ? ((-1.9247*x+2.1833)) : ((x>=0.5) ? ((-1.3457*x+1.7200)) : ((x>=0.0) ? ((-1.0472*x+1.5708)) : ((x>=-0.5) ? ((-1.0472*x+1.5708)) : ((x>=-0.8) ? ((-1.3457*x+1.4216)) : ((x>=-0.9) ? ((-1.9247*x+0.9583)) : ((x>=-1.0) ? ((-3.1336*x+-0.1297)) : ((x>=-1.0) ? ((-10.0167*x+-6.8751)) : 1 )))))))))));
        if (y < 0) phi = 2*PI - phi;
        if (phi<0) phi =  phi + PI;
        return phi;
        
}


class EnergyComputerBase
{

public:

    REAL *dataimg;
    const int *datasz;
    

  
    SphereInterpolator *sinterp;
    
    ParticleGrid<Particle> *pcontainer;
    
    int w,h,d;
    REAL voxsize_w;
    REAL voxsize_h;
    REAL voxsize_d;

    int w_sp,h_sp,d_sp;    
    REAL voxsize_sp_w;
    REAL voxsize_sp_h;
    REAL voxsize_sp_d;
    
    
    int nip; // number of data vertices on sphere

    
	REAL *spatprobimg;
	REAL *cumulspatprob;
	int *spatidx;
	int scnt; 
    
    
    
    
    
    REAL eigen_energy;

    EnergyComputerBase(REAL *data, const int *dsz, double *voxsize, SphereInterpolator *sp, ParticleGrid<Particle> *pcon, REAL *spimg, int spmult)
    {
        dataimg = data;
        datasz = dsz;
        sinterp = sp;   

        spatprobimg = spimg;
        
        nip = datasz[0];
        
        
        w = datasz[1];
        h = datasz[2];
        d = datasz[3];

        voxsize_w = voxsize[0];
        voxsize_h = voxsize[1];
        voxsize_d = voxsize[2];
       
        
        w_sp = datasz[1]*spmult;
        h_sp = datasz[2]*spmult;
        d_sp = datasz[3]*spmult;
        
        voxsize_sp_w = voxsize[0]/spmult;
        voxsize_sp_h = voxsize[1]/spmult;
        voxsize_sp_d = voxsize[2]/spmult;
        

        fprintf(stderr,"Data size (voxels) : %i x %i x %i\n",w,h,d);
        fprintf(stderr,"voxel size:  %f x %f x %f\n",voxsize_w,voxsize_h,voxsize_d);
        fprintf(stderr,"mask_oversamp_mult: %i\n",spmult);
        
        if (nip != sp->nverts)
        {
            fprintf(stderr,"EnergyComputer: error during init: data does not match with interpolation scheme\n");
        }

        pcontainer = pcon;
        
  
		int totsz = w_sp*h_sp*d_sp;
		cumulspatprob = (REAL*) malloc(sizeof(REAL) * totsz);
		spatidx = (int*) malloc(sizeof(int) * totsz);
        if (cumulspatprob == 0 || spatidx == 0)
        {
            fprintf(stderr,"EnergyCOmputerBase: out of memory!\n");
            return ;
        }
        
        
		scnt = 0;
		cumulspatprob[0] = 0;
		for (int x = 0; x < w_sp;x++)
			for (int y = 0; y < h_sp;y++)
				for (int z = 0; z < d_sp;z++)
				{
					int idx = x+(y+z*h_sp)*w_sp;
					if (spatprobimg[idx] > 0.5)
					{
						cumulspatprob[scnt+1] = cumulspatprob[scnt] + spatprobimg[idx];
						spatidx[scnt] = idx;
						scnt++;
					}
				}

		for (int k = 0; k < scnt; k++)
		{
			cumulspatprob[k] /= cumulspatprob[scnt];
		}		

		fprintf(stderr,"#active voxels: %i (in mask units) \n",scnt);        
       
        
    }
    
    ~EnergyComputerBase()
    {
        free(cumulspatprob);        
		free(spatidx);
    }

	virtual void setParameters()
	{
			mexPrintf("in set base\n");
	}
	

	
	void drawSpatPosition(pVector *R)
	{
			
		REAL r = mtrand.frand();		
		int j;
		int rl = 1;
		int rh = scnt;
		while(rh != rl)
		{
			j = rl + (rh-rl)/2;
			if (r < cumulspatprob[j])
				{
				rh = j;
				continue;
				}
			if (r > cumulspatprob[j])
				{
				rl = j+1;
				continue;
				}
			break;
		}	
		R->setXYZ(voxsize_sp_w*((REAL)(spatidx[rh-1] % w_sp)  + mtrand.frand()),
		          voxsize_sp_h*((REAL)((spatidx[rh-1]/w_sp) % h_sp)  + mtrand.frand()),
		          voxsize_sp_d*((REAL)(spatidx[rh-1]/(w_sp*h_sp))    + mtrand.frand()));	
	
	}
	
	
	REAL SpatProb(pVector R)
	{
			
		R.storeXYZ();
		int rx = int(pVector::store[0]/voxsize_sp_w);
		int ry = int(pVector::store[1]/voxsize_sp_h);
		int rz = int(pVector::store[2]/voxsize_sp_d);
		if (rx >= 0 && rx < w_sp && ry >= 0 && ry < h_sp && rz >= 0 && rz < d_sp)
			return spatprobimg[rx + w_sp* (ry + h_sp*rz)];
		else
			return 0;
	
	}

	


/*	
	inline REAL evaluateODF(pVector &R, pVector &N, REAL &len)
	{
		const int CU = 10;
		pVector Rs;
		REAL Dn = 0;
		int xint,yint,zint,spatindex;

		sinterp->getInterpolation(N);
		for (int i=-CU; i < CU;i++)
		{
			Rs = R + (N * len) * ((REAL)i/CU);
			xint = int(Rs.x);
			yint = int(Rs.y);
			zint = int(Rs.z);
			if (xint > 0 && xint < w-1 && yint > 0 && yint < h-1 && zint > 0 && zint < d-1)
			{
				spatindex = (xint + w*(yint+h*zint)) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2]);
			}
		}
		
		Dn /= (REAL)(2*CU+1); 
		return Dn;
	}
*/


    inline REAL evaluateODF(pVector &R, pVector &N, REAL &len)
	{
		const int CU = 10;
		pVector Rs;
		REAL Dn = 0;
		int xint,yint,zint,spatindex;

		sinterp->getInterpolation(N);

		for (int i=-CU; i <= CU;i++)
		{
			Rs = R + (N * len) * ((REAL)i/CU);
        
            Rs.storeXYZ();
            REAL Rx = pVector::store[0]/voxsize_w-0.5;
            REAL Ry = pVector::store[1]/voxsize_h-0.5;
            REAL Rz = pVector::store[2]/voxsize_d-0.5;


            xint = int(floor(Rx));
            yint = int(floor(Ry));
            zint = int(floor(Rz));
	
			
			if (xint >= 0 && xint < w-1 && yint >= 0 && yint < h-1 && zint >= 0 && zint < d-1)
			{
				REAL xfrac = Rx-xint;
				REAL yfrac = Ry-yint;
				REAL zfrac = Rz-zint;
			
				REAL weight;
							
				weight = (1-xfrac)*(1-yfrac)*(1-zfrac);
				spatindex = (xint + w*(yint+h*zint)) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
	
				weight = (xfrac)*(1-yfrac)*(1-zfrac);
				spatindex = (xint+1 + w*(yint+h*zint)) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (1-xfrac)*(yfrac)*(1-zfrac);
				spatindex = (xint + w*(yint+1+h*zint)) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (1-xfrac)*(1-yfrac)*(zfrac);
				spatindex = (xint + w*(yint+h*(zint+1))) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (xfrac)*(yfrac)*(1-zfrac);
				spatindex = (xint+1 + w*(yint+1+h*zint)) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (1-xfrac)*(yfrac)*(zfrac);
				spatindex = (xint + w*(yint+1+h*(zint+1))) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (xfrac)*(1-yfrac)*(zfrac);
				spatindex = (xint+1 + w*(yint+h*(zint+1))) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
				weight = (xfrac)*(yfrac)*(zfrac);
				spatindex = (xint+1 + w*(yint+1+h*(zint+1))) *nip;
				Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
			
					
			}


		}
		
		Dn *= 1.0/(2*CU+1); 
		return Dn;
	}


/*
	inline REAL evaluateODF(pVector &R, pVector &N, REAL &len)
	{
                
		R.storeXYZ();
        
		REAL Rx = pVector::store[0]/voxsize_w;
		REAL Ry = pVector::store[1]/voxsize_h;
		REAL Rz = pVector::store[2]/voxsize_d;


		int xint = int(Rx);
		int yint = int(Ry);
		int zint = int(Rz);
		
		if (xint >= 0 && xint < w-1 && yint >= 0 && yint < h-1 && zint >= 0 && zint < d-1)
		{
			REAL xfrac = Rx-xint;
			REAL yfrac = Ry-yint;
			REAL zfrac = Rz-zint;
			sinterp->getInterpolation(N);
		
			REAL weight;
			int spatindex;
			REAL Dn = 0;
		
			weight = (1-xfrac)*(1-yfrac)*(1-zfrac);
			spatindex = (xint + w*(yint+h*zint)) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;

			weight = (xfrac)*(1-yfrac)*(1-zfrac);
			spatindex = (xint+1 + w*(yint+h*zint)) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (1-xfrac)*(yfrac)*(1-zfrac);
			spatindex = (xint + w*(yint+1+h*zint)) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (1-xfrac)*(1-yfrac)*(zfrac);
			spatindex = (xint + w*(yint+h*(zint+1))) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (xfrac)*(yfrac)*(1-zfrac);
			spatindex = (xint+1 + w*(yint+1+h*zint)) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (1-xfrac)*(yfrac)*(zfrac);
			spatindex = (xint + w*(yint+1+h*(zint+1))) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (xfrac)*(1-yfrac)*(zfrac);
			spatindex = (xint+1 + w*(yint+h*(zint+1))) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			weight = (xfrac)*(yfrac)*(zfrac);
			spatindex = (xint+1 + w*(yint+1+h*(zint+1))) *nip;
			Dn += (dataimg[spatindex + sinterp->idx[0]]*sinterp->interpw[0] + dataimg[spatindex + sinterp->idx[1]]*sinterp->interpw[1] + dataimg[spatindex + sinterp->idx[2]]* sinterp->interpw[2])*weight;
		
			return Dn;
				
		}
		return 0;
	}

*/

	virtual inline REAL computeExternalEnergy(pVector &R, pVector &N, REAL &cap, REAL &len, Particle *dp) { return 0;}
	virtual	inline REAL computeInternalEnergy(Particle *p1) {return 0;}
	virtual	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1) {return 0;}
	virtual	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1, Particle *p2, int ep2) {return 0;}
        


};

#endif





#ifndef _ENCOMPINTERFACE
#define _ENCOMPINTERFACE

             

#include "SphereInterpolator.cpp"



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

    int w_sp,h_sp,d_sp;    // of mask
    REAL voxsize_sp_w;
    REAL voxsize_sp_h;
    REAL voxsize_sp_d;
    
    int w_vf,h_vf,d_vf;    
    
    int nip; // number of data vertices on sphere
    int nalpha; // number of expansions coeffs in LUT
    
    
	REAL *msmap; // avergae q-space signal (at native DWI resolution)
    
    REAL *vfmap; // vfmap[0] ~ vf_int  (at oversampled resolution)
                 // vfmap[w*h*d] ~ q-space average of the model
                 // vfmap[2*w*h*d] ~ the S0 of the model

    int *idxmap;
	REAL *spatprobimg;
	REAL *cumulspatprob;
	int *spatidx;
	int scnt; 
    int GM_mask; //>0 if GM_mask is present
    
    
    REAL sigma_s;
	REAL particle_weight;
    
    REAL prior_strength;
    REAL eigen_energy;
    REAL meanval_sq;
    
    REAL guide_weight;

    EnergyComputerBase(REAL *data, const int *dsz, double *voxsize, SphereInterpolator *sp, ParticleGrid<Particle> *pcon, REAL *vfmap, REAL *spimg, REAL *s2andidx, int spmult)
    {
        dataimg = data;
        datasz = dsz;
        sinterp = sp;   

        spatprobimg = spimg;
        this->vfmap = vfmap;
        
        nalpha = datasz[0];
        nip = datasz[1];
                
        w = datasz[2];
        h = datasz[3];
        d = datasz[4];

        voxsize_w = voxsize[0];
        voxsize_h = voxsize[1];
        voxsize_d = voxsize[2];
               
        w_sp = w*spmult;
        h_sp = h*spmult;
        d_sp = d*spmult;
        
        voxsize_sp_w = voxsize[0]/spmult;
        voxsize_sp_h = voxsize[1]/spmult;
        voxsize_sp_d = voxsize[2]/spmult;
        
        msmap = s2andidx;

        idxmap = (int*) malloc(sizeof(int)*w*h*d);
        for (int k = 0; k < w*h*d; k++)
            idxmap[k] = int(s2andidx[k+w*h*d])-1;
        
        
        fprintf(stderr,"Data size (voxels) : %i x %i x %i\n",w,h,d);
		fprintf(stderr,"Data directions (sizeLUT * dirs) : %i * %i\n",nalpha,nip);
        
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
        free(idxmap);
    }

	virtual void setParameters()
	{
			mexPrintf("in set base\n");
	}
	

    // defined in interactionLUTs.h
    REAL modelmodelCorr(REAL Dpara1,REAL Dorth1,REAL Dpara2,REAL Dorth2,REAL dot);
    REAL modelmodelCorrFixTis(REAL dot);
    REAL modelSignalCorr(REAL Dpara,REAL Dorth,REAL *alpha);
    REAL modelSignalCorrGuide(REAL *alpha);
    
    
    
    
    //-------------------- computes model-model correlation function
    inline REAL mminteract(REAL &q1,REAL &w1,REAL &Di1, REAL &Da1, REAL &Dp1, REAL fi1, REAL fw1,
                           REAL &q2,REAL &w2,REAL &Di2, REAL &Da2, REAL &Dp2, REAL fi2, REAL fw2, REAL dot,
                           REAL &res_sw_proj1,REAL &res_sw_proj2)
    {    
        
        
            REAL res = 0;
            
            // needed to keep track of positivity of sw fraction
            res_sw_proj1 = 0;
            res_sw_proj2 = 0;
            
    
            // modelmodel correlation
            REAL Mi1 = modelmodelCorr(Di1,0.0,0,0.0,dot);
            REAL Mi2 = modelmodelCorr(Di2,0.0,0,0.0,dot);
            REAL Ma1 = modelmodelCorr(Da1,Dp1,0,0,dot);
            REAL Ma2 = modelmodelCorr(Da2,Dp2,0,0,dot);            
            res = fi1*fi2*(modelmodelCorr(Di1,0.0,Di2,0.0,dot)- meanval_sq*Mi1*Mi2) + (1-fi1)*(1-fi2)*(modelmodelCorr(Da1,Dp1,Da2,Dp2,dot)-meanval_sq*Ma1*Ma2) +                      
                    fi1*(1-fi2)*(modelmodelCorr(Di1,0.0,Da2,Dp2,dot)-meanval_sq*Mi1*Ma2) + (1-fi1)*fi2*(modelmodelCorr(Da1,Dp1,Di2,0.0,dot)-meanval_sq*Ma1*Mi2) ;
            res_sw_proj1 = w1*(fi1*Mi1 + (1-fi1)*Ma1);
            res_sw_proj2 = w2*(fi2*Mi2 + (1-fi2)*Ma2);                                
            res *= (w1*w2);
            
            // tracking guide prior
            res += q1*q2*modelmodelCorrFixTis(dot);
                                
            return res;
                    
    }

    

    //------------------------- returns the model-signal correlation for a given particle
    inline REAL getModelSignalCorrelation(Particle *p, REAL volfrac_int)
	{
        
        REAL *alpha = (REAL*) malloc(sizeof(REAL)*nalpha); // windoof
		int xint,yint,zint,spatindex;

        // get voxel index
        REAL Rx = p->R.x/voxsize_w;
        REAL Ry = p->R.y/voxsize_h;
        REAL Rz = p->R.z/voxsize_d;        
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));
                

        // get barycentric interpolation weights/indices
        int i0,i1,i2;
        REAL w0,w1,w2;
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (SPHEREINTERPOL)
        #endif
        {
            sinterp->getInterpolation(p->N);
            i0 = sinterp->idx[0]*nalpha;
            i1 = sinterp->idx[1]*nalpha;
            i2 = sinterp->idx[2]*nalpha;
            w0 = sinterp->interpw[0];
            w1 = sinterp->interpw[1];
            w2 = sinterp->interpw[2];
        }
        
        
        
        // gather alpha vector needed for modelSignalCorr
        REAL weight;        
        for (int k = 0; k < nalpha; k++)            
        {
                alpha[k] = 0;
                
                if (xint >= 0 && xint < w && yint >= 0 && yint < h && zint >= 0 && zint < d)
                {
                    
                    int sindex = idxmap[(xint + w*(yint+h*zint))];
                    if (sindex >= 0)
                    {
                        spatindex =  sindex*nip*nalpha;                    
                        alpha[k] += (dataimg[spatindex+i0+k]*w0 + dataimg[spatindex+i1+k]*w1 + dataimg[spatindex+i2+k]* w2);
                    }
                }

		}
        
        REAL res = 0;
        
        // model signal correaltion
        res = p->w*(volfrac_int*modelSignalCorr(p->Di, 0.0, alpha) + (1-volfrac_int)*modelSignalCorr(p->Da, p->Dp, alpha));                
                
        // tracking guide prior
        res += p->q*modelSignalCorrGuide(alpha);

        free(alpha);
        
        return res;
	}
    
    
 

    
    //----------------- gets the q-space average at some position R
    inline REAL getMS(pVector &R)
	{
		int xint,yint,zint,spatindex;
        
        REAL Rx = R.x/voxsize_w;
        REAL Ry = R.y/voxsize_h;
        REAL Rz = R.z/voxsize_d;
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));
	
        REAL xfrac = Rx-xint;
        REAL yfrac = Ry-yint;
        REAL zfrac = Rz-zint;
                
        REAL weight;                
        REAL res = 0;
        
        if (xint >= 0 && xint < w && yint >= 0 && yint < h && zint >= 0 && zint < d)
        {

            spatindex = (xint + w*(yint+h*zint));                    
            res += msmap[spatindex];
        }

        return res;
        
	}

    
    //----------------- gets the intenral volume fraction at some position R
    inline void getVF_int(pVector &R, REAL &vfi)
	{
        
                
		int xint,yint,zint,spatindex;
        
        REAL Rx = R.x/(voxsize_w*sigma_s);
        REAL Ry = R.y/(voxsize_h*sigma_s);
        REAL Rz = R.z/(voxsize_d*sigma_s);
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));
	                
                
        vfi = 0;        
        if (xint >= 0 && xint < w_vf && yint >= 0 && yint < h_vf && zint >= 0 && zint < d_vf)
        {
            spatindex = (xint + w_vf*(yint+h_vf*zint));                    
            vfi = vfmap[spatindex];            
        }
        
	}
    
 
    
    //------------------- returns some random position inside WM-mask
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
	
	
    //---------------- returns the probability that drawSpatPosition draws this position
	REAL SpatProb(pVector R)
	{
			
		int rx = int(floor(R.x/voxsize_sp_w));
		int ry = int(floor(R.y/voxsize_sp_h));
		int rz = int(floor(R.z/voxsize_sp_d));
		if (rx >= 0 && rx < w_sp && ry >= 0 && ry < h_sp && rz >= 0 && rz < d_sp)
			return spatprobimg[rx + w_sp* (ry + h_sp*rz)];
		else
			return 0;	
	}

    
	int isGrayMatter(pVector R)
	{
        if (GM_mask>0)
        {	
            int rx = int(floor(R.x/voxsize_sp_w));
            int ry = int(floor(R.y/voxsize_sp_h));
            int rz = int(floor(R.z/voxsize_sp_d));
            if (rx >= 0 && rx < w_sp && ry >= 0 && ry < h_sp && rz >= 0 && rz < d_sp)
                return spatprobimg[w_sp*h_sp*d_sp + rx + w_sp* (ry + h_sp*rz)]>0;
            else
                return 0;	
        }
        else
            return (SpatProb(R) == 0);
	}
    
    
    
    inline REAL getPriorDirProb_t(pVector R,pVector N)
	{
        
		int xint,yint,zint,spatindex;

        // get voxel index
        REAL Rx = R.x/voxsize_w;
        REAL Ry = R.y/voxsize_h;
        REAL Rz = R.z/voxsize_d;        
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));

        
         // get barycentric interpolation weights/indices
        int i0,i1,i2;
        REAL w0,w1,w2;
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (SPHEREINTERPOL)
        #endif
        {
            sinterp->getInterpolation(N);
            i0 = sinterp->idx[0]*nalpha;
            i1 = sinterp->idx[1]*nalpha;
            i2 = sinterp->idx[2]*nalpha;
            w0 = sinterp->interpw[0];
            w1 = sinterp->interpw[1];
            w2 = sinterp->interpw[2];
        }
        
        
        
        int k = nalpha-1; 
        if (xint >= 0 && xint < w && yint >= 0 && yint < h && zint >= 0 && zint < d)
        {

            int sindex = idxmap[(xint + w*(yint+h*zint))];
            if (sindex >= 0)
            {                  
                return (dataimg[spatindex+i0+k]*w0 + dataimg[spatindex+i1+k]*w1 + dataimg[spatindex+i2+k]* w2);                
            }
            else
                return 1;
        }
        else
            return 1;
                
        
        
    }
    
    inline REAL getPriorDirProb(pVector R,pVector N)
    {
        //return 1;
        return 0.5*(getPriorDirProb_t( R, N)+getPriorDirProb_t( R, N*(-1)));
    }

    
    
    inline pVector drawPriorDir(pVector R)
	{
        
        
		int xint,yint,zint,spatindex;

        // get voxel index
        REAL Rx = R.x/voxsize_w;
        REAL Ry = R.y/voxsize_h;
        REAL Rz = R.z/voxsize_d;        
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));

        pVector dummy;
        dummy.rand_sphere();        
        
      //  return dummy;
        int k = nalpha-1; 
        if (xint >= 0 && xint < w && yint >= 0 && yint < h && zint >= 0 && zint < d)
        {

            int sindex = idxmap[(xint + w*(yint+h*zint))];
            if (sindex >= 0)
            {
                  
                spatindex =  sindex*nip*nalpha;                    
                SimpSamp<int> *RS = new SimpSamp<int>(nip);
                
                for (int j = 0; j < nip;j++)                    
                    RS->add(dataimg[spatindex+j*nalpha+k],j);
                int d = RS->drawObj();
                pVector N(sinterp->dirs[d*3],sinterp->dirs[d*3+1],sinterp->dirs[d*3+2]);      
                delete RS;
                return N;
            }
            else
                return dummy;
        }
        else
            return dummy;
        
        
        
    }
    
    
    //---------------- computes 1. the q-space average of the model and 2. the S0 of the model
    //---------------- and stores it for computing the vfsw (in matlab)
    void computeMeanProjMap()
    {
        for (int z = 0; z < d_vf;z++)
        for (int y = 0; y < h_vf;y++)
        for (int x = 0; x < w_vf;x++)
        {          
            int spatindex = (x + w_vf*(y+h_vf*z));                    
            vfmap[spatindex+w_vf*h_vf*d_vf] = 0;
            vfmap[spatindex+w_vf*h_vf*d_vf*2] = 0;
            
            int cnt1; 
            Particle **P1 = pcontainer->getCell(x,y,z,cnt1);
            for (int k = 0; k < cnt1; k++)
            {
                Particle *p1 = P1[k];
                vfmap[spatindex+w_vf*h_vf*d_vf] += particle_weight*p1->w*(modelmodelCorr(p1->Da,p1->Dp,0,0,1)*(1-vfmap[spatindex]) + modelmodelCorr(p1->Di,0,0,0,1)*vfmap[spatindex]);
                vfmap[spatindex+2*w_vf*h_vf*d_vf] += particle_weight*p1->w;               

            }            
        }
            
    }
    
        
	virtual inline REAL computeExternalEnergy(Particle *p,Particle *dp) {return 0;}
	virtual inline REAL computeExternalEnergy(pVector &R) {return 0;}	
	virtual inline REAL computeInternalEnergy(pVector &R) {return 0;}	    
    virtual	inline REAL computeInternalEnergy(Particle *p1) {return 0;}
	virtual	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1) {return 0;}
	virtual	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1, Particle *p2, int ep2) {return 0;}
        
    virtual inline REAL smoothnessVF(pVector &R) { return 0; }
    virtual inline REAL smoothnessSegs(Particle *p) {return 0;}
};

#include <interactionLUTs.h>


#endif



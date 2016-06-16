

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
    int lmax; // nalpha = lmax*numbshells
    int numbshells; // nalpha = lmax*numbshells
    
    
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
        
        numbshells = 0;
        lmax =  0;
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
		fprintf(stderr,"Data directions (sizeLUT * dirs) : %i,%i * %i\n",lmax-2, numbshells,nip);
        
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
    inline void mminteract_modelfree(REAL dot, REAL fac, REAL *alpha)
    {    

        dot = fabs(dot);
        
        if (dot==1)
        {
           for (int k = 0; k < lmax; k++)
                alpha[k] += fac;
           return;
        }
        
        alpha[0] += fac;
                
        REAL tmp1 = 1;
        REAL tmp2 = dot;
        REAL tmp3,tmp4;
        
        for (int k = 0; k < lmax-1; k++)
        {
            int l = 2*k+2;
            tmp3 = ((2*l-1)*dot*tmp2 - (l-1)*tmp1)/l;
            tmp4 = ((2*l+1)*dot*tmp3 - (l)*tmp2)/(l+1);
            alpha[k+1] += fac*tmp3;
            tmp1 = tmp3;
            tmp2 = tmp4;
            
        }
                
                                    
    }

    
    
    
  inline void getModelSignalCorrelation_modelfree(Particle *p, REAL *alpha)
	{
        
		int xint,yint,zint,spatindex;

        // get voxel index
        REAL Rx = p->R.x/voxsize_w;
        REAL Ry = p->R.y/voxsize_h;
        REAL Rz = p->R.z/voxsize_d;        
        xint = int(floor(Rx));
        yint = int(floor(Ry));
        zint = int(floor(Rz));
                

        int offs = (lmax-2)*numbshells;
        
        // get barycentric interpolation weights/indices
        int i0,i1,i2;
        REAL w0,w1,w2;
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (SPHEREINTERPOL)
        #endif
        {
            sinterp->getInterpolation(p->N);
            i0 = sinterp->idx[0]*offs;
            i1 = sinterp->idx[1]*offs;
            i2 = sinterp->idx[2]*offs;
            w0 = sinterp->interpw[0];
            w1 = sinterp->interpw[1];
            w2 = sinterp->interpw[2];
        }
        
       
        // gather alpha vector needed for modelSignalCorr
        REAL weight;        
        const REAL sq23 = 1/(sqrt(2)*3);
                
        if (xint >= 0 && xint < w && yint >= 0 && yint < h && zint >= 0 && zint < d)
        {

            int sindex = idxmap[(xint + w*(yint+h*zint))];
            if (sindex >= 0)
            {
                spatindex =  sindex*(nip*offs + numbshells*7);      
                REAL tmp ;
                
                for(int k = 0; k< numbshells;k++)
                    for(int j = 2; j < lmax;j++)
                    {
                        int idx = j+k*lmax;
                        int idx2 = j-2+k*(lmax-2);
//                         int idx = j+k*numbshells;
//                         int idx2 = j+(k-2)*numbshells;
                        alpha[idx] += p->w*(dataimg[spatindex+i0+idx2]*w0 + dataimg[spatindex+i1+idx2]*w1 + dataimg[spatindex+i2+idx2]* w2);
                    }
                
//                 for (int k = 2; k < nalpha; k++)
//                 {
//                     if (k%lmax > 1)
//                     {
//                     alpha[k-2] += p->w*(dataimg[spatindex+i0+k]*w0 + dataimg[spatindex+i1+k]*w1 + dataimg[spatindex+i2+k]* w2);
//                     }
//                 }
                for(int k = 0; k< numbshells;k++)
                {
                    pVector n = p->N;
                    REAL qxx = (2*n.x*n.x - n.y*n.y - n.z*n.z)*sq23;
                    REAL qyy = (2*n.y*n.y - n.x*n.x - n.z*n.z)*sq23;
                    REAL qzz = (2*n.z*n.z - n.y*n.y - n.x*n.x)*sq23;
                    REAL qxy = n.x*n.y;
                    REAL qxz = n.x*n.z;
                    REAL qyz = n.z*n.y;
                    REAL d2 = p->w*3*(
                     qxx*dataimg[spatindex+nip*offs +1 + k*7 ] +
                     qyy*dataimg[spatindex+nip*offs +2 + k*7] +
                     qzz*dataimg[spatindex+nip*offs +3 + k*7] +
                     qxy*dataimg[spatindex+nip*offs +4 + k*7] +
                     qxz*dataimg[spatindex+nip*offs +5 + k*7] +
                     qyz*dataimg[spatindex+nip*offs +6 + k*7]);
                    
//                      if (mtrand.frand() > 0.9999)
//                      {
//                         if (k == 0)
//                          fprintf(stderr,"%i) %.20f  \n",k,tmp/d2);
//                       //   fprintf(stderr,"%i)-- %f  \n",k,dataimg[spatindex+nip*nalpha ]/60/alpha[0+k*lmax]);
//                      }
//                     
                    alpha[0+k*lmax] += dataimg[spatindex+nip*offs +0 + k*7 ]*p->w;
                    alpha[1+k*lmax] += d2;
                }
                    
            
            
            }
        }
        
        
        
        
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
            return false; //(SpatProb(R) == 0);
	}
    
    
    
    
    //---------------- computes 1. the q-space average of the model and 2. the S0 of the model
    //---------------- and stores it for computing the vfsw (in matlab)
    void computeMeanProjMap(REAL *vf_fibs)
    {
        
        
        
          
        REAL *alphaSM = (REAL*) malloc(sizeof(REAL)*lmax*numbshells); 
        REAL *alphaMM = (REAL*) malloc(sizeof(REAL)*lmax); 
        for (int k = 0; k < lmax; k++)
            alphaMM[k] = 0;
        for (int k = 0; k < lmax*numbshells; k++)
            alphaSM[k] = 0;

   
        
        for (int z = 0; z < d_vf;z++)
        for (int y = 0; y < h_vf;y++)
        for (int x = 0; x < w_vf;x++)
        {          
            int spatindex = (x + w_vf*(y+h_vf*z));                    
          
            
            for (int k = 0; k < lmax; k++)
                alphaMM[k] = 0;
            for (int k = 0; k < lmax*numbshells; k++)
                alphaSM[k] = 0;
                        
            
            int cnt1; 
            Particle **P1 = pcontainer->getCell(x,y,z,cnt1);
            REAL sumW = 0;
            for (int k = 0; k < cnt1; k++)
            {
                Particle *p1 = P1[k];            
                sumW += p1->w;
                getModelSignalCorrelation_modelfree(p1, alphaSM)   ;         
                for (int j = k; j < cnt1; j++)
                {
                    Particle *p2 = P1[j];
                    REAL fac = ((k!=j)?2:1) * p1->w*p2->w;
                    mminteract_modelfree(p2->N*p1->N, fac,alphaMM );
                }                                            
            }

            if (vf_fibs != 0)
            {
                for (int k = 0; k < cnt1; k++)
                {
                    Particle *p1 = P1[k];            
                    int seg_idx = pcontainer->ID_2_index(p1->ID);

                    for (int k = 0; k < numbshells;k++)
                        for (int j = 0; j < lmax;j++)
                        {
                            int idx = j+(k+1)*lmax;
                            vf_fibs[seg_idx+pcontainer->pcnt*idx] = alphaSM[j+k*lmax]/sumW;                                        
                        }

                    for (int j = 0; j < lmax;j++)
                    {
                        vf_fibs[seg_idx+pcontainer->pcnt*j] = alphaMM[j]/(sumW*sumW);                                        
                    }
                }            
            }
            
            for (int k = 0; k < numbshells;k++)
                for (int j = 0; j < lmax;j++)
                {
                    int idx = j+(k+1)*lmax;
                    vfmap[spatindex+w_vf*h_vf*d_vf*idx] = alphaSM[j+k*lmax]/sumW;                                        
                }
            
            for (int j = 0; j < lmax;j++)
            {
                vfmap[spatindex+w_vf*h_vf*d_vf*j] = alphaMM[j]/(sumW*sumW);                                        
            }
        }

        
        free(alphaMM);
        free(alphaSM);
        
            
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

//#include <interactionLUTs.h>


#endif



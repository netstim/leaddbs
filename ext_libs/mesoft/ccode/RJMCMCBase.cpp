
#include "ParticleGrid.cpp"
#include <mex.h>
#include <vector>

#ifdef PARALLEL_OPENMP
#include <omp.h>
#endif

class RJMCMCBase
{
	public:


	ParticleGrid<Particle> pcontainer;
    REAL *vfmap;
    
	REAL *dataimg;
	const int *datasz;


  	EnergyComputerBase *enc;
    

	REAL itmax;
    REAL gammaDE;
	REAL width;
	REAL height;
	REAL depth;
	double *voxsize;

	int attrcnt;
    int accepted;
    int iterations_done;
    int numcores;
                
	REAL p_birth;
	REAL p_death;	
	REAL p_shift;	
    REAL p_Dmod;
    REAL p_vfmod;
    REAL p_conprob;

    
    
	RJMCMCBase(REAL *points,int numPoints, REAL *vfmap, REAL *dimg, const int *dsz, double *voxsize, double cellsize, REAL *cellsize2, REAL len, int numcores)
	{

		dataimg = dimg;
		datasz = dsz;
        this->vfmap = vfmap;
        
		width = datasz[2]*voxsize[0];
		height = datasz[3]*voxsize[1];   
		depth = datasz[4]*voxsize[2];

		this->voxsize = voxsize;
        this->numcores = numcores;

		fprintf(stderr,"Data dimensions (mm) :  %f x %f x %f\n",width,height,depth);
		fprintf(stderr,"Data dimensions (voxel) :  %i x %i x %i\n",datasz[2],datasz[3],datasz[4]);
		fprintf(stderr,"Data directions (sizeLUT * dirs) : %i * %i\n",datasz[0],datasz[1]);
		fprintf(stderr,"voxel size (mm) :  %f x %f x %f\n",voxsize[0],voxsize[1],voxsize[2]);
	

        ////// dimensions of the particle grid
        // the one to find connection partners
		REAL cellcnt_x = (int)((REAL)width/cellsize) +2;
		REAL cellcnt_y = (int)((REAL)height/cellsize) +2;
		REAL cellcnt_z = (int)((REAL)depth/cellsize) +2;
		int cell_capacity = 2048;
        // for finding partners whithn the same voxel
		REAL cellcnt_x2 = (int)((REAL)width/cellsize2[0]) +2;
		REAL cellcnt_y2 = (int)((REAL)height/cellsize2[1]) +2;
		REAL cellcnt_z2 = (int)((REAL)depth/cellsize2[2]) +2;
		int cell_capacity2 = 20; 

		fprintf(stderr,"grid dimensions :  %f x %f x %f\n",cellcnt_x,cellcnt_y,cellcnt_z);
		fprintf(stderr,"grid cell size (mm) :  %f^3\n",cellsize);
		fprintf(stderr,"grid2 dimensions :  %f x %f x %f\n",cellcnt_x2,cellcnt_y2,cellcnt_z2);
		fprintf(stderr,"grid2 cell size (mm) :  %fx%fx%f\n",cellsize2[0],cellsize2[1],cellsize2[2]);
    	fprintf(stderr,"cell capacity :  %i\n",cell_capacity);
    	fprintf(stderr,"#cells*cellcap :  %.1f K\n",cell_capacity*cellcnt_x*cellcnt_y*cellcnt_z/1000);
        #ifdef PARALLEL_OPENMP
		fprintf(stderr,"#threads :  %i/%i\n",numcores,omp_get_max_threads());
        #endif
        
        // allocate the particle grid
		int minsize = 100000;
		int err = pcontainer.allocate(((numPoints>minsize)? (numPoints+100000) : minsize), 
                cellcnt_x, cellcnt_y, cellcnt_z, cellsize,cell_capacity,
                cellcnt_x2, cellcnt_y2, cellcnt_z2, cellsize2,cell_capacity2,len);
        if (err == -1)
        {
            fprintf(stderr,"RJMCMCBase: out of Memory!\n");
            return;
        }
        
        // read the data from MATLAB into our structure
		attrcnt = 15;
		for (int k = 0; k < numPoints; k++)
			{
				Particle *p = pcontainer.newParticle(pVector(points[attrcnt*k  ], points[attrcnt*k+1],points[attrcnt*k+2]),
                                                     pVector(points[attrcnt*k+3], points[attrcnt*k+4],points[attrcnt*k+5]));
				if (p!=0) 
				{
					p->N = pVector(points[attrcnt*k+3],points[attrcnt*k+4],points[attrcnt*k+5]);				
					p->len =  points[attrcnt*k+7];	
					p->mID = (int) points[attrcnt*k+8];				
					p->pID = (int) points[attrcnt*k+9];	
					p->Di =  points[attrcnt*k+10];	
					p->Da =  points[attrcnt*k+11];	
					p->Dp =  points[attrcnt*k+12];	
					p->w =  points[attrcnt*k+13];	
					p->q =  points[attrcnt*k+14];	
                    p->label = 0;
                    p->active = true;
                    
					if (p->mID != -1)
                    {
                        pcontainer.removeFromGrid(1,p);
						pcontainer.concnt++;				
                    }
					if (p->pID != -1)
                    {
                        pcontainer.removeFromGrid(0,p);
						pcontainer.concnt++;		
                    }					
				}
				else
				{
					fprintf(stderr,"error: cannot allocate particle,  con. indices will be wrong! \n");
				}
			}
		pcontainer.concnt /= 2;
   


		itmax = 0;
        accepted = 0;

	}	

	~RJMCMCBase()
	{

	}



    // write out to MATLAB
	void writeout(REAL *npoints)
	{
	
		for (int k = 0; k < pcontainer.pcnt; k++)
		{
			Particle *p = &(pcontainer.particles[k]);
			npoints[attrcnt*k] = p->R.x;
			npoints[attrcnt*k+1] = p->R.y;
			npoints[attrcnt*k+2] = p->R.z;
			npoints[attrcnt*k+3] = p->N.x;
			npoints[attrcnt*k+4] = p->N.y;
			npoints[attrcnt*k+5] = p->N.z;
			npoints[attrcnt*k+6] = p->label;	
			npoints[attrcnt*k+7] = p->len;		
			npoints[attrcnt*k+8] = pcontainer.ID_2_index(p->mID);		
			npoints[attrcnt*k+9] = pcontainer.ID_2_index(p->pID);		
			npoints[attrcnt*k+10] = p->Di;
			npoints[attrcnt*k+11] = p->Da;
			npoints[attrcnt*k+12] = p->Dp;
			npoints[attrcnt*k+13] = p->w;
			npoints[attrcnt*k+14] = p->q;                    
                    

		}
	

	}






	void setEnergyComputer(EnergyComputerBase *e)
	{
		enc = e;
	}


	void setParameters(REAL Temp)
	{

	}	
    
    
    /////////////////////////////// the MAIN LOOP
	void iterate(const mxArray *Handle)
	{
		#ifdef TIMING
		tic(&total_time);
		#endif
                
        int frq = 500000;  // number iteration between MATLAB calls of 'drawnow'
        int numits;
        int outer = 0;
        REAL updownratio_break;
        
        bool adaptiveschedule = false;
        numits = (int) itmax;

        if (frq > numits)
        {
            frq = numits;
            outer = 1;
        }
        else
        {
            outer = numits/frq;
        }
        
        
        if (gammaDE != 0)  // then we have an adative temp-schedule
        {
            REAL itmax_d = gammaDE;
            stats.lam_deltaE = itmax_d;
            birthstats.lam_deltaE = itmax_d;
            deathstats.lam_deltaE = itmax_d;
            shiftstats.lam_deltaE = itmax_d;
            capstats.lam_deltaE = itmax_d;
            vfstats.lam_deltaE = itmax_d;
            connstats.lam_deltaE = itmax_d;
            
            adaptiveschedule = true;
        }
        else
        {
            adaptiveschedule = false;
        }
        
             
        #ifdef PARALLEL_OPENMP
        omp_set_dynamic(0);     // Explicitly disable dynamic teams
        omp_set_num_threads(numcores); // Use n threads for all consecutive parallel regions
        #endif

        #ifdef _WIN64
            static struct timeval2 timeS;
        #else
            static struct timeval timeS;
        #endif
        
        long time =0;
        
               
        
        int f;
		for (f = 0; f < outer;f++)
        {
            
            stats.clear();
            birthstats.clear();
            deathstats.clear();
            connstats.clear();
            shiftstats.clear();
            capstats.clear();
            vfstats.clear();
            //fprintf(stderr,"parain\n"); fflush(stderr);
            time = 0;
            gettimeofday( &timeS, NULL);
            time -= (timeS.tv_sec*1000000 + timeS.tv_usec);
            
            #ifdef PARALLEL_OPENMP
            #pragma omp parallel for 
            #endif
            for (int it = 0; it < frq;it++)
            {
                if (!pcontainer.needs_reallocation)
                    iterate_onestep();
            }            
            
            if (pcontainer.needs_reallocation)
            {
                if (pcontainer.reallocate() == -1)
                {
                    fprintf(stderr,"out of Memory!!\n");
                    return;
                }
            }
            
            gettimeofday( &timeS, NULL);
            time += (timeS.tv_sec*1000000 + timeS.tv_usec);	
   
            fprintf(stderr,"#freeslots: %i  pcnt:%i  #particles:%i\n",
                    (int) pcontainer.freeslots.size(),(int) pcontainer.pcnt,(int)(pcontainer.pcnt-pcontainer.freeslots.size()));
            pcontainer.defrag();

            fprintf(stderr,"#it: %ik  #p: %i #c: %i time: %.2fs  \n",  frq*(f+1)/1000, pcontainer.pcnt,pcontainer.concnt,double(time)/1000000);
            stats.report_movetype(" total:");
            
            char buf[256];
            snprintf(buf, sizeof(buf), "birth (%.3f):",p_birth);
            birthstats.report_movetype(buf);
            snprintf(buf, sizeof(buf), "death (%.3f):",p_death);
            deathstats.report_movetype(buf);
            snprintf(buf, sizeof(buf), "shift (%.3f):",p_shift);
            shiftstats.report_movetype(buf);
            snprintf(buf, sizeof(buf), "model (%.3f):",p_Dmod);
            capstats.report_movetype(buf);
            snprintf(buf, sizeof(buf), "volfr (%.3f):",p_vfmod);
            vfstats.report_movetype(buf);
            snprintf(buf, sizeof(buf), "track (%.3f):",p_conprob);
            connstats.report_movetype(buf);

            
            REAL cr = (REAL) pcontainer.celloverflows;
            if (cr > 0.01)
                fprintf(stderr,"warning: celloverflows : %.2f \n",cr);
            if (*mxGetPr(Handle) != 0)
            {
                mexEvalString("drawnow");
                if (mxGetProperty(Handle,0,"Tag") == 0)
                {
                    fprintf(stderr,"Termination by User.\n\n");
                	break;
                }
            }
            accepted = 0;

            
            if (adaptiveschedule && stats.deltaE < 0)
            {
                fprintf(stderr,"equilibrium threshold reached, breaking\n");
                break;
            }
            

        }
            
        iterations_done = (f+1)*frq; //it;

	}

    virtual void iterate_onestep()
	{

	}
    
    
    
    
    int label_fibers()
	{
		
		int cur_label = 1;

		for (int k = 0; k < pcontainer.pcnt;k++)
		{

			Particle *dp =  &(pcontainer.particles[k]);
            if (dp->active)
                if (dp->label == 0)
                {
                    dp->label = cur_label;
                    labelrecursivly(dp,0);
                    cur_label++;	
                }	
		}
		return cur_label-1;

	}

	void labelrecursivly(Particle *dp,int depth)
	{
	
		int label = dp->label;

		if (dp->mID != -1)
		{
			if (pcontainer.particles[dp->mID].label == 0)
			{	
				pcontainer.particles[dp->mID].label = label;
				labelrecursivly(&(pcontainer.particles[dp->mID]),depth+1);
			}
		}
		if (dp->pID != -1)
		{
			if (pcontainer.particles[dp->pID].label == 0)
			{
				pcontainer.particles[dp->pID].label = label;
				labelrecursivly(&(pcontainer.particles[dp->pID]),depth+1);
			}
		}

	}

    
    
    struct Cross { int i; int j; REAL wi; REAL wj; };
    vector<Cross> createCrossingList()
    {

        vector<Cross> crosslist;
        ParticleGrid<Particle>::GridArray *g = pcontainer.getVoxelGridArray();
        for (int i = 0;i < g->gridsize;i++)
        {
               if (g->occnt[i] > 0)
               {
                   for (int k = 0; k < g->occnt[i]; k++)
                   {
                       Particle *p1 = g->grid[g->csize*i+ k];
                       for (int j = k+1; j < g->occnt[i]; j++)
                       {
                            Particle *p2 = g->grid[g->csize*i+ j];
                            
                            Cross x;
                            x.i = p1->label;
                            x.j = p2->label;
                            x.wi = p1->w;
                            x.wj = p2->w;
                            crosslist.push_back(x);                            
                       }
                   }
               }
                
        }
        return crosslist;

            
    }
    
    
    
    
};




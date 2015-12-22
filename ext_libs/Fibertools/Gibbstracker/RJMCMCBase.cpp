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

#include "ParticleGrid.cpp"
#include <mex.h>

class RJMCMCBase
{
	public:


	ParticleGrid<Particle> pcontainer;
	REAL *dataimg;
	const int *datasz;


  	EnergyComputerBase *enc;
    

	REAL itmax;
	REAL width;
	REAL height;
	REAL depth;
	double *voxsize;

	int attrcnt;
    int accepted;
    int iterations_done;

	RJMCMCBase(REAL *points,int numPoints, REAL *dimg, const int *dsz, double *voxsize, double cellsize)
	{

		dataimg = dimg;
		datasz = dsz;
		width = datasz[1]*voxsize[0];
		height = datasz[2]*voxsize[1];   
		depth = datasz[3]*voxsize[2];

		this->voxsize = voxsize;

		fprintf(stderr,"Data dimensions (mm) :  %f x %f x %f\n",width,height,depth);
		fprintf(stderr,"Data dimensions (voxel) :  %i x %i x %i\n",datasz[1],datasz[2],datasz[3]);
		fprintf(stderr,"voxel size (mm) :  %f x %f x %f\n",voxsize[0],voxsize[1],voxsize[2]);
	

		REAL cellcnt_x = (int)((REAL)width/cellsize) +2;
		REAL cellcnt_y = (int)((REAL)height/cellsize) +2;
		REAL cellcnt_z = (int)((REAL)depth/cellsize) +2;
		int cell_capacity = 2048;

		fprintf(stderr,"grid dimensions :  %f x %f x %f\n",cellcnt_x,cellcnt_y,cellcnt_z);
		fprintf(stderr,"grid cell size (mm) :  %f^3\n",cellsize);
    	fprintf(stderr,"cell capacity :  %i\n",cell_capacity);
    	fprintf(stderr,"#cells*cellcap :  %.1f K\n",cell_capacity*cellcnt_x*cellcnt_y*cellcnt_z/1000);


		int minsize = 1000000;
		int err = pcontainer.allocate(((numPoints>minsize)? (numPoints+100000) : minsize), cellcnt_x, cellcnt_y, cellcnt_z, cellsize,cell_capacity);
//		pcontainer.allocate(1000000, cellcnt_x, cellcnt_y, cellcnt_y, cellsize,cell_capacity);


        if (err == -1)
        {
            fprintf(stderr,"RJMCMCBase: out of Memory!\n");
            return;
        }
        
		attrcnt = 10;
		for (int k = 0; k < numPoints; k++)
			{
				Particle *p = pcontainer.newParticle(pVector(points[attrcnt*k], points[attrcnt*k+1],points[attrcnt*k+2]));
				if (p!=0) 
				{
					p->N = pVector(points[attrcnt*k+3],points[attrcnt*k+4],points[attrcnt*k+5]);				
					p->cap =  points[attrcnt*k+6];		
					p->len =  points[attrcnt*k+7];	
					p->mID = (int) points[attrcnt*k+8];				
					p->pID = (int) points[attrcnt*k+9];	
					if (p->mID != -1)
						pcontainer.concnt++;				
					if (p->pID != -1)
						pcontainer.concnt++;		
					p->label = 0;		
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



	
	void writeout(REAL *npoints)
	{
	
		for (int k = 0; k < pcontainer.pcnt; k++)
		{
			Particle *p = &(pcontainer.particles[k]);
			p->R.storeXYZ();
			npoints[attrcnt*k] = pVector::store[0];
			npoints[attrcnt*k+1] = pVector::store[1];
			npoints[attrcnt*k+2] = pVector::store[2];
			p->N.storeXYZ();
			npoints[attrcnt*k+3] = pVector::store[0];
			npoints[attrcnt*k+4] = pVector::store[1];
			npoints[attrcnt*k+5] = pVector::store[2];
			npoints[attrcnt*k+6] = p->cap;		
			npoints[attrcnt*k+7] = p->len;		
			npoints[attrcnt*k+8] = pcontainer.ID_2_index(p->mID);		
			npoints[attrcnt*k+9] = pcontainer.ID_2_index(p->pID);		
		}
	

	}






	void setEnergyComputer(EnergyComputerBase *e)
	{
		enc = e;
	}


	void setParameters(REAL Temp)
	{

	}	

	void iterate(const mxArray *Handle)
	{
		#ifdef TIMING
		tic(&total_time);
		#endif
        int frq = 100000;
        int numits;
        REAL updownratio_break;
        
        bool adaptiveschedule = false;
        if (itmax <= 1)
        {
            updownratio_break = itmax;
            adaptiveschedule = true;
        }
        else
        {
            numits = (int) itmax;
            adaptiveschedule = false;
        }
             
        int it;
		for (it = 0; it < numits || adaptiveschedule;it++)
		{
			iterate_onestep();
			
			if (it%frq == frq-1)
			{
				fprintf(stderr,"#iterations: %ik   #particles: %ik   #connections: %ik   \n",
						it/1000, pcontainer.pcnt/1000,pcontainer.concnt/1000);
                stats.report_movetype(" props:");
                if (adaptiveschedule && stats.updownratio > updownratio_break)
                {
                    fprintf(stderr,"%f %f\n", stats.updownratio , updownratio_break);
                    break;
                }
                
                #ifdef TIMING
            // 	mexPrintf("\nEnergy\n------------------------\n");
            // 	externalenergy_time.report("external  ");
            // 	odfeval_time.report("odfeval   ");
            // 	internalenergy_time.report("internal  ");
            // 	total_time.report_time("total energy comp.  ");
            // 	
            // 	mexPrintf("overhead for proposals und stuff:%.1fms\n",
            // 	(total_time.time-(externalenergy_time.time+odfeval_time.time+internalenergy_time.time))/1000.0);

                mexPrintf("\nProposals\n------------------------\n");
                birthproposal_time.report("birth  ");
                deathproposal_time.report("death  ");
                shiftproposal_time.report("shift  ");
                connproposal_time.report("conne  ");
            //	capproposal_time.report("capch  ");
            //	lenproposal_time.report("length ");
                mexPrintf("\n");
                birthstats.report("birth  ");
                deathstats.report("death  ");
                shiftstats.report("shift  ");
                connstats.report("conne  ");
            //	lenstats.report("length ");
            //	capstats.report("capch  ");
                mexPrintf("\n");
                birthstats.clear();
                deathstats.clear();
                shiftstats.clear();
                connstats.clear();


                #endif                
                
                
                
                
				REAL cr = (REAL) pcontainer.celloverflows / it *100000;
				if (cr > 0.01)
					fprintf(stderr,"warning: celloverflows per 10??? it : %.2f \n",cr);
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

				
			}


		}
	
        iterations_done = it;

		#ifdef TIMING
		toc(&total_time);
		#endif

	}

        virtual void iterate_onestep()
	{

	}
};




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



#ifdef TIMING

static struct timeval timeS;

class PropStats
{
	int N;
	int accept;
	public:
	void clear() { N = 0; accept = 0;}
	void propose() {N++;}
	void accepted() {accept++;}

	void report(const char *s)
	{
		mexPrintf("%s  #proposals: %8.2fk  acceptratio: %.2f \% \n",s,1.0*N/1000.0,100.0*accept/N);
	}
};


class Timing
{
	public:
	Timing() { time = 0; ncalls = 0;}
	void clear() {time = 0; ncalls=0;}
		

	long time;
	int ncalls;

	void report(const char *s)
	{
		mexPrintf("%s total: %10.2fms calls: %10.1fk t/call: %10.3fms \n",s,time/1000.0,1.0*ncalls/1000.0,1.0*time/ncalls); 
	}

	void report_time(const char *s)
	{
		mexPrintf("%s: %.2fms \n",s,time/1000.0); 
	}

};

inline void tic(Timing *t)
{
	gettimeofday( &timeS, NULL);
	t->time -= (timeS.tv_sec*1000000 + timeS.tv_usec);
	t->ncalls++;
}
inline void toc(Timing *t)
{
	gettimeofday( &timeS, NULL);
	t->time += (timeS.tv_sec*1000000 + timeS.tv_usec);	
}


#endif



#include "MersenneTwister.h"
MTRand mtrand;




#include "ParticleGrid.cpp"


class ParticleL : public Particle
{
	public:
	int label;
	int numerator;
};


class CCAnalysis
{
	public:


		
	ParticleL *particles;
	int pcnt;

	int attrcnt;


	CCAnalysis(REAL *points,int numPoints)
	{

		particles = (ParticleL*) malloc(sizeof(ParticleL)*numPoints);
		pcnt = numPoints;
		attrcnt = 10;
		for (int k = 0; k < numPoints; k++)
			{
				ParticleL *p = &(particles[k]);
				p->R = pVector(points[attrcnt*k], points[attrcnt*k+1],points[attrcnt*k+2]);
				p->N = pVector(points[attrcnt*k+3],points[attrcnt*k+4],points[attrcnt*k+5]);				
				p->cap =  points[attrcnt*k+6];		
				p->len =  points[attrcnt*k+7];	
				p->mID = (int) points[attrcnt*k+8];				
				p->pID = (int) points[attrcnt*k+9];
				p->ID = k;
				p->label = 0;					
			}
        
		

	}	

	~CCAnalysis()
	{
		free(particles);
	}



        int iterate()
	{
		
		int cur_label = 1;

		for (int k = 0; k < pcnt;k++)
		{

			ParticleL *dp =  &(particles[k]);
			if (dp->label == 0)
			{
				dp->label = cur_label;
				labelrecursivly(dp,0);
				cur_label++;	
			}	
		}
		return cur_label-1;

	}

	void labelrecursivly(ParticleL *dp,int depth)
	{
	
		int label = dp->label;

		if (dp->mID != -1)
		{
			if (particles[dp->mID].label == 0)
			{	
				particles[dp->mID].label = label;
				labelrecursivly(&(particles[dp->mID]),depth+1);
			}
		}
		if (dp->pID != -1)
		{
			if (particles[dp->pID].label == 0)
			{
				particles[dp->pID].label = label;
				labelrecursivly(&(particles[dp->pID]),depth+1);
			}
		}

	}







};





static int cmpfloat2(const void *p1,const void *p2)
{
	if (((REAL*)p1)[1] > ((REAL*)p2)[1])
        return 1;
    else
        return -1;
}



typedef std::vector<int> vecint;



void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 2) {
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

	fprintf(stderr,"building fibers "); fflush(stderr);
	
	int pcnt = 0;
	const mxArray *Points;
	Points = prhs[pcnt++];       
	int numPoints = mxGetN(Points);
	REAL *points = (REAL*) mxGetData(Points);
	
	
	const mxArray *Params;
	Params = prhs[pcnt++];       
	double *params = (double*) mxGetPr(Params);


	int minnumelements = (int) params[0];
	int maxnumelements = (int) params[1];
	

	CCAnalysis ccana(points,numPoints);

	#ifdef TIMING
	
	#endif

	int numc = ccana.iterate();

	fprintf(stderr,"."); fflush(stderr);

	vector<vecint> components(numc);
    
    int i;

	for (i = 0;i < ccana.pcnt;i++)
	{
		components[ccana.particles[i].label-1].push_back(i);
	}

	fprintf(stderr,"."); fflush(stderr);

	for (i = 0; i < numc; i++)
	{
		ParticleL *last = &(ccana.particles[components[i][0]]);
		last->numerator = 0;
		ParticleL *next = (last->pID == -1)? 0 : &(ccana.particles[last->pID]);
		for (;;)
		{
			if (next == 0)
				break;
			next->numerator = last->numerator+1;
			int nextID = -1;
			if (next->pID != last->ID)
				nextID = next->pID;
			if (next->mID != last->ID)
				nextID = next->mID;
			last = next;
			next = (nextID == -1)? 0: &(ccana.particles[nextID]);
			if (last->numerator > components[i].size()) // circular
				break;
		}

		last = &(ccana.particles[components[i][0]]);		
		next = (last->mID == -1)? 0 : &(ccana.particles[last->mID]);
		for (;;)
		{
			if (next == 0)
				break;
			next->numerator = last->numerator-1;
			int nextID = -1;
			if (next->pID != last->ID)
				nextID = next->pID;
			if (next->mID != last->ID)
				nextID = next->mID;
			last = next;
			next = (nextID == -1)? 0: &(ccana.particles[nextID]);
			if (last->numerator < -components[i].size()) // circular
				break;
		}

	}

	fprintf(stderr,"."); fflush(stderr);


	#ifdef TIMING
	
	#endif

	int index = 0;
	for (i = 0; i < numc; i++)
	{
		if (components[i].size() >= minnumelements && components[i].size() <= maxnumelements)
		{
			index++;
		}
	}

	int cdims[] = {index};
	plhs[0] = mxCreateCellArray(1,cdims);
	
	index = 0;
	for (i = 0; i < numc; i++)
	{
		mxArray *ll = 0;
		if (components[i].size() >= minnumelements && components[i].size() <= maxnumelements)
		{
			ll = mxCreateNumericMatrix(2,components[i].size(),mxGetClassID(Points),mxREAL);
			REAL *dat = (REAL*) mxGetData(ll);
			for (int k = 0; k < components[i].size(); k++)
			{
				dat[2*k] = components[i][k];
				dat[2*k+1] = ccana.particles[components[i][k]].numerator;
			}
			qsort(dat,components[i].size(),sizeof(REAL)*2,cmpfloat2);
			mxSetCell(plhs[0],index++,ll);
		}
	}

	fprintf(stderr,".\n"); fflush(stderr);

}


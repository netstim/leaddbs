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
#include <string.h>
#include <sys/time.h>

using namespace std;

#define REAL float
#define PI 3.1415926536
//#define INFINITY 99999999999.0




#define ntype 100000
class PropStats
{
    int movetype[ntype];
    int movepos;
    
    
	public:
	void clear() { N = 0; accept = 0; for(int k = 0; k< ntype;k++) movetype[k] = 0; movepos = 0;}
	
    void propose() {N++;}
    void accepted() { accept++; }

    void accepted_up()   {  movetype[movepos%ntype] = 1; movepos++; accept++; N++;}
    void accepted_down() {  movetype[movepos%ntype] = -1;  movepos++; accept++; N++;}
    void rejected()      {  N++;}
	void report(const char *s)
	{
		mexPrintf("%s  #proposals: %8.2fk  acceptratio: %.2f   \n",s,1.0*N/1000.0,100.0*accept/N);
	}
    
    double updownratio;
    int up;
    int down;
	int N;
	int accept;
    
    void compute_updownratio()
    {
        up = 0;
        down  = 0;
        for (int k = 0; k < ntype; k++)
        {
            if (movetype[k] == 1)
                up++;
            if (movetype[k] == -1)
                down++;
        }
        updownratio = double(up)/(1+down);
    }
    
	void report_movetype(const char *s)     
	{
        compute_updownratio();
		fprintf(stderr,"%s  #proposals: %8.2fk  acceptratio: %.2f  up %i down %i  updown-ratio: %.3f   \n",s,1.0*N/1000.0,100.0*accept/N,up,down,updownratio );

	}
};





//#define TIMING



#ifdef TIMING

static struct timeval timeS;



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

Timing externalenergy_time;
Timing internalenergy_time;
Timing odfeval_time;
Timing total_time;

Timing shiftproposal_time;
Timing birthproposal_time;
Timing deathproposal_time;
Timing capproposal_time;
Timing lenproposal_time;
Timing connproposal_time;

PropStats deathstats;
PropStats birthstats;
PropStats connstats;
PropStats shiftstats;
PropStats capstats;
PropStats lenstats;


#endif

PropStats stats;


#include "MersenneTwister.h"
MTRand mtrand;
REAL *BESSEL_APPROXCOEFF;


#include "EnergyComputer_connec.cpp"
//#include "EnergyComputer_center.cpp"
//#include "RJMCMC_singlegradprop.cpp"
#include "RJMCMC_randshift.cpp"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 7 && nrhs != 8) {
	
	printf(" wrong usage!!!.\n\n");
    return;
	} else if(nlhs>3) {
	printf("Too many output arguments\n");
    return;
	}
	
	int pcnt = 0;
	const mxArray *Points;
	Points = prhs[pcnt++];       
	int numPoints = mxGetN(Points);
	REAL *points = (REAL*) mxGetData(Points);
	
	const mxArray *DataImg;
	DataImg = prhs[pcnt++];       
	REAL *dimg = (REAL*) mxGetData(DataImg);
	const int *dsize = mxGetDimensions(DataImg);
    
	const mxArray *Mask;
	Mask = prhs[pcnt++];       
	REAL *mask = (REAL*) mxGetData(Mask);
    const int *dmsize = mxGetDimensions(Mask);
    int mask_oversamp_mult = dmsize[0]/dsize[1];
    
    
    
	const mxArray *VoxSize;
	VoxSize = prhs[pcnt++];       
	double *voxsize = (double*) mxGetPr(VoxSize);
            
    

	const mxArray *Params;
	Params = prhs[pcnt++];       
	double *params = (double*) mxGetPr(Params);

	REAL Temp = (REAL) params[0];	
	REAL numit = (REAL) params[1];
	REAL conprob = (REAL) params[2];
	REAL particle_weight = (REAL) params[3];
	REAL particle_width = (REAL) params[4];
	REAL particle_len = (REAL) params[5];
	REAL chempot_connection = (REAL) params[6];
	REAL chempot_particle = (REAL) params[7];
	REAL inex_balance = (REAL) params[8];
	REAL chempot2 = (REAL) params[9];
	REAL meanval_sq = (REAL) params[10];

    const mxArray *BesselExpansion;
	BesselExpansion = prhs[pcnt++];       
	BESSEL_APPROXCOEFF = (REAL*) mxGetData(BesselExpansion);
    
    
    
	// read spherical-interpolator data

	const mxArray *sinterpstruct = prhs[pcnt++];
	mxArray *Indices = mxGetField(sinterpstruct,0,"indices");
	mxArray *BaryCoords = mxGetField(sinterpstruct,0,"barycoords");
	mxArray *Beta = mxGetField(sinterpstruct,0,"beta");
	mxArray *NumInterpPoints = mxGetField(sinterpstruct,0,"numpoints");
	
	REAL *indimg = (REAL*) mxGetData(Indices);
	const int *isize = mxGetDimensions(Indices);
	int totsz = isize[0]*isize[1]*isize[2]*isize[3];
	int *indeximg = (int*) malloc(sizeof(int)*totsz);
	for (int k =0;k < totsz;k++)
		indeximg[k] = int(indimg[k])-1;      
	REAL *barycoords = (REAL*) mxGetData(BaryCoords);
	REAL *beta = (REAL*) mxGetData(Beta);
	int nip = int(*((REAL*)mxGetData(NumInterpPoints)));
	
    SphereInterpolator *sinterp = new SphereInterpolator(barycoords,indeximg,nip,isize[2],beta[0]);
	
	double breakhandle = 0;	
	const mxArray *BreakHandle;	
	if (nrhs == 8)
	{
		BreakHandle = prhs[pcnt++];
		breakhandle = *mxGetPr(BreakHandle);
	}

	#ifdef TIMING
	externalenergy_time.clear();
	internalenergy_time.clear();
	odfeval_time.clear();
	total_time.clear();

	shiftproposal_time.clear();
	birthproposal_time.clear();
	deathproposal_time.clear();
	capproposal_time.clear();
	lenproposal_time.clear();
	connproposal_time.clear();

	deathstats.clear();
	birthstats.clear();
	connstats.clear();
	shiftstats.clear();
	capstats.clear();
	lenstats.clear();
	#endif
	
    stats.clear();


	REAL cellsize_pl = 2*particle_len;
	REAL cellsize_wi = 6*particle_width;
	
	REAL cellsize = (cellsize_pl > cellsize_wi)?cellsize_pl:cellsize_wi;
		
	
	REAL curv_hardthres = 0.7; // deactivated

    fprintf(stderr,"setting up MH-sampler \n"); fflush(stderr);
	RJMCMC sampler(points,numPoints, dimg, dsize, voxsize, cellsize);
    fprintf(stderr,"setting up Energy-computer \n"); fflush(stderr);
    EnergyComputer encomp(dimg,dsize,voxsize,sinterp,&(sampler.pcontainer),mask,mask_oversamp_mult);

    fprintf(stderr,"setting up parameters\n"); fflush(stderr);
	sampler.setParameters(Temp,numit,conprob,particle_len,curv_hardthres,chempot_particle);
   	sampler.setEnergyComputer(&encomp);
    encomp.setParameters(particle_weight,particle_width,chempot_connection*particle_len*particle_len,particle_len,curv_hardthres,inex_balance,chempot2,meanval_sq);
    
    fprintf(stderr,"starting to iterate\n"); fflush(stderr);
	sampler.iterate(breakhandle);
	
	int cnt = sampler.pcontainer.pcnt;


	
	int dims[] = {sampler.attrcnt, sampler.pcontainer.pcnt};
	plhs[0] = mxCreateNumericArray(2,dims,mxGetClassID(Points),mxREAL);
	REAL *npoints = (REAL*) mxGetData(plhs[0]);	
	sampler.writeout(npoints);
	
    
    
    
    
    stats.compute_updownratio();
    
    const char *field_names[] = {"up","down","accepted","rejected","iterations"};
    int dim = 1;
	mxArray* mxStat=mxCreateStructArray(1,&dim,5,field_names);
    plhs[1] = mxStat;
    
    int outcnt = 0;
    mxArray *out[20];
    double *s;
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = stats.up;
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = stats.down;
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = stats.accept;
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = stats.N - stats.accept  ;
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = sampler.iterations_done;

    for (int k = 0; k < outcnt; k++)
        mxSetField( mxStat,0,field_names[k], out[k]);
	

    
    
    
    
    
    
    
    
    
    
    
    

	delete sinterp;
	free(indeximg);


}


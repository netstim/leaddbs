
#include <mex.h>
#include <matrix.h>
#include <math.h>
#include <map>
#include <string.h>

using namespace std;

#define REAL float
#define PI 3.1415926536

////////////////////////////////////

#include "MersenneTwister.h"
MTRand mtrand;
int dbgflag;


#ifndef LINUX_MACHINE    
    #define INFINITY 9999999999999  // not defined on a windows pc
#else
    #include <sys/time.h>
#endif


////////// to monitor statistics
#define ntype 100000
class PropStats
{
    int movetype[ntype];
    int movepos;
    
	public:
   PropStats()
   {
        lam_deltaE = 0.99999;
        deltaE = 0;
   }
        
        
    double deltaE;
    double lam_deltaE;
	void clear() { N = 0; accept = 0; for(int k = 0; k< ntype;k++) movetype[k] = 0; movepos = 0; 
              }
	
    void propose() {N++;}
    void accepted(double dE) { 
        if (dE > 0)
            movetype[movepos%ntype] = -1;
        else
            movetype[movepos%ntype] = 1;
        if (dE != INFINITY)
            deltaE = dE*(1-lam_deltaE)+deltaE*lam_deltaE;
        accept++; N++; movepos++;
    }

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
    double indi;
    
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
        indi = double(down);
    }
    
	void report_movetype(const char *s)     
	{
        compute_updownratio();
		fprintf(stderr,"%s #prps: %8.2fk arat: %4.2f dE: %f \n ",s,1.0*N/1000.0,100.0*accept/N,deltaE*100000);

	}
};


PropStats deathstats;
PropStats birthstats;
PropStats connstats;
PropStats shiftstats;
PropStats capstats;
PropStats vfstats;
PropStats stats;



#include "EnergyComputer_connec.cpp"
#include "RJMCMC_randshift.cpp"




void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 8 && nrhs != 9) {
	
	printf(" wrong usage!!!.\n\n");
    return;
	} else if(nlhs>5) {
	printf("Too many output arguments\n");
    return;
	}

    
	int pcnt = 0;
    
    //////////////////////////// get all the arguments from MATLAB

    // the segments
    const mxArray *Points;
	Points = prhs[pcnt++];       
	int numPoints = mxGetN(Points);
	REAL *points = (REAL*) mxGetData(Points);
    
    // the volume fraction maps
	const mxArray *VFmap;
	VFmap = prhs[pcnt++];       
	REAL *vfmap_in = (REAL*) mxGetData(VFmap);
    const int *vfmapsize = mxGetDimensions(VFmap);
    
    // the data
	const mxArray *DataImg;
	DataImg = prhs[pcnt++];       
	REAL *dimg = (REAL*) mxGetData(DataImg);
	const int *dsqsize = mxGetDimensions(DataImg);
    
    // the WM-mask
	const mxArray *Mask;
	Mask = prhs[pcnt++];       
	REAL *mask = (REAL*) mxGetData(Mask);
    const int *dmsize = mxGetDimensions(Mask);
    
    // index array where the data is located and a mean-signal map
	const mxArray *S2andIDX;
	S2andIDX = prhs[pcnt++];       
	REAL *s2andidx = (REAL*) mxGetData(S2andIDX);
    const int *datasize = mxGetDimensions(S2andIDX);
    int mask_oversamp_mult = dmsize[0]/datasize[0];
    
    
    // the  voxel size (our coordinates are in mm)
	const mxArray *VoxSize;
	VoxSize = prhs[pcnt++];       
	double *voxsize = (double*) mxGetPr(VoxSize);
            
    // read parameters (ParamStruct -> auxiliary_classes)
	const mxArray *Params;
	Params = prhs[pcnt++];           
    ParamStruct info;
    info.readParameters(Params);
 
    
	// read spherical-interpolator data
	const mxArray *sinterpstruct = prhs[pcnt++];
    SphereInterpolator *sinterp = new SphereInterpolator(sinterpstruct);

 
   // a size vector of the full data [N_approx M_sphere w h d] (actually the data itself does not exist anymore in that way, because data it's kept sparsely)
 
   int lmax = int(info.lmax);
    int numbshells = int(info.numbvals);
    const int datacombisize[5] = {(lmax/2+1)*numbshells ,sinterp->nverts, datasize[0],datasize[1],datasize[2]};
     
//    const int datacombisize[5] = {dsqsize[0],dsqsize[1], datasize[0],datasize[1],datasize[2]};
     
    
    // if this matlab-handle disappears tracking is stopped
	double breakhandle = 0;	
	const mxArray *BreakHandle;	
	if (nrhs == 9)
	{
		BreakHandle = prhs[pcnt++];
		breakhandle = *mxGetPr(BreakHandle);
	}

	////////////////////////////////////
            
    
    // either initlize vfmap or just get it from previous step
    REAL *vfmap;
    {
    int vfmapsize[5] =  { int((REAL)datasize[0]/info.particle_width),
                          int((REAL)datasize[1]/info.particle_width),
                          int((REAL)datasize[2]/info.particle_width), 
                          int(lmax/2+1) , int(info.numbvals+1)};
    plhs[1] = mxCreateNumericArray(5,vfmapsize,mxGetClassID(Points),mxREAL);
    vfmap = (REAL*) mxGetData(plhs[1]);
    }       
       
	// initialize everything
	REAL cellsize = info.particle_len*2;
    REAL cellsize2[3] = {info.particle_width*voxsize[0] , info.particle_width*voxsize[1], info.particle_width*voxsize[2]};
	
    fprintf(stderr,"setting up MH-sampler \n"); fflush(stderr);
	RJMCMC sampler(points,numPoints, vfmap, dimg, datacombisize, voxsize, cellsize,cellsize2,info.particle_len,info.numcores);

    fprintf(stderr,"setting up Energy-computer \n"); fflush(stderr);
    EnergyComputer encomp(dimg,datacombisize,voxsize,sinterp,&(sampler.pcontainer),vfmap,mask,s2andidx,mask_oversamp_mult);

    fprintf(stderr,"setting up parameters\n"); fflush(stderr);
   	sampler.setEnergyComputer(&encomp);
    encomp.setParameters(info);	
    sampler.setParameters(info);
    
    
    //----------- start tracking
    fprintf(stderr,"starting to iterate\n"); fflush(stderr);
	sampler.iterate(breakhandle);
    //-------------------------
    
    
    // compute some more maps needed in matlab to calculate all model parameters
    
    
    REAL *vf_fibs = 0;
    if (nlhs >=5)
    {
        int vf_fibs_size[3] = {  sampler.pcontainer.pcnt,
                              int(lmax/2+1) , int(info.numbvals+1)};
        plhs[4] = mxCreateNumericArray(3,vf_fibs_size,mxGetClassID(Points),mxREAL);
        vf_fibs = (REAL*) mxGetData(plhs[4]);
    }
    
    
    encomp.computeMeanProjMap(vf_fibs);

    
    // label fibers uniqely.    
    sampler.label_fibers();
    if (nlhs >=4)
    {
        vector<RJMCMCBase::Cross> c = sampler.createCrossingList();
        plhs[3] = mxCreateNumericMatrix(4,c.size(),mxGetClassID(Points),mxREAL);
        REAL *arr = (REAL*) mxGetData(plhs[3]);	
        int cn = c.size();
        for (int k = 0; k< cn;k++)
        {
            arr[4*k] = c[k].i;
            arr[4*k+1] = c[k].j;
            arr[4*k+2] = c[k].wi;
            arr[4*k+3] = c[k].wj;            
        }
        
    }
    
    
	// write tracking state to matlab variable
	int cnt = sampler.pcontainer.pcnt;	
	int dims[] = {sampler.attrcnt, sampler.pcontainer.pcnt};
	plhs[0] = mxCreateNumericArray(2,dims,mxGetClassID(Points),mxREAL);
	REAL *npoints = (REAL*) mxGetData(plhs[0]);	
	sampler.writeout(npoints);
		        
    // gather statistitcs and put it into matlab structure
    stats.compute_updownratio();
    birthstats.compute_updownratio();
    deathstats.compute_updownratio();    
    const char *field_names[] = {"up","down","accepted","rejected","iterations","deltaE"};
    int dim = 1;
	mxArray* mxStat=mxCreateStructArray(1,&dim,6,field_names);
    plhs[2] = mxStat;    
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
    out[outcnt] = mxCreateNumericArray(1,&dim,mxDOUBLE_CLASS,mxREAL); s =  (double*) mxGetData(out[outcnt]); outcnt++;
    *s = stats.deltaE;
    for (int k = 0; k < outcnt; k++)
        mxSetField( mxStat,0,field_names[k], out[k]);
	
    
    // finally free some non-MATLAB mallocs
	delete sinterp;
	//free(indeximg);

}


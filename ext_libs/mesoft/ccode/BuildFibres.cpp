
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
		attrcnt = 15;
		for (int k = 0; k < numPoints; k++)
			{
				ParticleL *p = &(particles[k]);
				p->R = pVector(points[attrcnt*k], points[attrcnt*k+1],points[attrcnt*k+2]);
				p->N = pVector(points[attrcnt*k+3],points[attrcnt*k+4],points[attrcnt*k+5]);				
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
	printf("\nUsage: Df = STderivative(f)\n\n",nrhs);
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

	mwSize cdims[] = {index};
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


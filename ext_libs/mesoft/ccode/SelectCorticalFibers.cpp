
#include "mex.h"
#include "matrix.h"
#include <math.h>
#include <map>
#include <vector>
#include <string.h>
// #include <sys/time.h>

using namespace std;

#define REAL float
#define PI 3.1415926536








void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 3) {

    return;
	} else if(nlhs>2) {
	printf("Too many output Arguments %i\n",nlhs);
    return;
	}
    
	
    
    
  
    
	int pcnt = 0;
	const mxArray *Labeling;
	Labeling = prhs[pcnt++];       
	REAL *labeling = (REAL*) mxGetData(Labeling);
    const mwSize *dims = mxGetDimensions(Labeling);
    

    REAL Nsize = REAL(*mxGetPr(prhs[pcnt++]));
      
  	const mxArray *Conns;
	Conns = prhs[pcnt++];       
	int numConns = mxGetN(Conns);
	REAL *conns = (REAL*) mxGetData(Conns);  
    
    
     
    
    
    mexPrintf("total number fibers: %i \n",numConns/2);
    
    const mwSize dims2[] = {static_cast<mwSize>(numConns/2),1};
    plhs[0] = mxCreateNumericArray(2,dims2,mxGetClassID(Labeling),mxREAL);
	REAL *inc = (REAL*) mxGetData(plhs[0]);	
    
    double Nsize2 = Nsize*Nsize;
   
    int totsz = dims[0]*dims[1]*dims[2];
   
    for (int k = 0; k < numConns/2;k++)
    {
            
        if (k%(numConns/15) == 0)
        {
            mexPrintf("%.0f%%, ",100.0*(float)k/numConns*2);
            mexEvalString("drawnow;");
        }
                
        const REAL d = 0.5;
        int idx = int(d+conns[k*6])+int(d+conns[k*6+1])*dims[0]+int(d+conns[k*6+2])*dims[0]*dims[1];
        if (idx < 0 || idx >= totsz)
            continue;
        if (labeling[idx] == 0)
            continue;
        idx = int(d+conns[k*6+3])+int(d+conns[k*6+4])*dims[0]+int(d+conns[k*6+5])*dims[0]*dims[1];
        if (idx < 0 || idx >= totsz)
            continue;
        if (labeling[idx] == 0)
            continue;

         inc[k] = 1;
            
        
    }
    
    
}



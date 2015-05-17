
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 1) {	
	printf(" wrong usage!!!.\n\n");
    return; }
	
	int pcnt = 0;
	const mxArray *Str;
	Str = prhs[0];       

    char *str;

    str = (char *) mxCalloc(mxGetN(Str)+1, sizeof(char));
    mxGetString(Str, str, mxGetN(Str)+1);


      

    fprintf(stderr,"%s\n",str);
    
    
}
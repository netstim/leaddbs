/* MEX file to take a set of grey points and a set of underlying mesh points and */
/* compute which grey point maps to which mesh point - no interpolation! */
/* v0.1 A.Wade 11/2000 */


#include <math.h>
#include "mex.h"

void getNearest(double *inpPtr,double *refPtr,double *outPtr, double *distPtr, unsigned int inpRows,unsigned int refRows)
{

	// Run through each member of inpPtr and take the distance between it and
	// each member of refPtr. Return the indices of the nearest points in refPtr

	unsigned int RefCounter,InpCounter;
	double cLowest,thisDist,distX,distY,distZ, xInp,yInp,zInp; 
	unsigned int twiceRefRows,twiceInpRows,indLowest;
	
	twiceRefRows=2*refRows;
	twiceInpRows=2*inpRows;

    for (InpCounter=0;InpCounter<inpRows;InpCounter++) {
		
		
		xInp=inpPtr[InpCounter]; // Matlab orders its data so that the row dimension increases fastest
		yInp=inpPtr[InpCounter+inpRows];
		zInp=inpPtr[InpCounter+twiceInpRows];

		cLowest=9999999;
		indLowest=0;

		for (RefCounter=0;RefCounter<refRows;RefCounter++) { // Should reverse loop for speed improvement
			distX=refPtr[RefCounter]-xInp;
			distY=refPtr[RefCounter+refRows]-yInp;
			distZ=refPtr[RefCounter+twiceRefRows]-zInp;
			thisDist=(distX*distX+distY*distY+distZ*distZ);

			if (thisDist<cLowest) {
				cLowest=thisDist;
				indLowest=RefCounter;
			} // end if

		} // next RefCounter

		*(outPtr++)=indLowest+1; // add the one because matlab references arrays from 1
		*(distPtr++)=cLowest; // Returns the squared distances
	} // next InpCounter


} // end of function

void mexFunction( int nlhs, 
		  mxArray *plhs[], 
		  int nrhs, 
		  const mxArray *prhs[] )
{ 
    double *yp; 
    double *t,*y, *refPtr, *inpPtr, *outPtr, *distPtr; 
    unsigned int refRows,refCols,inpRows,inpCols; 

    /* Check for proper number of arguments */

    if (nrhs != 2) { 
		mexErrMsgTxt("Two input vectors required."); 
    } else if (nlhs > 2) {
		mexErrMsgTxt("Error! Too many output arguments."); 
    } 

    
  /*  Get the dimensions and a pointer to the reference matrix. */
	refRows = mxGetM(prhs[0]);
	refCols = mxGetN(prhs[0]);
	refPtr = mxGetPr(prhs[0]);

  /*  Get the dimensions and a pointer to the input matrix. */
	inpRows = mxGetM(prhs[1]);
	inpCols = mxGetN(prhs[1]);
	inpPtr = mxGetPr(prhs[1]);

  
 	// Check for type and size

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) { 
		mexErrMsgTxt("Argument 1 must be a real double"); 
    } 
    
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) { 
		mexErrMsgTxt("Argument 2 must be a real double"); 
    } 

	if ((refCols!=3) || (inpCols!=3)) {
		mexErrMsgTxt("Both input arguments must have 3 columns"); 
	}

    //printf("Got here ... assigning");
    /* Create a matrix for the return argument */ 
    plhs[0] = mxCreateDoubleMatrix(inpRows, 1, mxREAL); 
    plhs[1] = mxCreateDoubleMatrix(inpRows, 1, mxREAL); 

  
    outPtr=mxGetPr(plhs[0]);
    distPtr=mxGetPr(plhs[1]);
    /* Assign pointers to the various parameters */ 

    if ((!outPtr) || (!distPtr)) { 
		mexErrMsgTxt("Could not assign output matrices. Out of memory?"); 
    } 

  //  printf("Got here 232");
    /* Do the actual computations in a subroutine */
    getNearest(inpPtr,refPtr,outPtr,distPtr,inpRows,refRows); 
//  printf("Got here dfjdkj");
  
    
	// Don't compute the flops
  return;
    
}






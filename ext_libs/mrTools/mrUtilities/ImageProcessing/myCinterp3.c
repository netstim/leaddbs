/*  Matlab MEX file.

05/04/99 - ON - Corrected bug in interpolation (switched wa and wb)
                and some changes to improve speed (register variables,
                and avoiding calling function modf in the loop).
*/

/* #define Storage */
/* #define CAP */
/* #include "machine.h" */
/* #include "debug.h" */

#include "mex.h"
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>

void myCinterp3(double *vol,int nslice,int size1,int size2,int npts,
            double *x, double *f,double badval)
	    
{
   register int  i, saglen,index;
   register double wa,wb,wc;
   register double a,b,c;
   register double xM, yM, zM;
   register double *volR = vol;
   register double *xR   = x;
   register double *fR   = f;
   register int    nptsR = npts;
   register double badvalR = badval;
   register int size1R = size1;
   register int size2R = size2;

   saglen=size1*size2;

   for (i=0;i<nptsR;i++) {
     xM = *(xR+i);
     yM = *(xR+nptsR+i);
     zM = *(xR+2*nptsR+i);

     // handle boundary conditions properly. if you are
     // right exactly on the edge, then move the coordinate
     // over by just a little, so that you don't try to
     // dereference a value outside the volume (even though
     // you would be multiplying that value by a weight of 0).
     if (xM == size2R) xM=xM-0.00000001;
     if (yM == size1R) yM=yM-0.00000001;
     if (zM == nslice) zM=zM-0.00000001;

     /* distances to nearest lower integer coordinate */
     a = (int) xM; wa = xM - a; a--;
     b = (int) yM; wb = yM - b; b--;
     c = (int) zM; wc = zM - c; c--;

     /*     wa=modf(*(xR+i),&a);a--;
     wb=modf(*(xR+nptsR+i),&b);b--;
     wc=modf(*(xR+2*nptsR+i),&c);c--;*/
   if (xM<1 || xM>size2R || 
       yM<1 || yM>size1R ||
       zM<1 || zM>nslice ) fR[i]=badvalR;
     else { 
      
     index=c*saglen+a*size1R+b;

     fR[i] = (1-wc) * ( (1-wa) * ( (1-wb) * (*(volR+index)) 
                                 +   wb  * (*(volR+index+1)) )
                      +   wa  * ( (1-wb) * (*(volR+index+size1R))
                                 +   wb  * (*(volR+index+size1R+1))  ))
           +   wc  * ( (1-wa) * ( (1-wb) * (*(volR+index+saglen)) 
                                 +  wb  * (*(volR+index+saglen+1)) ) 
                      +   wa  * ( (1-wb) * (*(volR+index+saglen+size1R))
                                 +  wb  * (*(volR+index+saglen+size1R+1)) ) ); 


     }
 }

}

/* #define output plhs[0]; */

void mexFunction(int nlhs,   /* number of arguments on lhs */
		 mxArray	*plhs[],   /* Matrices on lhs      */
		 int nrhs,	   /* no. of mat on rhs    */
		 const mxArray	*prhs[]    /* Matrices on rhs      */
		 )
{
  int npts = 0;
  double *vol,*nslice,*size1,*size2,*x,*f,badval;
  double *sagSize;

  /* Check for proper number of arguments */

  if (nrhs==0) { /* help */
   printf("myCinterp3(volume,sagSize,numSlices,samp,badval)\n");
   printf("\n  volume: volume anatomy data with size 1xprod(sagSize)*numSlices\n");
   printf("  sagSize: size of sagittal anatomy slices\n");
   printf("  numSlices: number of sagittal anatomy slices\n");
   printf("  samp:  data points to interpolate with size prod(sagSize)x3\n");
   printf("  badval (optional): returned for points outside anatomy cube.  0.0 is default\n");
   
   
 }
 
  else {
    if (nrhs < 4) {
      mexErrMsgTxt("CInterp3.c needs at least four arguments.");
  }
  if( (x = mxGetPr(prhs[0])) == NULL) {
    mexErrMsgTxt("myCinterp3: null ptr for input matrix.");
  }
    
  /* Call our routine to read the cap stuff */
  vol=mxGetPr(prhs[0]);
  sagSize=mxGetPr(prhs[1]);
  nslice=mxGetPr(prhs[2]);
  x=mxGetPr(prhs[3]);
  npts=mxGetM(prhs[3]);
  if (nrhs==5) {
    badval=*mxGetPr(prhs[4]);
  }
  else {
    badval=0.0;
  }

  plhs[0] = mxCreateDoubleMatrix(1,npts,mxREAL);

  f = mxGetPr(plhs[0]);

  myCinterp3(vol,(int) *nslice,(int) *sagSize,(int) *(sagSize+1),npts,x,f,badval);
  }      

}


#include "mex.h"
#include "matrix.h"
#include <math.h>

#define REAL double


  double BESSI0(double X);
  double BESSI1(double X);

// ---------------------------------------------------------------------
  double BESSI(int N, double X) {
/*----------------------------------------------------------------------
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows. 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
------------------------------------------------------------------------*/

      int IACC = 40; 
	  double BIGNO = 1e10, BIGNI = 1e-10;
      double TOX, BIM, BI, BIP, BSI;
      int J, M;

      if (N==0)  return (BESSI0(X));
      if (N==1)  return (BESSI1(X));
      if (X==0.0) return 0.0;

      TOX = 2.0/X;
      BIP = 0.0;
      BI  = 1.0;
      BSI = 0.0;
      M = (int) (2*((N+floor(sqrt(IACC*N)))));
      for (J = M; J>0; J--) {
        BIM = BIP+J*TOX*BI;
        BIP = BI;
        BI  = BIM;
        if (fabs(BI) > BIGNO) {
          BI  = BI*BIGNI;
          BIP = BIP*BIGNI;
          BSI = BSI*BIGNI;
        }
        if (J==N)  BSI = BIP;
      }
      return (BSI*BESSI0(X)/BI);
  }

// ----------------------------------------------------------------------
//  Auxiliary Bessel functions for N=0, N=1
  double BESSI0(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=1.0; P2=3.5156229; P3=3.0899424; P4=1.2067492;
      P5=0.2659732; P6=0.360768e-1; P7=0.45813e-2;
      Q1=0.39894228; Q2=0.1328592e-1; Q3=0.225319e-2;
      Q4=-0.157565e-2; Q5=0.916281e-2; Q6=-0.2057706e-1;
      Q7=0.2635537e-1; Q8=-0.1647633e-1; Q9=0.392377e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return (P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

// ---------------------------------------------------------------------
  double BESSI1(double X) {
      double Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX;
      P1=0.5; P2=0.87890594; P3=0.51498869; P4=0.15084934;
      P5=0.2658733e-1; P6=0.301532e-2; P7=0.32411e-3;
      Q1=0.39894228; Q2=-0.3988024e-1; Q3=-0.362018e-2;
      Q4=0.163801e-2; Q5=-0.1031555e-1; Q6=0.2282967e-1;
      Q7=-0.2895312e-1; Q8=0.1787654e-1; Q9=-0.420059e-2;
      if (fabs(X) < 3.75) {
        Y=(X/3.75)*(X/3.75);
        return(X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))));
      }
      else {
        AX=fabs(X);
        Y=3.75/AX;
        BX=exp(AX)/sqrt(AX);
        AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))));
        return (AX*BX);
      }
  }

REAL ricfit(REAL *a, REAL *w,int kvol)
{
    REAL X = 0;
    for (int k = 0;k < kvol;k++)
        X += w[k]*a[k];
    
    for (int it = 0; it < 100; it++)
    {
        REAL fs = 0;
        REAL fss = 0;
        for (int k = 0; k < kvol; k++)
        {
            REAL R;
            if (X*a[k]>100)
                R = 1;
            else
                R = BESSI1(X*a[k])/BESSI0(X*a[k]);
            fs += w[k]*(X-a[k]*R);
            fss += w[k]*( (a[k]*R)*(a[k]*R) - (a[k]*a[k] - a[k]*R/X) + 1);
        }
        REAL delta = fs/fss;
        X = X - delta;
        if (delta < 0.001)
            break;
    }
    return X;
}











void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
 	
    if(nrhs != 2) {
	
	printf(" wrong usage!!!.\n\n");
    return;
	} else if(nlhs>1) {
	printf("Too many output arguments\n");
    return;
	}
	
	int pcnt = 0;
	const mxArray *Data;
	Data = prhs[pcnt++];       
	REAL *data = (REAL*) mxGetData(Data);
    const int *dims_data = mxGetDimensions(Data);

	const mxArray *Ker;
	Ker = prhs[pcnt++];       
	REAL *ker = (REAL*) mxGetData(Ker);
    const int *dims_ker = mxGetDimensions(Ker);
    
    int w = dims_data[0];
    int h = dims_data[1];
    int d = dims_data[2];
    
    
    int dims[] = {w, h, d};
	plhs[0] = mxCreateNumericArray(3,dims,mxGetClassID(Data),mxREAL);
	REAL *out = (REAL*) mxGetData(plhs[0]);	
    
    int ksz = (dims_ker[0]-1)/2;
    
    int kvol = dims_ker[0]*dims_ker[0]*dims_ker[0];
    
    REAL dtmp[kvol];
    
    for (int z = ksz; z < d-ksz; z++)
    {
        mexPrintf(".");
        mexEvalString("drawnow");
        for (int y = ksz; y < h-ksz; y++)
        {
            
            for (int x = ksz; x < w-ksz; x++)
                {
                    int cnt = 0;
                    for (int az = -ksz; az <= ksz; az++)
                        for (int ay = -ksz; ay <= ksz; ay++)
                            for (int ax = -ksz; ax <= ksz; ax++)
                                {
                                    dtmp[cnt] = data[x+ax + (y+ay)*w + (z+az)*w*h];
                                    cnt++;
                                }
                    out[x + y*w + z*w*h] = ricfit(dtmp,ker,kvol);
                }
        }       
    }
}    


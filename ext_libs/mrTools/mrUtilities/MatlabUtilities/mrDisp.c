#ifdef documentation
=========================================================================

     program: mydisp.c
          by: justin gardner
     purpose: print w/out newline character
        date: 07/08/03
     compile: mex mydisp.c
     
=========================================================================
#endif

////////////////////
// include section
////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "mex.h"

///////////////////
// define section
///////////////////
#define STRSIZE 2048

///////////////////
// function decls
///////////////////
void usageError();

//////////////////////////////////////////
// function mexFunction called by matlab
//////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  char str[STRSIZE];

  // check input arguments
  if (nrhs != 1){
    usageError();
    return;
  }

  // get string
  mxGetString(prhs[0], str, mxGetN(prhs[0])+1);

  // print string
  printf("%s",str);
  fflush(stdout);
}

////////////////////////
// function usageError
////////////////////////
void usageError()
{
  printf("USAGE: mydisp('string')\n");
}


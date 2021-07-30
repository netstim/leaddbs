// isMatlabVer.c
// Compare Matlab version to specified number
// Match = isMatlabVer(Relop, N)
// INPUT:
//   Relop: Comparison operator as string: '<', '<=', '>', '>=', '=='.
//   N:     Number to compare with as DOUBLE vector with 1 to 4 elements.
//
// OUTPUT:
//   Match: Locical scalar, TRUE for matching comparison, FALSE otherwise.
//
// EXAMPLES:
//   version ==> '7.8.0.342 (R2009a)'  (different results for other version!)
//   isMatlabVer('<=', 7)               % ==> TRUE
//   isMatlabVer('>',  6)               % ==> TRUE
//   isMatlabVer('<',  [7, 8])          % ==> FALSE
//   isMatlabVer('<=', [7, 8])          % ==> TRUE
//   isMatlabVer('>',  [7, 8, 0, 342])  % ==> FALSE
//   isMatlabVer('==', 7)               % ==> TRUE
//   isMatlabVer('==', [7, 10])         % ==> FALSE
//   isMatlabVer('>',  [7, 8, 0])       % ==> FALSE (the 342 is not considered!)
//
// NOTES: The C-Mex function takes about 0.6% of the processing time needed
//   by Matlab's VERLESSTHAN, which can check other toolboxes also.
//
// Compile:   mex isMatlabVer.c
//
// Compiler: GCC 9.3.0, Clang 12.0, MSVC 2019 (M/T)
// Author: Jan Simon, Heidelberg, (C) 2010 matlab.THISYEAR(a)nMINUSsimon.de
// License: BSD - use, copy, modify on own risk, mention the author.
//
// See also: VER, VERLESSTHAN.

/*
% $JRev: R0d V:003 Sum:PESD4Qe2aN2S Date:12-Apr-2010 16:13:38 $
% $License: BSD $
% $File: Tools\Mex\Source\isMatlabVer.c $
% History:
% 001: 12-Apr-2010 03:48, Replace weak SSCANF(VERSION, '%f', 1) < X.
*/

/*
% Rev 14-April-2021 (Ningfei Li, ningfei.li(a)gmail.com):
%   Clean up.
%   Use "u" instead of "L" for string literal (char16_t). 
%   Update compiled binaries.
*/

#include "mex.h"
#include <math.h>
#include <string.h>

static int MatlabV[4] = {0, 0, 0, 0};

void Init(void);
void BadRelopError(void);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxChar *Relop;
  mwSize nData, i;
  int    c;
  double *Data;
  
  // Proper number of arguments:
  if (nrhs != 2) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadNInput",
                       "2 inputs required.");
  }
  if (nlhs > 1) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadNOutput",
                       "1 output allowed.");
  }
  
  // Type of inputs:
  if (!mxIsChar(prhs[0])) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadInput1",
                       "1st input must be a string.");
  }
  if (!mxIsDouble(prhs[1])) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadInput2",
                       "2nd input must be DOUBLE.");
  }
  nData = mxGetNumberOfElements(prhs[1]);
  if (nData == 0 || nData > 4) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadInput2Len",
                       "2nd input must have 1 to 4 elements.");
  }
  
  // Initialize - get Matlab version:
  if (MatlabV[0] == 0) {
     Init();
  }
  
  // Get input data:
  Data = mxGetPr(prhs[1]);
  Relop  = (mxChar *) mxGetData(prhs[0]);
       
  // Check first element to be integer:
  if (*Data != floor(*Data)) {
     mexErrMsgIdAndTxt("JSim:isMatlabVer:BadInput2Value",
                       "Input N must contain integer values.");
  }

  // Start the comparison:
  switch (mxGetNumberOfElements(prhs[0])) {
     case 1:
        if (!memcmp(Relop, u"<", sizeof(mxChar))) {         // MatlabV < N
           for (i = 0; i < nData; i++) {
              if ((c = MatlabV[i] - (int) Data[i]) != 0) {
                 plhs[0] = mxCreateLogicalScalar(c < 0);
                 return;
              }
           }
           plhs[0] = mxCreateLogicalScalar(c < 0);
    
        } else if (!memcmp(Relop, u">", sizeof(mxChar))) {  // MatlabV > N
           for (i = 0; i < nData; i++) {
              if ((c = MatlabV[i] - (int) Data[i]) != 0) {
                 plhs[0] = mxCreateLogicalScalar(c > 0);
                 return;
              }
           }
           plhs[0] = mxCreateLogicalScalar(c > 0);
     
        } else {
           BadRelopError();
        }
        break;
        
     case 2:
        if (!memcmp(Relop, u"<=", 2 * sizeof(mxChar))) {         // MatlabV <= N
           for (i = 0; i < nData; i++) {
              if ((c = MatlabV[i] - (int) Data[i]) != 0) {
                 plhs[0] = mxCreateLogicalScalar(c < 0);
                 return;
              }
           }
           plhs[0] = mxCreateLogicalScalar(c <= 0);
           
        } else if (!memcmp(Relop, u">=", 2 * sizeof(mxChar))) {  // MatlabV >= N
           for (i = 0; i < nData; i++) {
              if ((c = MatlabV[i] - (int) Data[i]) != 0) {
                 plhs[0] = mxCreateLogicalScalar(c > 0);
                 return;
              }
           }
           plhs[0] = mxCreateLogicalScalar(c >= 0);
           
        } else if (!memcmp(Relop, u"==", 2 * sizeof(mxChar))) {  // MatlabV == N
           for (i = 0; i < nData; i++) {
              if (MatlabV[i] != (int) Data[i]) {
                 plhs[0] = mxCreateLogicalScalar(false);
                 return;
              }
           }
           plhs[0] = mxCreateLogicalScalar(true);
           
        } else {
           BadRelopError();
        }
        
        break;
        
     default:
         BadRelopError();
  }
  
  return;
}

// =============================================================================
void Init(void)
{
  mxArray *Arg[3], *Ver[1];
  double  *VerP;

  // Call "version", store result in Arg[0]:
  mexCallMATLAB(1, Arg, 0, NULL, "version");
  
  // Call SSCANF(version, '%d.', 4), reuse Arg[0]:
  Arg[1] = mxCreateString("%d.");
  Arg[2] = mxCreateDoubleScalar(4);
  mexCallMATLAB(1, Ver, 3, Arg, "sscanf");
  
  // Store result locally:
  VerP       = mxGetPr(Ver[0]);
  MatlabV[0] = (int) VerP[0];
  MatlabV[1] = (int) VerP[1];
  MatlabV[2] = (int) VerP[2];
  MatlabV[3] = (int) VerP[3];
  
  // Cleanup:
  mxDestroyArray(Arg[0]);
  mxDestroyArray(Arg[1]);
  mxDestroyArray(Arg[2]);
  mxDestroyArray(Ver[0]);
  
  return;
}

// =============================================================================
void BadRelopError(void)
{
  mexErrMsgIdAndTxt("JSim:isMatlabVer:BadRelop",
                    "Operator must be: '<=', '<', '>=', '>', '=='.");
}

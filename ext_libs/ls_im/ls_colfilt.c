/**
 * imfilter() replacement
 * Author: Stanislaw Adaszewski, 2012
 */
 
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const mxArray *A, *B, *C;
    const mwSize *A_dim, *B_dim;
    mwSize x, y;
    mwSize A_ofs = 0;
    double *A_dat, *B_dat;
    mwSize M, N;
    mwSize m, n;
    mxArray *R;
    double *R_dat;
    mwSize X, Y;
    char *funcname;
    mwSize T_dim[2];
    mxArray *T;
    double *T_dat;
    
    if (nlhs != 1 || nrhs != 3) {
        mexErrMsgTxt("Usage: R=ls_im2col(A, B, C), where:\n\nA - input image,\nB - block size,\nC - function name,\nR - result.");
        return;
    }
    
    A = prhs[0];
    B = prhs[1];
    C = prhs[2];
    
    if (!mxIsDouble(A) || mxIsComplex(A) || mxIsSparse(A)) {
        mexErrMsgTxt("A has to be float, real, full.");
        return;
    }
    
    if (!mxIsDouble(B) || mxIsComplex(B) || mxIsSparse(B)) {
        mexErrMsgTxt("B has to be float, real, full.");
        return;
    }
    
    if (!mxIsChar(C) || mxIsEmpty(C)) {
        mexErrMsgTxt("C has to be the name of a function.");
        return;
    }
    
    A_dim = mxGetDimensions(A);
    B_dim = mxGetDimensions(B);
    
    if (B_dim[0] < 2 && B_dim[1] < 2) {
        mexErrMsgTxt("B has to contain at least two elements.");
        return;
    }
    
    funcname = (char*) malloc(mxGetM(C) * mxGetN(C) * sizeof(mxChar) + 1);
    mxGetString(C, funcname, mxGetM(C) * mxGetN(C) * sizeof(mxChar) + 1);
    
    A_dat = mxGetPr(A);
    B_dat = mxGetPr(B);
    
    X = (mwSize)A_dim[0];
    Y = (mwSize)A_dim[1];
    
    M = (mwSize)B_dat[0];
    N = (mwSize)B_dat[1];
    
    T_dim[0] = M * N;
    T_dim[1] = 1;
    T = mxCreateNumericArray(2, T_dim, mxDOUBLE_CLASS, mxREAL);
    T_dat = mxGetPr(T);
   
    R = mxCreateNumericArray(2, A_dim, mxDOUBLE_CLASS, mxREAL);
    R_dat = mxGetPr(R);
   
    A_ofs = 0;
    for (y = 0; y < A_dim[1]; y++)
        for (x = 0; x < A_dim[0]; x++) {
            mwSize cnt = 0;
            mxArray *lhs[] = {0};
            mxArray *rhs[] = {T};
            mxArray *exc;
            for (n = 0; n < N; n++)
                for (m = 0; m < M; m++) {
                    mwSize x1 = x + m - (M >> 1);
                    mwSize y1 = y + n - (N >> 1);
                    if (x1 >= 0 && x1 < A_dim[0] && y1 >= 0 && y1 < A_dim[1]) {
                        mwSize A_ofs1 = y1 * A_dim[0] + x1;
                        T_dat[cnt] = A_dat[A_ofs1];
                    } else {
                        T_dat[cnt] = 0;
                    }
                    cnt++;
                }
            
            exc = mexCallMATLABWithTrap(1, lhs, 1, rhs, funcname);
            if (exc) {
                mxDestroyArray(exc);
                mexErrMsgTxt("Error calling specified function.");
                mxDestroyArray(T);
                mxDestroyArray(R);
                free(funcname);
                return;
            } else {
                if (!mxIsEmpty(lhs[0]) && mxIsDouble(lhs[0]) && !mxIsComplex(lhs[0]) && !mxIsSparse(lhs[0])) {
                    R_dat[A_ofs] = mxGetPr(lhs[0])[0];
                } else {
                    R_dat[A_ofs] = 0;
                }
                mxDestroyArray(lhs[0]);
            }
            
            // mexPrintf("Here!\n");
            
            A_ofs++;
        }
    
    mxDestroyArray(T);    
    free(funcname);
    plhs[0] = R;
}
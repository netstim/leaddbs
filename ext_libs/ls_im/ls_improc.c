/**
 * imfilter() replacement
 * Author: Stanislaw Adaszewski, 2012
 */

#include "mex.h"
#include <float.h>
#include <pthread.h>

#define BM_FIXEDVAL 0
#define BM_SYMMETRIC 1

#define MODE_DILATE 0
#define MODE_ERODE 1
#define MODE_FILTER 2

typedef int mwSize;

typedef struct ls_improc_t_ {
	int mode;
	const mxArray *A, *B;
	mxArray *R;
	double *A_dat, *B_dat, *R_dat;
	mwSize A_ndim, B_ndim;
	const mwSize *A_dim, *B_dim;
	mwSize A_dim_x, A_dim_y, A_dim_z, B_dim_x, B_dim_y, B_dim_z;
	mwSize A_ofs;
	mwSize n_workers;
    int boundary_mode;
    double boundary_val;
} ls_improc_t;

typedef struct ls_improc_worker_t_ {
    ls_improc_t *g;
    mwSize n;
} ls_improc_worker_t;

void* worker(void *arg) {
	mwSize A_x, A_y, A_z, B_x, B_y, B_z;
	
	ls_improc_worker_t *w = (ls_improc_worker_t*) arg;
	ls_improc_t *g = w->g;
	mwSize n = w->n;
	
	mwSize start = n * g->A_dim_x / g->n_workers;
	mwSize stop = (n + 1) * g->A_dim_x / g->n_workers;
	
	if (stop > g->A_dim_x)
	stop = g->A_dim_x;

	for (A_z = 0; A_z < g->A_dim_z; A_z++)
		for (A_y = 0; A_y < g->A_dim_y; A_y++)
			for (A_x = start; A_x < stop; A_x++) {
				mwSize B_ofs = 0;
				mwSize A_ofs = A_x + g->A_dim_x * (A_y + g->A_dim_y * A_z);
				double sum = 0;
				for (B_z = 0; B_z < g->B_dim_z; B_z++)
					for (B_y = 0; B_y < g->B_dim_y; B_y++)
						for (B_x = 0; B_x < g->B_dim_x; B_x++) {
							mwSize R_x = A_x + B_x - (g->B_dim_x >> 1);
							mwSize R_y = A_y + B_y - (g->B_dim_y >> 1);
							mwSize R_z = A_z + B_z - (g->B_dim_z >> 1);
                            double B_val = g->B_dat[B_ofs];
                            
                            if (g->boundary_mode == BM_SYMMETRIC) {
                                if (R_x < 0)
                                    R_x += g->A_dim_x;
                                
                                if (R_y < 0)
                                    R_y += g->A_dim_y;
                                
                                if (R_z < 0)
                                    R_z += g->A_dim_z;
                                
                                if (R_x >= g->A_dim_x)
                                    R_x -= g->A_dim_x;
                                
                                if (R_y >= g->A_dim_y)
                                    R_y -= g->A_dim_y;
                                
                                if (R_z >= g->A_dim_z)
                                    R_z -= g->A_dim_z;
                            }
			
							if (R_x >= 0 && R_x < g->A_dim_x &&
								R_y >= 0 && R_y < g->A_dim_y &&
								R_z >= 0 && R_z < g->A_dim_z) {
                                mwSize R_ofs = R_x + g->A_dim_x * (R_y + g->A_dim_y * R_z);
                                if (g->mode == MODE_FILTER) {
                                    sum += g->A_dat[R_ofs] * B_val;
                                } else {
                                    double p = g->A_dat[A_ofs] * B_val;
                                    if (B_val == 1 && ((p > g->R_dat[R_ofs]) ^ g->mode)) {
                                        g->R_dat[R_ofs] = p;
                                    }
                                }
                            } else if (g->mode == MODE_FILTER) {
                                sum += g->boundary_val * B_val;
                            }
									
							B_ofs++;
						}
			
				if (g->mode == MODE_FILTER)
					g->R_dat[A_ofs] = sum;
			}

	return NULL;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	mwSize n;
	ls_improc_t g;

	if (nrhs < 4 || nlhs != 1) {
		mexErrMsgTxt("Usage: R = ls_improc(A, B, mode, nthreads, boundary), where:\n\nA is the input image,\nB is the structural element for morphological operations or a filter,\nmode is 0, 1, 2 for dilation, erosion and filter respectively,\nnthreads specifies number of threads used,\nboundary is boundary mode for filter,\nR is is the result image.");
		return;
	}

	if (mxIsEmpty(prhs[2]) || !mxIsDouble(prhs[2])) {
		mexErrMsgTxt("erode parameter should be 0 or 1, real");
		return;
	}

	if (mxIsEmpty(prhs[3]) || !mxIsDouble(prhs[3])) {
		mexErrMsgTxt("nthreads parameter should be >=1, real, float");
	}
    
    if (nrhs > 4 && !mxIsEmpty(prhs[4])) {
        if (mxIsDouble(prhs[4])) {
            g.boundary_mode = BM_FIXEDVAL;
            g.boundary_val = mxGetPr(prhs[4])[0];
        } else if (mxIsChar(prhs[4])) {
            mwSize buflen = mxGetM(prhs[4]) * mxGetN(prhs[4]) * sizeof(mxChar) + 1;
            char *buf = (char*) malloc(buflen);
            mxGetString(prhs[4], buf, buflen);
            if (!strcmp(buf, "symmetric")) {
                g.boundary_mode = BM_SYMMETRIC;
                g.boundary_val = 0;
            } else {
                mexErrMsgTxt("Unsupported boundary mode");
                return;
            }
            free(buf);
        } else {
            mxErrMsgTxt("boundary has to be either float or char");
        }
    } else {
        g.boundary_mode = BM_FIXEDVAL;
        g.boundary_val = 0;
    }
    
    // mexPrintf("boundary_mode: %d, boundary_val: %f\n", g.boundary_mode, g.boundary_val);
        
	g.mode = (int) mxGetPr(prhs[2])[0];
	g.n_workers = (mwSize) mxGetPr(prhs[3])[0];

	g.A = prhs[0];
	g.B = prhs[1];

	if (!mxIsDouble(g.A) || !mxIsDouble(g.B)) {
		mexErrMsgTxt("Parameters to ls_imdilate have to be real matrices.");
		return;
	}

	g.A_dat = mxGetPr(g.A);
	g.B_dat = mxGetPr(g.B);

	g.A_ndim = mxGetNumberOfDimensions(g.A);
	g.B_ndim = mxGetNumberOfDimensions(g.B);

	if (g.A_ndim > 3) {
		mexErrMsgTxt("Only up to 3 dimensions are supported.");
		return;
	}

	g.A_dim = mxGetDimensions(g.A);
	g.B_dim = mxGetDimensions(g.B);

	g.A_dim_x = g.A_dim[0];
	g.A_dim_y = g.A_dim[1];

	g.B_dim_x = g.B_dim[0];
	g.B_dim_y = g.B_dim[1];

	if (g.A_ndim == 3)
		g.A_dim_z = g.A_dim[2];
	else
		g.A_dim_z = 1;

	if (g.B_ndim == 3)
		g.B_dim_z = g.B_dim[2];
	else
		g.B_dim_z = 1;

	if (g.mode < 2) {
		for (n = 0; n < g.B_dim_x * g.B_dim_y * g.B_dim_z; n++)
			if (g.B_dat[n] != 0 && g.B_dat[n] != 1) {
				mexErrMsgTxt("Structural element must contain only 0s and 1s.");
				return;
			}
	}

	g.R = mxCreateNumericArray(g.A_ndim, g.A_dim, mxDOUBLE_CLASS, mxREAL);

	g.R_dat = mxGetPr(g.R);

	if (g.mode < 2) {
		for (n = 0; n < g.A_dim_x * g.A_dim_y * g.A_dim_z; n++)
			g.R_dat[n] = g.mode ? DBL_MAX : -DBL_MAX;
	}

	plhs[0] = g.R;

	if (g.n_workers < 1)
		g.n_workers = 1;
	{
		pthread_t *thr = malloc(sizeof(pthread_t) * g.n_workers);
		ls_improc_worker_t *w = malloc(sizeof(ls_improc_worker_t) * g.n_workers);
	
		for (n = 0; n < g.n_workers; n++) {
			w[n].g = &g;
			w[n].n = n;
			pthread_create(&thr[n], NULL, worker, &w[n]);
		}
		for (n = 0; n < g.n_workers; n++) {
			pthread_join(thr[n], NULL);
		}
	
		free(thr);
		free(w);
	}
}

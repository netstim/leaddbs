
#include "auxilary_classes.cpp"

class SphereInterpolator
{
public:
    REAL *barycoords;
    int *indices;
    int size;  // size of LUT
    int sN;   // (sizeofLUT-1)/2
    int nverts; // number of data vertices

    REAL *dirs;

    REAL beta;
    REAL inva;
    REAL b;
    
    
    int *idx;
    REAL *interpw;

    
    
    SphereInterpolator(const mxArray *sinterpstruct)
    {
        mxArray *Indices = mxGetField(sinterpstruct,0,"indices");
        mxArray *BaryCoords = mxGetField(sinterpstruct,0,"barycoords");
        mxArray *Beta = mxGetField(sinterpstruct,0,"beta");
        mxArray *NumInterpPoints = mxGetField(sinterpstruct,0,"numpoints");
        mxArray *Dirs = mxGetField(sinterpstruct,0,"bDir");	
        REAL *indimg = (REAL*) mxGetData(Indices);
        const int *isize = mxGetDimensions(Indices);
        int totsz = isize[0]*isize[1]*isize[2]*isize[3];
        int *indeximg = (int*) malloc(sizeof(int)*totsz);
        for (int k =0;k < totsz;k++)
            indeximg[k] = int(indimg[k])-1;      
        REAL *barycoords = (REAL*) mxGetData(BaryCoords);
        REAL *beta = (REAL*) mxGetData(Beta);
        int nip = int(*((REAL*)mxGetData(NumInterpPoints)));	
        
        init(barycoords,indeximg,nip,isize[2],beta[0],(REAL*) mxGetData(Dirs));
    }
    
    ~SphereInterpolator()
    {
        free(indices);
    }
    
    void init(REAL *barycoords, int *indices, int numverts, int sizeLUT, REAL beta,REAL *dirs)
    {
        
        this->barycoords = barycoords;
        this->indices = indices;
        this->size = sizeLUT;
        this->sN = (sizeLUT-1)/2;
        this->nverts = numverts;
        this->beta = beta;
        this->dirs = dirs;

        inva = (sqrt(1+beta)-sqrt(beta));
        b = 1/(1-sqrt(1/beta + 1));    
       
    }

    
    inline void getInterpolation(pVector &N)
    {
	    REAL nx = N.x;
	    REAL ny = N.y;
	    REAL nz = N.z;

            if (nz > 0.5)
            {
                int x = real2int(nx);
                int y = real2int(ny);
                int i = 3*6*(x+y*size);  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nz < -0.5)
            {
                int x = real2int(nx);
                int y = real2int(ny);
                int i = 3*(1+6*(x+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nx > 0.5)
            {
                int z = real2int(nz);
                int y = real2int(ny);
                int i = 3*(2+6*(z+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (nx < -0.5)
            {
                int z = real2int(nz);
                int y = real2int(ny);
                int i = 3*(3+6*(z+y*size));  // (:,2,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }
            if (ny > 0)
            {
                int x = real2int(nx);
                int z = real2int(nz);
                int i = 3*(4+6*(x+z*size));  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }     
            else
            {
                int x = real2int(nx);
                int z = real2int(nz);
                int i = 3*(5+6*(x+z*size));  // (:,1,x,y)
                idx = indices+i;  
                interpw = barycoords +i; 
                return;
            }            
                        
    }
    
    
    inline REAL invrescale(REAL f)
    {
        REAL x = (fabs(f)-b)*inva;
        if (f>0)            
            return (x*x-beta);
        else
            return beta - x*x;
    }
    
    inline int real2int(REAL x)
    {
       return int((invrescale(x)+1)*sN-0.5);
    
    }

    
};




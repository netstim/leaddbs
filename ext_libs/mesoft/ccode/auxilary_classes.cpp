
#ifndef _AUXCLASS
#define _AUXCLASS

#ifndef INFINITY
#define INFINITY 1000000000
#endif

//#define INLINE __attribute__ ((always_inline))
#define INLINE inline

//#define __SSE2


class ParamStruct
{
public:
    ParamStruct()
    {}

	REAL Temp;
	REAL numit;
	REAL gammaDE;
    REAL GM_mask;
	REAL particle_weight;
	REAL particle_width;
	REAL particle_len;
	REAL particle_expnump;
	REAL chempot_connection;
    REAL kappa_connection;
    REAL bound_attract;
	REAL chempot_particle;
	REAL num_qspace_samples;	
	REAL meanval_sq;
    REAL numcores;
    REAL voxelSmooth_vf;
    REAL voxelSmooth_Ds;
    REAL voxelSmooth_Ds_range;
    REAL fiberSmooth;
    REAL conpot_lambda;
    REAL data_lambda;
    REAL penalty_SW;    
    REAL trace_equality_constraint;
    REAL curv_hardthres;
    REAL trackingguide_strength;
    REAL directional_propsal_distrib;
    REAL restrictions;
    REAL prop_p_birth;
    REAL prop_p_death;		
    REAL prop_p_shift;
    REAL prop_p_Dmod;     
    REAL prop_p_vfmod; 
    REAL prop_p_conprob;
    REAL numbvals;
    REAL lmax;
    
    
    void readParameters(const mxArray *mxarr)
    {
        Temp = read("currentTemp",mxarr);
        numit = read("maxit",mxarr);
        gammaDE = read("lam_DE",mxarr);
        GM_mask = read("GM_mask",mxarr);
        particle_weight = read("p_weight",mxarr); 
        particle_width = read("p_wid",mxarr); particle_width = 1/particle_width;
        particle_len = read("p_len",mxarr);
        chempot_connection = read("c_likli",mxarr);
        kappa_connection = read("c_kappa",mxarr);
        bound_attract = read("c_bound_attract",mxarr);
        particle_expnump = read("p_expnump",mxarr);
        chempot_particle = read("p_chempot",mxarr);
        num_qspace_samples = read("num_qspace_samples",mxarr);             
        meanval_sq = read("alpha",mxarr);
        numcores = read("numcores",mxarr);
        voxelSmooth_vf = read("voxelSmooth_vf",mxarr);
        voxelSmooth_Ds = read("voxelSmooth_Ds",mxarr);
        voxelSmooth_Ds_range = read("voxelSmooth_Ds_range",mxarr);
        fiberSmooth = read("fiberSmooth",mxarr);
        conpot_lambda = read("conpot_lambda",mxarr);
        data_lambda = read("data_lambda",mxarr);
        penalty_SW = read("penalty_SW",mxarr);
        trace_equality_constraint = read("trace_equality_constraint",mxarr);
        curv_hardthres = read("curv_hardthres",mxarr);
        trackingguide_strength= read("trackingguide_strength",mxarr);
        directional_propsal_distrib= read("directional_propsal_distrib",mxarr);
        restrictions= read("restrictions",mxarr);
        prop_p_birth = read("prop_p_birth",mxarr);
        prop_p_death = read("prop_p_death",mxarr);		
        prop_p_shift = read("prop_p_shift",mxarr);
        prop_p_Dmod = read("prop_p_Dmod",mxarr);     
        prop_p_vfmod = read("prop_p_vfmod",mxarr); 
        prop_p_conprob = read("prop_p_conprob",mxarr);
        lmax = read("lmax",mxarr); 
        numbvals = read("numbvals",mxarr);

    }
    
    REAL read(const char *name, const mxArray *mxarr)
    {
        mxArray *data = mxGetField(mxarr,0,name);
        REAL val = *mxGetPr(data);
        fprintf(stderr,"Parameter %s set to %f\n",name,val);
        return val;
    }
    
    
   
};








class pVector
{
       
    public:

    REAL x;
    REAL y;
    REAL z;

    pVector()
    {
	
    }    

    INLINE void setXYZ(REAL sx,REAL sy, REAL sz) 
    {
	x = sx;
	y = sy;
	z = sz;
    }

    
    pVector(REAL x,REAL y,REAL z)
    {
       this->x = x;
       this->y = y;       
       this->z = z;
    }
    
    INLINE void rand(REAL w,REAL h,REAL d)
    {
        this->x = mtrand.frand()*w;
        this->y = mtrand.frand()*h;
        this->z = mtrand.frand()*d;

    }
    
    INLINE void rand_sphere()
    {
        this->x = mtrand.frandn();
        this->y = mtrand.frandn();
        this->z = mtrand.frandn();
        normalize();
    }
    
    INLINE void normalize()
    {
        REAL norm = sqrt(x*x+y*y+z*z)+ 0.00000001;
        *this /= norm;
       
    }
    
    INLINE REAL norm_square()
    {
        return x*x + y*y + z*z;
    }
    
    INLINE void distortn(REAL sigma)
    {
        x += sigma*mtrand.frandn();
        y += sigma*mtrand.frandn();
        z += sigma*mtrand.frandn();
    }
    
   
    INLINE pVector operator*(REAL s)
    {
        return pVector(s*x,s*y,s*z);
    }
    
    INLINE void operator*=(REAL &s)
    {
         x *= s;
         y *= s;
         z *= s;
    }
    
    INLINE pVector operator+(pVector R)
    {
        return pVector(x+R.x,y+R.y,z+R.z);
    }
    
    INLINE void operator+=(pVector R)
    {
         x += R.x;
         y += R.y;
         z += R.z;
    }
    
    INLINE pVector operator-(pVector R)
    {
        return pVector(x-R.x,y-R.y,z-R.z);
    }

    INLINE void operator-=(pVector R)
    {
         x -= R.x;
         y -= R.y;
         z -= R.z;
    }
    
    INLINE pVector operator/(REAL &s)
    {
        return pVector(x/s,y/s,z/s);
    }
    
    INLINE void operator/=(REAL &s)
    {
         x /= s;
         y /= s;
         z /= s;
    }
    
    INLINE bool operator==(pVector R)
    {
        return (x == R.x) && (y == R.y) && (z == R.z);
    }
       
    INLINE REAL operator*(pVector R)
    {
        return x*R.x+y*R.y+z*R.z;
    }
    
    
};


//----------------- the Particle
class Particle 
{
	public:

	Particle()
	{
		pID = -1;
		mID = -1;
	}

	pVector R; // position
	pVector N; // direction
	REAL len;  // length (fixed, so basically unsued)
    
    REAL Di,Da,Dp,w,q; // diffusion parameters and weights

	int gridindex[3]; // index in the grids where it is living	
	int ID;
	int pID;
	int mID;
    
    bool active; // is this a real particle, or just an empty slot in the particle list
    int label;
};



//------------ a class to draw a object from a numerically given distribution (by bisectional search)
template <class T>
class SimpSamp
{
	
	
public:
	T *objs;
	int cnt;
	REAL *P;


	SimpSamp(int maxobj)
	{
		P = (REAL*) malloc(sizeof(REAL)*(maxobj+1));
		objs = (T*) malloc(sizeof(T)*(maxobj+1));
        clear();
	}
	~SimpSamp()
	{
		free(P);
		free(objs);
	}

	INLINE void clear()
	{
		cnt = 1;
		P[0] = 0;
	}

	INLINE void add(REAL p, T obj)
	{
		P[cnt] = P[cnt-1] + p;
		objs[cnt-1] = obj;
		cnt++;
	}


	INLINE int draw()
	{
		REAL r = mtrand.frand()*P[cnt-1];		
		int j;
		int rl = 1;
		int rh = cnt-1;
		while(rh != rl)
		{
			j = rl + (rh-rl)/2;
			if (r < P[j])
				{
				rh = j;
				continue;
				}
			if (r > P[j])
				{
				rl = j+1;
				continue;
				}
			break;
		}	
		return rh-1;
	}





	INLINE T drawObj()
	{
		return objs[draw()];
	}

	INLINE bool isempty()
	{
		if (cnt == 1)
			return true;
		else
			return false;
	}


	REAL probFor(int idx)
	{
		return (P[idx+1]-P[idx])/P[cnt-1];
	}

	REAL probFor(T &t)
	{
		for (int i = 1; i< cnt;i++)
		{
			if (t == objs[i-1])
				return probFor(i-1);
		}
		return 0;
	}



};


class EndPoint 
{
	public:
	EndPoint()
	{}

	EndPoint(Particle *p,int ep)
	{
		this->p = p;
		this->ep = ep;
	}
	Particle *p;
	int ep; 

	inline bool operator==(EndPoint P)
	{
		return (P.p == p) && (P.ep == ep);
	}
};


#endif



/***************************************************************************************************************
    Copyright (c) 2011, Marco Reisert, Valerij G. Kiselev, Medical Physics, University Medical Center Freiburg
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
	* Redistributions of source code must retain the above copyright
	  notice, this list of conditions and the following disclaimer.
	* Redistributions in binary form must reproduce the above copyright
	  notice, this list of conditions and the following disclaimer in the
	  documentation and/or other materials provided with the distribution.
	* Neither the name of the <organization> nor the
	  names of its contributors may be used to endorse or promote products
	  derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
    ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
****************************************************************************************************************/

#ifndef _AUXCLASS
#define _AUXCLASS

#ifndef INFINITY
#define INFINITY 1000000000
#endif

#define REAL float
#define PI 3.1415926536

//#define INLINE __attribute__ ((always_inline))
#define INLINE inline

//#define __SSE2



#ifdef __SSE2
#include <xmmintrin.h>

class pVector
{

    private:
    __m128 r;


    public:

    static REAL store[4];    


    pVector()
    {
	

    }
    
    pVector(__m128 q)
    {
	r = q;
    }

    
    pVector(REAL x,REAL y,REAL z)
    {
	r = _mm_set_ps(0,z,y,x);
    }


    INLINE void storeXYZ()
    {
	_mm_store_ps(store,r);
    }
      
    INLINE void setXYZ(REAL sx,REAL sy,REAL sz)
    {
	r = _mm_set_ps(0,sz,sy,sx);
    }   

   

    
    INLINE void rand(REAL w,REAL h,REAL d)
    {
        REAL x = mtrand.frand()*w;
        REAL y = mtrand.frand()*h;
        REAL z = mtrand.frand()*d;
	r = _mm_set_ps(0,z,y,x);

    }
    
    INLINE void rand_sphere()
    {
	r = _mm_set_ps(0,mtrand.frandn(),mtrand.frandn(),mtrand.frandn());
        normalize();
    }

    INLINE void normalize()
    {
	__m128 q = _mm_mul_ps(r,r);
	_mm_store_ps(store,q);
        REAL norm = sqrt(store[0]+store[1]+store[2]) + 0.00000001;
	q = _mm_set_ps1(1/norm);
	r = _mm_mul_ps(r,q);
       
    }

    INLINE REAL norm_square()
    {
	__m128 q = _mm_mul_ps(r,r);
	_mm_store_ps(store,q);
        return store[0]+store[1]+store[2];
    }

    INLINE void distortn(REAL sigma)
    {
	__m128 s = _mm_set_ps(0,mtrand.frandn(),mtrand.frandn(),mtrand.frandn());
	__m128 q = _mm_set_ps1(sigma);
	r = _mm_add_ps(r,_mm_mul_ps(q,s));
    }

    INLINE pVector operator*(REAL s)
    {
	 __m128 q = _mm_set_ps1(s);
	 return pVector(_mm_mul_ps(q,r));
    }
    
    INLINE void operator*=(REAL &s)
    {
	 __m128 q = _mm_set_ps1(s);
	 r = _mm_mul_ps(q,r);
    }
    
    INLINE pVector operator+(pVector R)
    {
        return pVector(_mm_add_ps(R.r,r));
    }
    
    INLINE void operator+=(pVector R)
    {
         r = _mm_add_ps(r,R.r);
    }
    
    INLINE pVector operator-(pVector R)
    {
        return pVector(_mm_sub_ps(r,R.r));
    }

    INLINE void operator-=(pVector R)
    {
         r = _mm_sub_ps(r,R.r);
    }
    
    INLINE pVector operator/(REAL &s)
    {
	 __m128 q = _mm_set_ps1(s);
	 return pVector(_mm_div_ps(q,r));
    }
    
    INLINE void operator/=(REAL &s)
    {
	 __m128 q = _mm_set_ps1(s);
         r = _mm_div_ps(r,q);
    }
    
    INLINE REAL operator*(pVector R)
    {
	__m128 q = _mm_mul_ps(r,R.r);
	_mm_store_ps(store,q);
        return store[0]+store[1]+store[2];
    }
    
    
};

REAL pVector::store[4];    

#else

class pVector
{

    private:
    REAL x;
    REAL y;
    REAL z;
    
    public:

    static REAL store[4];    

    pVector()
    {
	

    }
    

    INLINE void storeXYZ()
    {
	store[0] = x;
	store[1] = y;
	store[2] = z;
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
    
       
    INLINE REAL operator*(pVector R)
    {
        return x*R.x+y*R.y+z*R.z;
    }
    
    
};

REAL pVector::store[4];    

#endif


class Particle 
{
	public:

	Particle()
	{
		label = 0;
		pID = -1;
		mID = -1;
	}

	pVector R;
	pVector N;
	REAL cap;
	REAL len;

	int gridindex; // index in the grid where it is living	
	int ID;
	int pID;
	int mID;

	int label; 
};


class EnergyGradient
{
    public:
        pVector gR;
        pVector gN;

    INLINE REAL norm2()
    {
        return gR.norm_square() + gN.norm_square();
    }
} ;


template <class T>
class SimpSamp
{
	
	REAL *P;
	int cnt;
	
public:
	T *objs;


	SimpSamp()
	{
		P = (REAL*) malloc(sizeof(REAL)*1000);
		objs = (T*) malloc(sizeof(T)*1000);
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

// 	INLINE int draw()
// 	{
// 		REAL r = mtrand.frand()*P[cnt-1];
// 		for (int i = 1; i < cnt; i++)
// 		{
// 			if (r <= P[i])
// 				return i-1; 
// 		}
// 		return cnt-2;
// 	}

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

class Track
{
	public:
	EndPoint track[1000];
	REAL energy;
	REAL proposal_probability;
	int length;

	void clear()
	{
		length = 0;
		energy = 0;
		proposal_probability = 1;
	}


	bool isequal(Track &t)
	{
		for (int i = 0; i < length;i++)
		{
			if (track[i].p != t.track[i].p || track[i].ep != t.track[i].ep)
				return false;
		}
		return true;
	}

};

REAL getMax(REAL *arr, int cnt)
{
	REAL max = arr[0];
	for (int i = 1; i < cnt; i++)
	{
		if (arr[i] > max) 
			max = arr[i];
	}
	return max;
}



REAL getMin(REAL *arr, int cnt)
{
	REAL min = arr[0];
	for (int i = 1; i < cnt; i++)
	{
		if (arr[i] < min) 
			min = arr[i];
	}
	return min;
}


int getArgMin(REAL *arr, int cnt)
{
	REAL min = arr[0];
	int idx = 0;
	for (int i = 1; i < cnt; i++)
	{
		if (arr[i] < min) 
		{
			min = arr[i];
			idx = i;
		}
	}
	return idx;
}


REAL PosOrientGaussDistrib(REAL dR, REAL C,REAL sigma1, REAL sigma2)
{
    REAL p2 =
    exp(-1/(2*sigma2))*pow(1/(sigma2*PI),1.5)/(4*sqrt(2.0)) *
         ( 2*C*sigma1 + exp(C*C/2*sigma2)*sqrt(2*PI)*sqrt(sigma2)*(C*C+sigma2)*(1+erff((float) C/sqrt(2*sigma2))));
    REAL p1 = sqrt(1/(2*PI*sigma1)) * exp(-dR/(2*sigma1));
                     
    return p1*p2;    
}





inline REAL distLseg(pVector &R1,pVector &N1,pVector &R2,pVector &N2,REAL &len)
{

	pVector D = R1-R2;
	REAL beta  = N1*N2;
	REAL divisor = 1.001-beta*beta;
	REAL gamma1 = N1*D;
	REAL gamma2 = N2*D;
	REAL t,u;
	REAL EPdist[4];

	pVector Q;
	REAL dist = 102400000.0;

	while(true)
	{

		t = -(gamma1+beta*gamma2) / divisor;
		u =  (gamma1*beta+gamma2) / divisor;
		if (fabs(t) < len && fabs(u) < len)
		{
			Q = D +  N1*t - N2*u;
			dist = Q*Q;
			break;
		}

		beta = len*beta;

		t = beta - gamma1;
		if (fabs(t) < len) 
		{
			Q = D +  N1*t - N2*len;
			REAL d = Q*Q;
			if (d < dist) dist = d;
		}

		t = -beta - gamma1;
		if (fabs(t) < len) 
		{
			Q = D +  N1*t + N2*len;
			REAL d = Q*Q;
			if (d < dist) dist = d;
		}

		u = beta + gamma2;
		if (fabs(u) < len)
		{
			Q = D +  N1*len - N2*u;
			REAL d = Q*Q;
			if (d < dist) dist = d;
		}					

		u = -beta + gamma2;
		if (fabs(u) < len)
		{
			Q = D -  N1*len - N2*u;
			REAL d = Q*Q;
			if (d < dist) dist = d;
		}					

		if (dist != 102400000.0)
			break;


		EPdist[0] =  beta + gamma1 - gamma2;
		EPdist[1] = -beta + gamma1 + gamma2;
		EPdist[2] = -beta - gamma1 - gamma2;
		EPdist[3] =  beta - gamma1 + gamma2;
		int c = getArgMin(EPdist,4);
		if (c==0) {t = +len; u = +len; }
		if (c==1) {t = +len; u = -len; }
		if (c==2) {t = -len; u = +len; }
		if (c==3) {t = -len; u = -len; }
		Q = D +  N1*t - N2*u;
		dist = Q*Q;
		break;

	}


	return dist; 

}


#endif



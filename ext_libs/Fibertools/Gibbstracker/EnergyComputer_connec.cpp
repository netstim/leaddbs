/***************************************************************************************************************
    Copyright (c) 2011, Marco Reisert, Valerij G. Kiselev, Medical Physics, University Medical Center Freiburg
    All rights reserved.

    Reference: Reisert M, Mader I, Anastasopoulos C, Weigel M, Schnell S, Kiselev V. 
    Global fiber reconstruction becomes practical.
    Neuroimage.  2011 Jan 15;54(2):955-62. Epub 2010 Sep 18.


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



/// bessel for wid = 1

//#define mbesseli0(x) ((x>=1.0) ? ((-0.2578*x+0.7236)*exp(x)) : ((x>=0.9) ? ((-0.2740*x+0.7398)*exp(x)) : ((x>=0.8) ? ((-0.3099*x+0.7720)*exp(x)) : ((x>=0.7) ? ((-0.3634*x+0.8149)*exp(x)) : ((x>=0.5) ? ((-0.4425*x+0.8663)*exp(x)) : ((x>=0.3) ? ((-0.5627*x+0.9264)*exp(x)) : ((x>=0.2) ? ((-0.6936*x+0.9657)*exp(x)) : ((x>=0.1) ? ((-0.8016*x+0.9873)*exp(x)) : ((x>=0.0) ? ((-0.9290*x+1.0000)*exp(x)) : 1 )))))))))
//#define mbesseli0(x) ((x>=1.0) ? ((0.5652*x+0.7009)) : ((x>=0.8) ? ((0.4978*x+0.7683)) : ((x>=0.6) ? ((0.3723*x+0.8686)) : ((x>=0.4) ? ((0.2582*x+0.9371)) : ((x>=0.2) ? ((0.1519*x+0.9796)) : ((x>=0.0) ? ((0.0501*x+1.0000)) : 1 ))))))
// // 

inline REAL mbesseli0(REAL x)
{
    REAL y = x*x;
    REAL erg = BESSEL_APPROXCOEFF[0]; 
    erg += y*BESSEL_APPROXCOEFF[1];
    erg += y*y*BESSEL_APPROXCOEFF[2]; 
    erg += y*y*y*BESSEL_APPROXCOEFF[3];
    return erg;
}


// 
// 
// inline REAL mbesseli0(REAL x)
// {
//     REAL y = x*x;
//     REAL erg = BESSEL_APPROXCOEFF[0]; 
//     erg += y*BESSEL_APPROXCOEFF[1];
//     erg += y*x*BESSEL_APPROXCOEFF[2]; 
//     erg += x*x*BESSEL_APPROXCOEFF[3];
//     return erg;
// }


inline REAL mexp(REAL x)
{

	return((x>=7.0) ? 0 : ((x>=5.0) ? (-0.0029*x+0.0213) : ((x>=3.0) ? (-0.0215*x+0.1144) : ((x>=2.0) ? (-0.0855*x+0.3064) : ((x>=1.0) ? (-0.2325*x+0.6004) : ((x>=0.5) ? (-0.4773*x+0.8452) : ((x>=0.0) ? (-0.7869*x+1.0000) : 1 )))))));
//	return exp(-x);

}


inline REAL mdoty(REAL x)
{

	return ((x>=0.98) ? 1 : 0);
//	return exp(-x);

}
#include "ParticleGrid.cpp"

#include "EnergyComputerBase.cpp"


class EnergyComputer : public EnergyComputerBase
{

public:


	REAL eigencon_energy;

    REAL chempot2;
    REAL meanval_sq;
	
	REAL gamma_s;
	REAL gamma_reg_s;

	REAL particle_weight;
	REAL ex_strength;
	REAL in_strength;

	REAL particle_length_sq;
	REAL curv_hard;

    EnergyComputer(REAL *data, const int *dsz,  double *cellsize, SphereInterpolator *sp, ParticleGrid<Particle> *pcon, REAL *spimg, int spmult) : EnergyComputerBase(data,dsz,cellsize,sp,pcon,spimg,spmult)
	{}


	void setParameters(REAL pwei,REAL pwid,REAL chempot_connection, REAL length,REAL curv_hardthres, REAL inex_balance, REAL chempot2, REAL meanval_sq)
	{


        this->chempot2 = chempot2;
        this->meanval_sq = meanval_sq;
        
		eigencon_energy = chempot_connection;	
		eigen_energy = 0;	
		particle_weight = pwei;
        
        REAL bal = 1/(1+exp(-inex_balance));
		ex_strength = 2*bal;                         // cleanup (todo)
		in_strength = 2*(1-bal)/length/length;         // cleanup (todo)
//		in_strength = 0.64/length/length;         // cleanup (todo)
	
		particle_length_sq = length*length;
		curv_hard = curv_hardthres;

		REAL sigma_s = pwid;
		gamma_s = 1/(sigma_s*sigma_s);
        gamma_reg_s =1/(length*length/4);

    	
	}  
    
    

	////////////////////////////////////////////////////////////////////////////
	////// External Energy 
	////////////////////////////////////////////////////////////////////////////



	inline REAL computeExternalEnergy(pVector &R, pVector &N, REAL &cap, REAL &len, Particle *dp)
	{

		REAL m = SpatProb(R);
		if (m == 0)
		{
			return -INFINITY;
		}

		#ifdef TIMING
		tic(&odfeval_time);
		#endif

	    REAL Dn = evaluateODF(R,N,len);

		#ifdef TIMING
		toc(&odfeval_time);
		#endif

                #ifdef TIMING
		tic(&externalenergy_time);
		#endif


		REAL Sn = 0;
		REAL Pn = 0;
        REAL Pplus = 0;
        REAL Pmin = 0;
     
        
//         
//        REAL d = 2.0;
//        
//        pVector Rp = R+N*(len*d);
// 	    pcontainer->computeNeighbors(Rp);		
// 		for (;;)
// 		{
// 			Particle *p =  pcontainer->getNextNeighbor();
// 			if (p == 0) break;
// 			if (dp != p)
// 			{
// 				REAL dot = fabs(N*p->N);
//                REAL bw = mdoty(dot);
//                REAL dpos = (p->R-Rp).norm_square();
//                REAL w = mexp(dpos*gamma_reg_s);	                
//                Pplus += w*bw;
//            }
// 		}			
//    
//        pVector Rm = R-N*(len*d);
// 	    pcontainer->computeNeighbors(Rm);		
// 		for (;;)
// 		{
// 			Particle *p =  pcontainer->getNextNeighbor();
// 			if (p == 0) break;
// 			if (dp != p)
// 			{
// 				REAL dot = fabs(N*p->N);
//                REAL bw = mdoty(dot);
//                REAL dpos = (p->R-Rm).norm_square();
//                REAL w = mexp(dpos*gamma_reg_s);	                
//                Pmin += w*bw;
//            }
// 		}			
//    
//        REAL ss = mexp(4*len*d*len*d*gamma_reg_s);
//        
//    
//        
//        
       
		pcontainer->computeNeighbors(R);		
		for (;;)
		{
			Particle *p =  pcontainer->getNextNeighbor();
			if (p == 0) break;
			if (dp != p)
			{
				REAL dot = fabs(N*p->N);
                REAL bw = mbesseli0(dot);
                REAL dpos = (p->R-R).norm_square();
                REAL w = mexp(dpos*gamma_s);	                
                Sn += w*(bw+chempot2)*p->cap ;
                w = mexp(dpos*gamma_reg_s);
                Pn += w*mdoty(dot);
            }
		}			

		#ifdef TIMING
		toc(&externalenergy_time);
		#endif
        REAL energy = 0;
        energy += (2*(Dn/particle_weight-Sn) - (mbesseli0(1.0)+chempot2)*cap)*cap;
      
 //       energy -= chempot2*(mbesseli0(1.0) + 2*Pn - (ss+Pplus+Pmin));
  //      energy -= 5*chempot2*(2*Pn - Pplus- Pmin + 1 - ss);
            
        
        return energy*ex_strength;

        
        




	}


	////////////////////////////////////////////////////////////////////////////
	////// Internal Energy 
	////////////////////////////////////////////////////////////////////////////

	inline REAL computeInternalEnergy(Particle *dp)
	{


		REAL energy = eigen_energy;
	

		if (dp->pID != -1)
			energy += computeInternalEnergyConnection(dp,+1);
		
		if (dp->mID != -1)
			energy += computeInternalEnergyConnection(dp,-1);


		return energy;


	}

	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1)
	{
	
		Particle *p2 = 0;
		int ep2;
		if (ep1 == 1)
			p2 = pcontainer->ID_2_address[p1->pID];
		else 		
			p2 = pcontainer->ID_2_address[p1->mID];
		if (p2->mID == p1->ID)
			ep2 = -1;
		else if (p2->pID == p1->ID)	
			ep2 = 1;
		else 
			fprintf(stderr,"EnergyComputer_connec: Connections are inconsistent!\n");

		if (p2 == 0)
			fprintf(stderr,"bug2");

	
		return computeInternalEnergyConnection(p1,ep1,p2,ep2);

	}        

	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1, Particle *p2, int ep2)
	{

		#ifdef TIMING
		tic(&internalenergy_time);
		#endif

		if ((p1->N*p2->N)*ep1*ep2 > -curv_hard)
			return -INFINITY;

		pVector R1 = p1->R + (p1->N * (p1->len * ep1));
		pVector R2 = p2->R + (p2->N * (p2->len * ep2));

		if ((R1-R2).norm_square() > particle_length_sq)
			return -INFINITY;		

		pVector R = (p2->R + p1->R)*0.5;
        
        if (SpatProb(R) == 0)
            return -INFINITY;
        
        
		REAL norm1 = (R1-R).norm_square();
		REAL norm2 = (R2-R).norm_square();
		


		REAL energy = (eigencon_energy-norm1-norm2)*in_strength;


		#ifdef TIMING
		toc(&internalenergy_time);
		#endif


		return energy;
	}        









};


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

#include "ParticleGrid.cpp"
#include "RJMCMCBase.cpp"

class RJMCMC : public RJMCMCBase
{
	public:

	REAL T_in ;
	REAL T_ex ;
	REAL dens;
	

	REAL p_birth;
	REAL p_death;	
	REAL p_shift;
	REAL p_shiftopt;
    REAL p_cap;



	REAL sigma_g; 
	REAL sigma_g_dir; 
    
	REAL gamma_g; 
	REAL Z_g; 

	REAL dthres;
	REAL nthres;
	REAL T_prop;
	REAL stopprobability;
    REAL del_prob;

    
    
	REAL len_def;
	REAL len_sig;

	REAL cap_def;
	REAL cap_sig;


   
    
	Track TrackProposal, TrackBackup;


	SimpSamp<EndPoint> simpsamp;




	RJMCMC(REAL *points,int numPoints,  REAL *dimg, const int *dsz, double *voxsz, double cellsz) : RJMCMCBase(points,numPoints,dimg,dsz,voxsz,cellsz)
	{
	}	
	
	
	void setParameters(REAL Temp, REAL numit, REAL conprob, REAL plen, REAL curv_hardthres, REAL chempot_particle)
	{

		itmax = numit;

	
		p_birth = 0.35;
		p_death = 0.15;		
		p_shift = 0.25;
		p_shiftopt = 0.1;
        p_cap = 0.0;

		REAL sum = p_birth+p_death+p_shift+p_shiftopt+conprob;
		p_birth /= sum; p_death /= sum; p_shift /= sum; p_shiftopt /= sum;		


	
		T_in = Temp;
		T_ex = 0.01;
        dens = exp(-chempot_particle/T_in)*1; //width*height*depth;	
//        dens = width*height*depth;	
	
       
		len_def = plen;
		len_sig = 0.0;
		cap_def = 1.0;
		cap_sig = 0.0;
				

		// shift proposal
		
		sigma_g = len_def/8.0;	
		sigma_g_dir = 0.25;	
		gamma_g = 1/(sigma_g*sigma_g*2);
		Z_g = pow(2*PI*sigma_g,3.0/2.0)*(PI*sigma_g/len_def);

		// conn proposal
		dthres = len_def;
		nthres = curv_hardthres;
		T_prop = 0.5;
		dthres *= dthres;
		stopprobability = exp(-1/T_prop);
        del_prob = 0.1;
        
        
       
	}	


        void iterate_onestep()
	{
			REAL randnum = mtrand.frand();
	
			///////////////////////////////////////////////////////////////
			//////// Birth Proposal
			///////////////////////////////////////////////////////////////
			if (randnum < p_birth)
			{
    
				#ifdef TIMING
				tic(&birthproposal_time);
				birthstats.propose();
				#endif

				pVector R;
				enc->drawSpatPosition(&R);

				pVector N; N.rand_sphere();
				REAL cap =  cap_def - cap_sig*mtrand.frand();
				REAL len =  len_def;// + len_sig*(mtrand.frand()-0.5);
				Particle prop;
				prop.R = R;
				prop.N = N;
				prop.cap = cap;
				prop.len = len;


				REAL prob =  dens * p_death /((p_birth)*(pcontainer.pcnt+1));
		
				REAL ex_energy = enc->computeExternalEnergy(R,N,cap,len,0);
				REAL in_energy = enc->computeInternalEnergy(&prop);

				prob *= exp((in_energy/T_in+ex_energy/T_ex)) ;

				if (prob > 1 || mtrand.frand() < prob)
				{
					Particle *p = pcontainer.newParticle(R);
					if (p!=0)
					{
						p->R = R;
			            p->N = N;
						p->cap = cap;						
						p->len = len;						
						#ifdef TIMING
						birthstats.accepted();
						#endif
                        if ((in_energy/T_in+ex_energy/T_ex) >0) 
                            stats.accepted_down();
                        else
                            stats.accepted_up();
                        
					}
				}	
                else
                    stats.rejected();
				#ifdef TIMING
				toc(&birthproposal_time);
				#endif
			}
			///////////////////////////////////////////////////////////////
			//////// Death Proposal
			///////////////////////////////////////////////////////////////
			else if (randnum < p_birth+p_death)
			{
				if (pcontainer.pcnt > 0)
				{
					#ifdef TIMING
					tic(&deathproposal_time);	
					deathstats.propose();
					#endif

					int pnum = rand()%pcontainer.pcnt;
					Particle *dp = &(pcontainer.particles[pnum]);
                    if (dp->pID == -1 && dp->mID == -1)
                    {

                        REAL ex_energy = enc->computeExternalEnergy(dp->R,dp->N,dp->cap,dp->len,dp);
                        REAL in_energy = enc->computeInternalEnergy(dp);

                        REAL prob = pcontainer.pcnt * (p_birth) /(dens*p_death); //*SpatProb(dp->R);
                        prob *= exp(-(in_energy/T_in+ex_energy/T_ex)) ;
                        if (prob > 1 || mtrand.frand() < prob)
                        {
                            pcontainer.remove(pnum);
                            #ifdef TIMING
                            deathstats.accepted();
                            #endif
                            if (-(in_energy/T_in+ex_energy/T_ex) >0) 
                                stats.accepted_down();
                            else
                                stats.accepted_up();
                            
                        }
                    }
                    else
                        stats.rejected();
                    
					#ifdef TIMING
					toc(&deathproposal_time);
					#endif
				}			

			}
			///////////////////////////////////////////////////////////////
			//////// Cap change Proposal
			///////////////////////////////////////////////////////////////
// 			else  if (randnum < p_birth+p_death+p_cap)
// 			{
// 				REAL energy = 0;
// 				if (pcontainer.pcnt > 0)
// 				{
// 
// 					int pnum = rand()%pcontainer.pcnt;
// 					Particle *p =  &(pcontainer.particles[pnum]);							
// 					Particle prop_p = *p;
// 		
// 					prop_p.cap = cap_def - cap_sig*mtrand.frand();
// 		
// 					REAL ex_energy = enc->computeExternalEnergy(prop_p.R,prop_p.N,prop_p.cap,p->len,p)
// 							- enc->computeExternalEnergy(p->R,p->N,p->cap,p->len,p);							 
// 					//REAL in_energy = enc->computeExternalEnergy(prop_p.R,prop_p.N,p->cap,p->len,p)
// 					//		- enc->computeExternalEnergy(p->R,p->N,p->cap,p->len,p);							 
// 					REAL prob = exp(ex_energy/T_ex);
// 					// * SpatProb(p->R) / SpatProb(prop_p.R);
// 					if (mtrand.frand() < prob)
// 					{
// 						p->cap = prop_p.cap;
//                         accepted++;
// 					}
// 				
// 				}
// 
// 			}	

			///////////////////////////////////////////////////////////////
			//////// Shift Proposal
			///////////////////////////////////////////////////////////////
			else  if (randnum < p_birth+p_death+p_shift+p_cap)
			{
				REAL energy = 0;
				if (pcontainer.pcnt > 0)
				{
					#ifdef TIMING
					tic(&shiftproposal_time);
					shiftstats.propose();
					#endif

					int pnum = rand()%pcontainer.pcnt;
					Particle *p =  &(pcontainer.particles[pnum]);							
					Particle prop_p = *p;
		
					prop_p.R.distortn(sigma_g);
					prop_p.N.distortn(sigma_g_dir);
					prop_p.N.normalize();
					
		
					REAL ex_energy = enc->computeExternalEnergy(prop_p.R,prop_p.N,p->cap,p->len,p)
							- enc->computeExternalEnergy(p->R,p->N,p->cap,p->len,p);
					REAL in_energy = enc->computeInternalEnergy(&prop_p) - enc->computeInternalEnergy(p);
							 
					REAL prob = exp(ex_energy/T_ex+in_energy/T_in);
					// * SpatProb(p->R) / SpatProb(prop_p.R);
					if (mtrand.frand() < prob)
					{
                        pVector Rtmp = p->R;
                        pVector Ntmp = p->N;
						p->R = prop_p.R;
						p->N = prop_p.N;
                        if (!pcontainer.tryUpdateGrid(pnum))
                        {
                            p->R = Rtmp;
                            p->N = Ntmp;
                        }
						#ifdef TIMING
						shiftstats.accepted();
						#endif
                        if ((in_energy/T_in+ex_energy/T_ex) >0) 
                            stats.accepted_down();
                        else
                            stats.accepted_up();
					}
                    else
                        stats.rejected();

					#ifdef TIMING
					toc(&shiftproposal_time);
					#endif
				
				}

			}	
			else  if (randnum < p_birth+p_death+p_shift+p_shiftopt+p_cap)
			{
				REAL energy = 0;
				if (pcontainer.pcnt > 0)
				{

					int pnum = rand()%pcontainer.pcnt;
					Particle *p =  &(pcontainer.particles[pnum]);							

					bool no_proposal = false;
					Particle prop_p = *p;
					if (p->pID != -1 && p->mID != -1)
					{
						Particle *plus = pcontainer.ID_2_address[p->pID];
						int ep_plus = (plus->pID == p->ID)? 1 : -1;
						Particle *minus = pcontainer.ID_2_address[p->mID];
						int ep_minus = (minus->pID == p->ID)? 1 : -1;
						prop_p.R = (plus->R + plus->N * (plus->len * ep_plus)  + minus->R + minus->N * (minus->len * ep_minus))*0.5;
						prop_p.N = plus->R - minus->R;
						prop_p.N.normalize();
					}
					else if (p->pID != -1)
					{
						Particle *plus = pcontainer.ID_2_address[p->pID];
						int ep_plus = (plus->pID == p->ID)? 1 : -1;
						prop_p.R = plus->R + plus->N * (plus->len * ep_plus * 2);
						prop_p.N = plus->N;
					}
					else if (p->mID != -1)
					{
						Particle *minus = pcontainer.ID_2_address[p->mID];
						int ep_minus = (minus->pID == p->ID)? 1 : -1;
						prop_p.R = minus->R + minus->N * (minus->len * ep_minus * 2);
						prop_p.N = minus->N;
					}
					else 
						no_proposal = true;
		
                    
                    
                
                    
                    
                    
					if (!no_proposal)
					{		
                        
                        Particle prop_r = prop_p;
                        REAL sm =0.2;
                    	prop_r.R.distortn(sigma_g*sm);
                        prop_r.N.distortn(sm*sigma_g_dir);					
                        prop_r.N.normalize();
				                        
                        
						REAL cos = prop_p.N*p->N;				
                        REAL p_rev = PosOrientGaussDistrib((prop_p.R-p->R).norm_square(),cos,
                                                         sigma_g*sigma_g, sigma_g_dir*sigma_g_dir);
						REAL cos_rp = prop_p.N*prop_r.N;				
                        REAL p_rp = PosOrientGaussDistrib((prop_p.R-prop_r.R).norm_square(),cos_rp,
                                                         sigma_g*sigma_g*sm*sm, sm*sm*sigma_g_dir*sigma_g_dir);

						REAL ex_energy = enc->computeExternalEnergy(prop_p.R,prop_p.N,p->cap,p->len,p)
								- enc->computeExternalEnergy(p->R,p->N,p->cap,p->len,p);
						REAL in_energy = enc->computeInternalEnergy(&prop_p) - enc->computeInternalEnergy(p);
                 
					//	REAL prob =0; // exp(ex_energy/T_ex+in_energy/T_in)*p_shift*p_rev/(p_shiftopt+p_shift*p_rev); 
						REAL prob = exp(ex_energy/T_ex+in_energy/T_in)*1/(1+p_shiftopt*p_rp/(p_shift*p_rev)); 
								 //* SpatProb(p->R) / SpatProb(prop_p.R);
              
						if (mtrand.frand() < prob)
						{
                            pVector Rtmp = p->R;
                            pVector Ntmp = p->N;	                            
							p->R = prop_p.R;
							p->N = prop_p.N;
							if (!pcontainer.tryUpdateGrid(pnum))
                            {
                                p->R = Rtmp;
                                p->N = Ntmp;                                     
                            }                                                       
                            
                            if ((in_energy/T_in+ex_energy/T_ex) >0) 
                                stats.accepted_down();
                            else
                                stats.accepted_up();
                      
						}
                        else
                            stats.rejected();
                        
					}
				}

			}	
			else
			{


				if (pcontainer.pcnt > 0)
				{

					#ifdef TIMING
					tic(&connproposal_time);
					connstats.propose();
					#endif

					int pnum = rand()%pcontainer.pcnt;
					Particle *p = &(pcontainer.particles[pnum]);

					EndPoint P;
					P.p = p;
					P.ep = (mtrand.frand() > 0.5)? 1 : -1;

					RemoveAndSaveTrack(P);
					if (TrackBackup.proposal_probability != 0)
					{
    					makeTrackProposal(P);
                        REAL edif = (TrackProposal.energy-TrackBackup.energy)/T_in;                       
						REAL prob = exp(edif)*(TrackBackup.proposal_probability * pow(del_prob,TrackProposal.length))
                                        /(TrackProposal.proposal_probability * pow(del_prob,TrackBackup.length));
						if (mtrand.frand() < prob)
						{
							implementTrack(TrackProposal);
                            if (!TrackProposal.isequal(TrackBackup))
                            {

                                #ifdef TIMING
                                connstats.accepted();
                                #endif
                                if (edif>0) 
                                    stats.accepted_down();
                                else
                                    stats.accepted_up();
                            }
						}
						else
						{
                            stats.rejected();
							implementTrack(TrackBackup);
						}
					}
					else
						implementTrack(TrackBackup);

					#ifdef TIMING
					toc(&connproposal_time);
					#endif
						

				}				



			}

			
	

	}


	void implementTrack(Track &T)
	{
		for (int k = 1; k < T.length;k++)
		{
			pcontainer.createConnection(T.track[k-1].p,T.track[k-1].ep,T.track[k].p,-T.track[k].ep);
		}
	}



	void RemoveAndSaveTrack(EndPoint P)
	{
		EndPoint Current = P;

		int cnt = 0;
		REAL energy = 0;
		REAL AccumProb = 1.0;
		TrackBackup.track[cnt] = Current;

		EndPoint Next;		
        
        

		for (;;)
		{
			Next.p = 0;
			if (Current.ep == 1)
			{
				if (Current.p->pID != -1)
				{
					Next.p = pcontainer.ID_2_address[Current.p->pID];
					Current.p->pID = -1;
					pcontainer.concnt--;
				}
			}
			else if (Current.ep == -1)
			{
				if (Current.p->mID != -1)
				{
					Next.p = pcontainer.ID_2_address[Current.p->mID];
					Current.p->mID = -1;
					pcontainer.concnt--;
				}
			}
			else
				{ fprintf(stderr,"RJMCMC_randshift: Connection inconsistent 3\n"); break; }

			if (Next.p == 0) // no successor
			{
				Next.ep = 0; // mark as empty successor
				break;
			}
			else
			{
				if (Next.p->pID == Current.p->ID)
				{
					Next.p->pID = -1;					
					Next.ep = 1;
				}
				else if (Next.p->mID == Current.p->ID)
				{
					Next.p->mID = -1;					
					Next.ep = -1;
				}
				else
					{ fprintf(stderr,"RJMCMC_randshift: Connection inconsistent 4\n"); break; }
			}
                      

    		computeEndPointProposalDistribution(Current);
          
            AccumProb *= (simpsamp.probFor(Next));
			
            if (Next.p == 0) // no successor -> break
				break;

			energy += enc->computeInternalEnergyConnection(Current.p,Current.ep,Next.p,Next.ep);			 

			Current = Next;	
			Current.ep *= -1;
			cnt++;
			TrackBackup.track[cnt] = Current;
            
                  
	        if (mtrand.rand() > del_prob)
            {                
                break;
            }            

		}		
		TrackBackup.energy = energy;
		TrackBackup.proposal_probability = AccumProb;
		TrackBackup.length = cnt+1;

	}



	void makeTrackProposal(EndPoint P)
	{
		EndPoint Current = P;
		int cnt = 0;
		REAL energy = 0;
		REAL AccumProb = 1.0;
		TrackProposal.track[cnt++] = Current;
		Current.p->label = 1;

		for (;;)
		{

			// next candidate is already connected
			if ((Current.ep == 1 && Current.p->pID != -1) || (Current.ep == -1 && Current.p->mID != -1))
				break;

			// track too long
			if (cnt > 250)
				break;		

			computeEndPointProposalDistribution(Current);

// 			// no candidates anymore
// 			if (simpsamp.isempty())
// 				break;

			int k = simpsamp.draw();

			// stop tracking proposed
			if (k==0)
				break;

			EndPoint Next = simpsamp.objs[k];
			REAL probability = simpsamp.probFor(k);

			// accumulate energy and proposal distribution
			energy += enc->computeInternalEnergyConnection(Current.p,Current.ep,Next.p,Next.ep);			
			AccumProb *= probability;

			// track to next endpoint
			Current = Next;
			Current.ep *= -1;

			Current.p->label = 1;  // put label to avoid loops
			TrackProposal.track[cnt++] = Current;

           
            
		}

		TrackProposal.energy = energy;
		TrackProposal.proposal_probability = AccumProb;
		TrackProposal.length = cnt;

		// clear labels
		for (int j = 0; j < TrackProposal.length;j++)
		{
			TrackProposal.track[j].p->label = 0;
		}

	}




	void computeEndPointProposalDistribution(EndPoint P)
	{
		Particle *p = P.p;
		int ep = P.ep;

		REAL dist,dot;
		pVector R = p->R + (p->N * ep*p->len);
		pcontainer.computeNeighbors(R);
		simpsamp.clear();

		simpsamp.add(stopprobability,EndPoint(0,0));		

		for (;;)
		{
			Particle *p2 =  pcontainer.getNextNeighbor();
			if (p2 == 0) break;
			if (p!=p2 && p2->label == 0)
			{
				if (p2->mID == -1)
				{							
					dist = (p2->R - p2->N * p2->len - R).norm_square();
					if (dist < dthres)
					{
						dot = p2->N*p->N * ep;
//						if (dot > nthres)
						{
							REAL en = enc->computeInternalEnergyConnection(p,ep,p2,-1);
							simpsamp.add(exp(en/T_prop),EndPoint(p2,-1));
						}
					}
				}
 				if (p2->pID == -1)
				{
					dist = (p2->R + p2->N * p2->len - R).norm_square();
					if (dist < dthres)
					{
						dot = p2->N*p->N * (-ep);
//						if (dot > nthres)
						{
							REAL en = enc->computeInternalEnergyConnection(p,ep,p2,+1);
							simpsamp.add(exp(en/T_prop),EndPoint(p2,+1));
						}
					}
				}
			}
		}
	}


};




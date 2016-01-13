
#include "ParticleGrid.cpp"
#include "RJMCMCBase.cpp"
#include <list>


class RJMCMC : public RJMCMCBase
{
	public:

	REAL Temp;
	REAL expected_nump;
    REAL Vol;
	


	REAL sigma_g; 
	REAL sigma_g_dir; 
    

	REAL dthres;
	REAL nthres;
	REAL T_prop;
    
	REAL len_def;

    REAL maxD;

    REAL bias;
    REAL modelpropwid;
    
    int restrictions;
    REAL directional_propsal_distrib;
            
    list<pVector> activeSpots;
    REAL safty_distance;


	RJMCMC(REAL *points,int numPoints,  REAL *vfmap, REAL *dimg, const int *dsz, double *voxsz, double cellsz,REAL *cellsz2,REAL len, int numcores) : RJMCMCBase(points,numPoints,vfmap, dimg,dsz,voxsz,cellsz,cellsz2,len,numcores)
	{
	}	
	
	void setParameters(ParamStruct &info)
	{

		itmax = info.numit;
        gammaDE = info.gammaDE;


        // distribution of proposal types
		p_birth = info.prop_p_birth;
		p_death = info.prop_p_death;
		p_shift = info.prop_p_shift;
        p_Dmod = info.prop_p_Dmod;
        p_vfmod = info.prop_p_vfmod;
        p_conprob = info.prop_p_conprob;
        
  
		REAL sum = p_birth+p_death+p_shift+p_Dmod+p_conprob+p_vfmod;
		p_birth /= sum; p_death /= sum; p_shift /= sum; p_Dmod /= sum; p_conprob /= sum; p_vfmod /=sum;
	
		Temp = info.Temp;
        
        
        // strength of tracking guide is proportional to temperature
        enc->prior_strength = Temp*info.trackingguide_strength;

        directional_propsal_distrib = info.directional_propsal_distrib;
        
        REAL voxelvol_mask = enc->voxsize_sp_w*enc->voxsize_sp_h*enc->voxsize_sp_d;
        REAL voxelvol_data = enc->voxsize_w*enc->voxsize_h*enc->voxsize_d;
        Vol = double(enc->scnt)*voxelvol_mask;
        
        // expected number of segments of the underlying Poison process
        expected_nump = info.particle_expnump*Vol/( enc->sigma_s* enc->sigma_s* enc->sigma_s)/voxelvol_data; 
       
        fprintf(stderr,"#vox mask: %i  expected_nump: %f\n",enc->scnt,expected_nump);
        fprintf(stderr,"voxvol mask: %.2f  voxvol data: %.2f \n",voxelvol_mask,voxelvol_data);
             
        
        restrictions = (int) info.restrictions;
        
        fprintf(stderr,"Ball: %i   Unique parallel diff: %i \n",restrictions&1,restrictions&2);
        
		len_def = info.particle_len;
				
        //--------- modelchange prop. parameters
        modelpropwid = sqrt(Temp);  // sigma of density ~ sqrt(Temp)
        maxD = 5;                   // maximal D 
        
		//--------- shift proposal parameters
		sigma_g = enc->sigma_s*1;  // spatial shift
		sigma_g_dir = 0.1;         // angular shift

		//--------- connection proposal parameters
		dthres = len_def*sqrt(2);  // hard threshold for distance between connected particles
		dthres *= dthres;
		nthres = info.curv_hardthres;   // hard threshold for angle between connected particles
		T_prop = Temp*1;           // temperature of proposal distribution
        
        //--------- parallel impl. safty distance
        safty_distance = len_def*7; // safty dist. for parallel proposals
        safty_distance *= safty_distance;
        
       
	}	

    void iterate_onestep()
	{
            REAL randnum = mtrand.frand();
        
        
        
			///////////////////////////////////////////////////////////////
			//////// Birth Proposal
			///////////////////////////////////////////////////////////////
			if (randnum < p_birth)
			{

                // get a random position and block neighborhood (for openMP)
				pVector R ;
                R = getRandomPosition();                
                
                
                // generate random proposal parameters
                REAL Di,Da,Dp,w,q;
                pVector N;
                if (directional_propsal_distrib>0)
                    N = enc->drawPriorDir(R);
                else
                    N.rand_sphere();        
                
                
                chooseDap(Di,Da,Dp);                
                w = chooseWeight();                             
                q = chooseWeight_guide();
                        
                // wrap it into a particle
                Particle proposal;
                proposal.R = R;
                proposal.N = N;
                proposal.Di = Di;
                proposal.Da = Da;
                proposal.Dp = Dp;
                proposal.w = w;
                proposal.q = q;
                
                // compute Energy differences
                REAL ex_energy;
				ex_energy = enc->computeExternalEnergy(&proposal,0);                
                ex_energy += enc->smoothnessSegs(&proposal);
                ex_energy += enc->computeInternalEnergy(&proposal) ;
                
                // compute Green's ratio
                int parcnt;
                
                parcnt = pcontainer.getNumParticles();
				REAL prob =  expected_nump * p_death /((p_birth)*(parcnt+1));
                
                
           //     pcontainer.getCell(R,parcnt);
           //     REAL prob =  3* p_death /((p_birth)*(parcnt+1));
                
				prob *= exp(ex_energy/Temp) ;

                if (prob > 1 || mtrand.frand() < prob)
                {                        
                    Particle *p = pcontainer.newParticle(R,N);
                    if (p!=0)
                    {
                        p->R = R;
                        p->N = N;
                        p->len = len_def;		
                        p->Di = Di;
                        p->Da = Da;
                        p->Dp = Dp;
                        p->w = w;
                        p->q = q;

                        stats.accepted(ex_energy/Temp);
                        birthstats.accepted(ex_energy/Temp);
                    }
                }
                else
                {
                    stats.rejected();
                    birthstats.rejected();
                }
                
                // release neighborhood
                bool removed = remove_from_activeSpots(R);
                if (!removed)
                      fprintf(stderr,"warning not released -------------------------------------- Birth\n");
                

            }
			///////////////////////////////////////////////////////////////
			//////// Death Proposal
			///////////////////////////////////////////////////////////////
			else if (randnum < p_birth+p_death)
			{
				if (pcontainer.pcnt > 0)
				{

                    // choose a random particle and block neighborhood
                    Particle *dp;
                    int pnum;
                    getRandomParticle(&dp,&pnum);
                    if (dp==0) return;
                    pVector spotR = dp->R;

                    
                    if (dp->pID == -1 && dp->mID == -1) // detah in only possible if there is no connection
                    {

                        // cmopute Energy differences
                        REAL ex_energy;
                        ex_energy = enc->computeExternalEnergy(dp,dp);
                        ex_energy += enc->smoothnessSegs(dp);
                        ex_energy += enc->computeInternalEnergy(dp);
                        
                        // compute Green's ratio
                        int parcnt;
                        parcnt = pcontainer.getNumParticles();      
                        REAL prob = parcnt * (p_birth) /(expected_nump*p_death);

//                        pcontainer.getCell(dp->R,parcnt);
 //                       REAL prob =  parcnt*(p_birth) /(3*p_death);                        
                        
                        prob *= exp(-(ex_energy/Temp));
                        
                        if (directional_propsal_distrib>0)
                            prob *= enc->getPriorDirProb(dp->R,dp->N);
                                                
                        if (prob > 1 || mtrand.frand() < prob)
                        {
                            pcontainer.remove(pnum);
                            stats.accepted(-ex_energy/Temp);
                            birthstats.accepted(-ex_energy/Temp);
                        }
                    }
                    else
                    {
                         stats.rejected();
                         deathstats.rejected();
                    }
                
                    if (!remove_from_activeSpots(spotR))
                          fprintf(stderr,"warning not releases -------------------------------------- death\n");
                    
				}			

			}
        
        
			///////////////////////////////////////////////////////////////
			//////// Diffusion parameter change Proposal
			///////////////////////////////////////////////////////////////
			else  if (randnum < p_birth+p_death+p_Dmod)
			{
        	
            
				REAL energy = 0;
				if (pcontainer.pcnt > 0)
				{

                    // choose a random particle and block neighborhood
                    Particle *p;
                    int pnum;
                    getRandomParticle(&p,&pnum);
                    if (p==0) return;                    
					Particle prop_p = *p;
                                            
                    
                    REAL sel = mtrand.frand();                   
                    if (sel >= 0.1) // either propose new diffusion parameters
                    {                    
                        distortDap(prop_p.Di,prop_p.Da,prop_p.Dp);               
                        prop_p.w =chooseWeight();
                    }
                    else // or just a new tracking guide weight
                        prop_p.q =chooseWeight_guide();
                    
                    
                    // compute Energy differences                    
					REAL ex_energy = enc->computeExternalEnergy(&prop_p,p) - enc->computeExternalEnergy(p,p);			
                    ex_energy += (enc->smoothnessSegs(&prop_p) - enc->smoothnessSegs(p));
                    REAL in_energy = (enc->computeInternalEnergy(&prop_p) - enc->computeInternalEnergy(p));
                    
                    REAL prob = exp(ex_energy/Temp + in_energy/Temp);
                    
                            
					if (mtrand.frand() < prob)
					{                                                
                        p->Da = prop_p.Da;
                        p->Dp = prop_p.Dp;
                        p->Di = prop_p.Di;
                        p->w = prop_p.w;
                        p->q = prop_p.q;
                                                
                        stats.accepted((in_energy/Temp+ex_energy/Temp));
                        capstats.accepted((in_energy/Temp+ex_energy/Temp));
                        
                        accepted++;
                    }
                    else
                    {
                        stats.rejected();
                        capstats.rejected();
                    }
                    
                    remove_from_activeSpots(p->R);                                                
		
				}

			}	
			///////////////////////////////////////////////////////////////
			//////// VF change Proposal
			///////////////////////////////////////////////////////////////        
			else  if (randnum < p_birth+p_death+p_vfmod+p_Dmod)
            {
                
                // get a random position                                    
				pVector R ;
                R = getRandomPosition();                
                if (enc->SpatProb(R) > 0)
                {

                    // find the corresponding voxel
                    REAL Rx = R.x/((REAL)enc->voxsize_w*enc->sigma_s);
                    REAL Ry = R.y/((REAL)enc->voxsize_h*enc->sigma_s);
                    REAL Rz = R.z/((REAL)enc->voxsize_d*enc->sigma_s);
                    int xint = int(floor(Rx));
                    int yint = int(floor(Ry));
                    int zint = int(floor(Rz));
                    int spatindex = (xint + enc->w_vf*(yint+enc->h_vf*zint));                    

                    
                    // compute energy before the change
                    REAL ex_energy = -enc->computeExternalEnergy(R);
                    REAL in_energy = -enc->computeInternalEnergy(R);
                    REAL sm_energy = enc->smoothnessVF(R);
                                                            
                    // do a random change
                    REAL vfisav = vfmap[spatindex];
                    distortWeight(vfmap[spatindex]);
                                
                    // compute energy after change
                    sm_energy -= enc->smoothnessVF(R);                    
                    ex_energy += enc->computeExternalEnergy(R);
                    in_energy += enc->computeInternalEnergy(R);
                    ex_energy += sm_energy;
                    
                    
                    
                    REAL prob = exp(ex_energy/Temp + in_energy/Temp);
                    if (ex_energy/Temp + in_energy/Temp != 0)
                    {
                        if (mtrand.frand() < prob)
                        {  
                            stats.accepted((ex_energy/Temp + in_energy/Temp));
                            vfstats.accepted((ex_energy/Temp + in_energy/Temp));
                            accepted++;                            
                        }
                        else
                        {
                            stats.rejected();
                            vfstats.rejected();

                            vfmap[spatindex] = vfisav;
                            
                        }
                    }    

                }
                remove_from_activeSpots(R);                                                            
            
           }

			///////////////////////////////////////////////////////////////
			//////// Shift Proposal
			///////////////////////////////////////////////////////////////
			else  if (randnum < p_birth+p_death+p_vfmod+p_shift+p_Dmod)
			{
                REAL energy = 0;
				if (pcontainer.pcnt > 0)
				{

                    Particle *p;
                    int pnum;
                    getRandomParticle(&p,&pnum);
                    if (p==0) return;
                                        
					Particle prop_p = *p;
                    pVector oldR = p->R;
                    
                    if (mtrand.frand() > 0.5) // shift position
                    {
                        prop_p.R.distortn(sigma_g ); 
                    }
                    else                      // change direction
                    {
                         prop_p.N.distortn(sigma_g_dir );
                         prop_p.N.normalize(); 
                    }
		
					REAL ex_energy = enc->computeExternalEnergy(&prop_p,p) - enc->computeExternalEnergy(p,p);
                    ex_energy += enc->smoothnessSegs(&prop_p)-enc->smoothnessSegs(p);
					REAL in_energy = enc->computeInternalEnergy(&prop_p) - enc->computeInternalEnergy(p);
							 
                    
                    
					REAL prob = exp(ex_energy/Temp+in_energy/Temp);
                    if (directional_propsal_distrib>0)
                    {
                       REAL pprior_prop = enc->getPriorDirProb(prop_p.R,prop_p.N)/ enc->getPriorDirProb(p->R,p->N);
                       prob *= pprior_prop;
                    }
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
                        else
                        {
                            if ((in_energy/Temp+ex_energy/Temp) != 0)
                            {
                                stats.accepted((in_energy/Temp+ex_energy/Temp));
                                shiftstats.accepted((in_energy/Temp+ex_energy/Temp));
                            }
                        }
                    }
                    else
                    {
                        stats.rejected();
                        shiftstats.rejected();
                    }
                    
                    remove_from_activeSpots(oldR);                                                
				
				}

			}	
        
			///////////////////////////////////////////////////////////////
			//////// Tracking Proposal
			///////////////////////////////////////////////////////////////
			else 
			{


				if (pcontainer.pcnt > 0)
				{
                    

                    Particle *p;
                    int pnum;
                    getRandomParticle(&p,&pnum);
                    if (p==0) return;

                    
                    
					EndPoint P;
					P.p = p;
					P.ep = (mtrand.frand() > 0.5)? 1 : -1;
                    EndPoint Current(0,0);
                                          

                    if (P.ep == 1)                                    
                    {
                        if (P.p->pID != -1)
                        {
                            Current.p = pcontainer.ID_2_address[P.p->pID];
                            if (Current.p->pID == p->ID)     
                                Current.ep = 1;
                            else
                                Current.ep = -1;
                            pcontainer.destroyConnection(p,+1);
                        }
                    }
                    else   
                    {
                        if (P.p->mID != -1)
                        {
                            Current.p = pcontainer.ID_2_address[P.p->mID];   
                            if (Current.p->pID == p->ID)     
                                Current.ep = 1;
                            else
                                Current.ep = -1;
                            pcontainer.destroyConnection(p,-1);
                        }
                    }
                                    
                    
                    
                  	SimpSamp<EndPoint> simpsamp = computeEndPointProposalDistribution(P);
        			int k = simpsamp.draw();

                    EndPoint New = simpsamp.objs[k];
                    REAL probability_cur = simpsamp.probFor(Current);
                    REAL probability_new = simpsamp.probFor(New);

                    
                    
                    REAL energy = 0;
                    energy += enc->computeInternalEnergyConnection(P.p,P.ep,New.p,New.ep);			
                    if (Current.p != 0)
                        energy += enc->computeInternalEnergyConnection(Current.p,Current.ep,0,0);
                    
                    
                    energy -= enc->computeInternalEnergyConnection(P.p,P.ep,Current.p,Current.ep);			
                    if (New.p != 0)
                        energy -= enc->computeInternalEnergyConnection(New.p,New.ep,0,0);
                    
                    
                    REAL prob = exp(energy/Temp) * probability_cur /probability_new;
                        if (mtrand.frand() < prob)
                        {
                            if (New.p != 0)
                            {
                                pcontainer.createConnection(P.p,P.ep,New.p,New.ep);
                            }

                            if (energy != 0)
                            {
                                stats.accepted(energy/Temp);
                                connstats.accepted(energy/Temp);
                            }
                        }
                        else
                        {
                        
                            if (Current.p != 0)
                            {
                                pcontainer.createConnection(P.p,P.ep,Current.p,Current.ep);
                            }
                            stats.rejected();
                            connstats.rejected();
                        }
                    
                    remove_from_activeSpots(p->R);                                                                        

						

				}				
            }

	}

    
    
    
    
    
    
    
	SimpSamp<EndPoint> computeEndPointProposalDistribution(EndPoint P)
	{
		Particle *p = P.p;
		int ep = P.ep;
        NeighborTracker nbtrack_p;
        NeighborTracker nbtrack_m;
        
		REAL dist,dot;
		pVector R = p->R + (p->N * ep*p->len);
        
		pcontainer.computeNeighbors(R,nbtrack_p,+1);
		pcontainer.computeNeighbors(R,nbtrack_m,-1);
        
        int ncnt_p = pcontainer.getNeighborCount(nbtrack_p);
        int ncnt_m = pcontainer.getNeighborCount(nbtrack_m);
        
        SimpSamp<EndPoint> simpsamp(ncnt_p+ncnt_m+1);
        
		simpsamp.clear();
        REAL nonconprob = exp( enc->computeInternalEnergyConnection(p,ep,0,0)  /T_prop);
		simpsamp.add(nonconprob,EndPoint(0,0));		
        
		for (;;)
		{
			Particle *p2 =  pcontainer.getNextNeighbor(nbtrack_m);
			if (p2 == 0) break;
			if (p!=p2)
			{
				if (p2->mID == -1)
				{							
					dist = (p2->R - p2->N * p2->len - R).norm_square();
					if (dist < dthres)
					{
						dot = p2->N*p->N * ep;
						if (dot > nthres)
						{
							REAL en = enc->computeInternalEnergyConnection(p,ep,p2,-1);
							simpsamp.add(exp(en/T_prop),EndPoint(p2,-1));
						}
					}
				}
                
            }
        }
            
		for (;;)
		{
			Particle *p2 =  pcontainer.getNextNeighbor(nbtrack_p);
			if (p2 == 0) break;
			if (p!=p2)
			{
            
 				if (p2->pID == -1)
				{
					dist = (p2->R + p2->N * p2->len - R).norm_square();
					if (dist < dthres)
					{
						dot = p2->N*p->N * (-ep);
						if (dot > nthres)
						{
							REAL en = enc->computeInternalEnergyConnection(p,ep,p2,+1);
							simpsamp.add(exp(en/T_prop),EndPoint(p2,+1));
						}
					}
				}
			}
		}
        
        
        return simpsamp;
	}

    void getRandomParticle(Particle **p, int *num)
    {
        Particle *dp = 0;
        int pnum = -1;
        bool safe = false;
        
        for(;;)
        {

            safe = true;
            if (pcontainer.getNumParticles() > 0)
            {                    

                #ifdef PARALLEL_OPENMP
                #pragma omp critical(PARTMOD2)
                #endif
                {

                    for(int k = 0; k < 100; k++) // try 100 times to get a particle
                    {
                        pnum = rand()%pcontainer.pcnt;
                        dp = &(pcontainer.particles[pnum]);
                        if (dp->active)
                            break;
                    }
                    if (!dp->active)
                    {
                        dp = 0;
                        pnum = -1;
                    }
                    else
                    {
                        list<pVector>::iterator i; 
                        for(i=activeSpots.begin(); i != activeSpots.end(); ++i)    
                        {
                            safe = safe && ((*i-dp->R).norm_square() > safty_distance);
                        }
                        
                        if (safe)
                            activeSpots.push_back(dp->R);                    
                        else
                        {                    
                            dp = 0;
                            pnum = -1;
                        }
                    }
                }
            }			
            else
            {
                dp = 0;
                pnum = -1;
            }

            if (safe)
                break;
            
        }        

        *p = dp;
        *num = pnum;
        
    }
    
    
    
    pVector getRandomPosition()
    {
        pVector R;
        
        bool safe = false;
        while (!safe)
        {
            #ifdef PARALLEL_OPENMP            
            #pragma omp critical(PARTMOD2)
            #endif
            {
            
                enc->drawSpatPosition(&R);
                safe = true;
                list<pVector>::iterator i; 
                for(i=activeSpots.begin(); i != activeSpots.end(); ++i)                        
                    safe = safe && ((*i-R).norm_square() > safty_distance);
                if (safe)
                    activeSpots.push_back(R);
            }
                
        }
        return R;
    }
    
    bool remove_from_activeSpots(pVector R)
    {
        bool sucess = false;
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (PARTMOD2)
        #endif
        {
             list<pVector>::iterator i; 
             for(i=activeSpots.begin(); i != activeSpots.end(); ++i)                        
                 if (*i == R)         
                 {
                    activeSpots.erase(i);        
                    sucess = true;
                    break;
                 }
        }
        
        return sucess;          
    }
    
    
    void chooseDap(REAL &Di, REAL &Da, REAL &Dp)
    {
        
        REAL D1,D2;
        
        Da = mtrand.frand()*maxD;
        
        if (restrictions&1) // ball
            Dp = Da;
        else
            Dp = mtrand.frand()*maxD;
        
        if (restrictions&2) // unique para
            Di = Da;
        else
            Di = mtrand.frand()*maxD;
  
//        Di = Da+2*Dp;
        //if (Di > maxD)
        //    Di = maxD;
    }



    REAL chooseWeight()
    {
        
        REAL t = 0.1;
        REAL s = 80;
//        REAL t = 0.2;
  //      REAL s = 40;
        
        return (t + log(exp(mtrand.frand()*s*(1-t)) - (1-exp(-s*t)) )/s);
        
    }
    
    REAL chooseWeight_guide()
    {
        
        REAL t = 0.1;
        REAL s = 80;
        
        return (t + log(exp(mtrand.frand()*s*(1-t)) - (1-exp(-s*t)) )/s);
                
        //REAL r = powf(mtrand.frand(),1);
        //return r;        
    }
        
    
    void distortWeight(REAL &w)
    {

        while(true)
        {
            REAL nw = w + mtrand.frandn()*modelpropwid*3;            
            if (nw >= 0 && nw <= 1)
            {
                w = nw;
                break;
            }
        }        
        
    }

   

    void distortDap(REAL &Di,REAL &Da, REAL &Dp)
    {   
        REAL Da_,Dp_,Di_;
        REAL min = 0.0;
        REAL mpc = 1;
        while(true)
        {            
            Di_ = Di + mtrand.frandn()*modelpropwid*maxD*mpc; 
            if (Di_ > min && Di_ < maxD)
                break;
        }
        while(true)
        {            
            Da_ = Da + mtrand.frandn()*modelpropwid*maxD*mpc; 
            if (Da_ > min && Da_ < maxD)
                break;
        }        
        while(true)
        {            
            Dp_ = Dp + mtrand.frandn()*modelpropwid*maxD*mpc; 
            if (Dp_ > min && Dp_ < maxD)
                break;
        }
        
        
       if (restrictions & 2) // unique para
               Di_ = Da_;
       if (restrictions & 1) // ball
               Dp_ = Da_;
               
        
//        Di_ = Da_+2*Dp_;
     //   if (Di_ > maxD)
     //       Di_ = maxD;
        
        Di = Di_;
        Da = Da_;
        Dp = Dp_;                   

        
    }
        
 
           
 
    
    
    

};


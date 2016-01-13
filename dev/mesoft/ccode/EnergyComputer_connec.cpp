
#include "ParticleGrid.cpp"
#include "EnergyComputerBase.cpp"


class EnergyComputer : public EnergyComputerBase
{

public:


    REAL kappa;
	REAL eigencon_energy;
    REAL bound_attract;
    
	REAL eigen_energy;

	REAL ex_strength;
	REAL in_strength;

	REAL particle_length_sq;
	REAL curv_hard;
    
    REAL vf_smooth;
    REAL fiber_smooth;
    REAL penalty_sw;
    REAL segment_smooth;
    REAL segment_smooth_range;
    
    int restrictions;
    REAL trace_equality_constraint;

    EnergyComputer(REAL *data, const int *dsz,  double *cellsize, SphereInterpolator *sp, ParticleGrid<Particle> *pcon, REAL *vfmap, REAL *spimg, REAL *s2andidx,int spmult) : EnergyComputerBase(data,dsz,cellsize,sp,pcon,vfmap, spimg,s2andidx,spmult)
	{}

    
    void setParameters(ParamStruct &info)
    {
		
		particle_weight = info.particle_weight;
        
		ex_strength = info.num_qspace_samples;
        particle_length_sq = (info.particle_len*info.particle_len);
		in_strength = info.conpot_lambda/particle_length_sq;
		eigencon_energy = info.chempot_connection*particle_length_sq	;
        kappa = info.kappa_connection;
        bound_attract = info.bound_attract;
		
		curv_hard = info.curv_hardthres;

		sigma_s = info.particle_width;

        w_vf = int((REAL)w/sigma_s);
        h_vf = int((REAL)h/sigma_s);
        d_vf = int((REAL)d/sigma_s);
        
              
        vf_smooth =  info.voxelSmooth_vf;
        fiber_smooth = info.fiberSmooth;
        

        eigen_energy = info.chempot_particle;

        
        penalty_sw = info.penalty_SW;
        segment_smooth = info.voxelSmooth_Ds;
        segment_smooth_range = info.voxelSmooth_Ds_range;
               
        restrictions = (int) info.restrictions;
        trace_equality_constraint = info.trace_equality_constraint;
        
        this->meanval_sq = 2*info.meanval_sq -info.meanval_sq*info.meanval_sq;
    
        GM_mask = info.GM_mask > 0;
        
        
        if (meanval_sq == 0) 
             penalty_sw = 0;
        
	}  
    
    

	////////////////////////////////////////////////////////////////////////////
	////// External Energy 
	////////////////////////////////////////////////////////////////////////////



    //--------- computes data likelihood at a specific location
	inline REAL computeExternalEnergy(pVector &R)
	{

        
        REAL vfi = 0;      
        getVF_int(R,vfi);        
        
        
        REAL ss = 1/sigma_s;

		REAL MM = 0;
        REAL SM = 0;
                    
        REAL Ms = getMS(R);
    
        REAL res_dummy,res_sw,total_proj_sw;
        total_proj_sw = 0;
        int cnt1;
        Particle **P1 = pcontainer->getCell(R,cnt1);
        for (int k = 0; k < cnt1; k++)
        {
            Particle *p1 = P1[k];
            SM += getModelSignalCorrelation(p1,vfi);
            for (int j = k; j < cnt1; j++)
            {
                Particle *p2 = P1[j];
                MM += ((k!=j)?2:1)*
                        mminteract(p2->q,p2->w,p2->Di,p2->Da,p2->Dp,vfi,0,
                                   p1->q,p1->w,p1->Di,p1->Da,p1->Dp,vfi,0,  p2->N*p1->N, res_dummy,res_sw);
            }
            total_proj_sw += res_sw;
		}			

        
        REAL alpha_sw = Ms/particle_weight - total_proj_sw;                           

        if (meanval_sq > 0) // only if sw is modelled imlpicilty
        {
            if (alpha_sw < 0)
            {                        
                return -INFINITY;
            }
        }
        
        REAL M0 = 1/particle_weight;
        REAL energy =   2*SM*M0 - MM - penalty_sw*alpha_sw;
         
    
        return energy*ex_strength;
        
	}
    
 
    //--------- computes data likelihood for specific particle    
	inline REAL computeExternalEnergy(Particle *tp ,Particle *dp)
	{

		REAL m = SpatProb(tp->R);
		if (m == 0)
		{
			return -INFINITY;
		}

        REAL vfi;
        getVF_int(tp->R,vfi);        
        
               
	    REAL Dn = getModelSignalCorrelation(tp,vfi);
        

		REAL Mn = 0;
        
        REAL Ms = getMS(tp->R);
        
        
        REAL res_cur_sw_proj = 0;
        REAL res_oth_sw_proj = 0;
        
        REAL res_tmp;
        
        int cnt1;
        Particle **P1 = pcontainer->getCell(tp->R,cnt1);
                
        for (int k = 0; k < cnt1; k++)
        {
            Particle *p = P1[k];
			if (dp == 0 || dp->ID != p->ID)
			{                              
                REAL dot = fabs(tp->N*p->N);

                Mn +=  mminteract(tp->q,tp->w,tp->Di,tp->Da,tp->Dp,vfi,0,
                                   p->q, p->w, p->Di, p->Da, p->Dp,vfi,0,
                                              dot,res_cur_sw_proj,res_tmp);
                res_oth_sw_proj += res_tmp;
            }          

        }
                
        REAL selfinteract =  mminteract(tp->q,tp->w,tp->Di,tp->Da,tp->Dp,vfi,0,
                                        tp->q,tp->w,tp->Di,tp->Da,tp->Dp,vfi,0, 1.0,res_cur_sw_proj,res_cur_sw_proj) ;
                
        REAL alpha_sw = Ms/particle_weight - (res_oth_sw_proj+res_cur_sw_proj);                
        
        if (meanval_sq > 0) // only if sw is modelled imlpicilty
        {        
            if (alpha_sw < 0)
            {
               return -INFINITY;
            }
        }
        
        REAL M1 = 1/particle_weight;
         
        REAL energy = 0;
              
        energy +=  (2*(M1*Dn-Mn) - selfinteract + penalty_sw*res_cur_sw_proj)*ex_strength;
          
        energy -= eigen_energy *ex_strength;
        
        REAL dd = tp->Di- (tp->Da+2*tp->Dp);
        energy -= trace_equality_constraint*dd*dd;
        
        
        return energy;
  
	}
    
    


	////////////////////////////////////////////////////////////////////////////
	////// Internal Energy 
	////////////////////////////////////////////////////////////////////////////

    
    //--------------------------------------------------------------------------
    //-- connection potentials
    //--------------------------------------------------------------------------
	inline REAL computeInternalEnergy(Particle *dp)
	{
        REAL energy = 0;
		if (dp->pID != -1)
        {
            Particle *dp2 = pcontainer->ID_2_address[dp->pID];
            if (dp2->mID == dp->ID)
                energy += computeInternalEnergyConnection(dp,+1,dp2,-1);    
            else if (dp2->pID == dp->ID)
                energy += computeInternalEnergyConnection(dp,+1,dp2,+1);    
            else
            {
                fprintf(stderr,"EnergyComputer_connec: Connections are inconsistent!\n");
            }
        }
        else
            energy += computeInternalEnergyConnection(dp,+1,0,0);    
        
		if (dp->mID != -1)
        {
            Particle *dp2 = pcontainer->ID_2_address[dp->mID];
            if (dp2->mID == dp->ID)
                energy += computeInternalEnergyConnection(dp,-1,dp2,-1);    
            else if (dp2->pID == dp->ID)
                energy += computeInternalEnergyConnection(dp,-1,dp2,+1);    
            else
            {
                fprintf(stderr,"EnergyComputer_connec: Connections are inconsistent!\n");
            }
        }
        else
            energy += computeInternalEnergyConnection(dp,-1,0,0);    
        
        return energy;
		
	}
    
    
	inline REAL computeInternalEnergyConnection(Particle *p1,int ep1, Particle *p2, int ep2)
	{


        if (bound_attract>0)
            if (p2 == 0) 
            {
                // if endpoint is inside graymatter we have some negative cost, but
                // only if other port is occuppied
                if (ep1 == +1 && p1->mID == -1)
                    return 0;
                if (ep1 == -1 && p1->pID == -1)
                    return 0;
                pVector R = p1->R + (p1->N * (p1->len * ep1));
                if (isGrayMatter(R) > 0)
                    return eigencon_energy*in_strength*bound_attract; 
                else
                    return 0;
            }
        
        if (p2 == 0)
            return 0;
        
                
        REAL dot = -1*(p1->N*p2->N)*ep1*ep2;
                
		if (dot < curv_hard)
			return -INFINITY;

		pVector R1 = p1->R + (p1->N * (p1->len * ep1));
		pVector R2 = p2->R + (p2->N * (p2->len * ep2));

                
		if ((R1-R2).norm_square() > particle_length_sq)
			return -INFINITY;		

		pVector R = (p2->R + p1->R)*0.5;
        pVector dR = (p2->R - p1->R);
        dR.normalize();
        
        if (SpatProb(R) == 0)
            return -INFINITY;
        
                        
		REAL norm1 = (R1-R).norm_square();
		REAL norm2 = (R2-R).norm_square();
        
		REAL norm1_ang = (R1-R)*dR; norm1_ang *= norm1_ang;
		REAL norm2_ang = (R2-R)*dR; norm2_ang *= norm2_ang;
        
        REAL bendEnergy = (kappa * (norm1_ang + norm2_ang) + (1-kappa)*(norm1-norm1_ang + norm2-norm2_ang))*2;
        
	
        REAL energy;
        if (fiber_smooth>0)
        {
            REAL vfi1=0;
            REAL vfi2=0;
            getVF_int(p1->R,vfi1);
            getVF_int(p2->R,vfi2);
            REAL dfi = vfi1-vfi2; dfi *= dfi;
            REAL dDi = p1->Di-p2->Di; dDi *= dDi;
            REAL dDa = p1->Da-p2->Da; dDa *= dDa;
            REAL dDp = p1->Dp-p2->Dp; dDp *= dDp;
            
            
            REAL f_smooth = 0;
            if (restrictions&1 && restrictions&2)
                f_smooth = dDi;
            else if (restrictions&1)
                f_smooth = dDi+dDa;
            else if (restrictions&2)
                f_smooth = dDi+dDp;
            else 
                f_smooth = dDi+dDa+dDp;
                
            
    		energy = (eigencon_energy-bendEnergy - fiber_smooth*f_smooth )*in_strength; 
        }
        else
            energy = (eigencon_energy-bendEnergy)*in_strength; 


		return energy;
	}        
    
	inline REAL computeInternalEnergy(pVector &R)
    {
        REAL intEnergy = 0;
        
        int cnt1;
        Particle **P1 = pcontainer->getCell(R,cnt1);
        for (int k = 0; k < cnt1; k++)
        {
            Particle *p1 = P1[k];            
            intEnergy += computeInternalEnergy(p1);
        }
        return intEnergy;
    }    
    


    
    //--------------------------------------------------------------------------
    //-- smoothness prior for internal volume fraction
    //--------------------------------------------------------------------------
    
    inline REAL smoothnessVF(pVector &R)
    {
  
        if (vf_smooth > 0)
        {
            REAL dx = (sigma_s*voxsize_w); 
            REAL dy = (sigma_s*voxsize_h); 
            REAL dz = (sigma_s*voxsize_d); 


            REAL devi = 0;
            REAL devw = 0;

            REAL vfi;
            getVF_int(R, vfi);

            pVector Rs;
            REAL vfiN;

            Rs = R; Rs.x += dx;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }
            Rs = R; Rs.x -= dx;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }
            Rs = R; Rs.y += dy;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }
            Rs = R; Rs.y -= dy;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }
            Rs = R; Rs.z += dz;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }
            Rs = R; Rs.z -= dz;
            if (SpatProb(Rs)>0)
            {
                getVF_int(Rs,vfiN);        
                devi += (vfiN-vfi)*(vfiN-vfi);
            }

            return vf_smooth*devi;
        }
        else
            return 0;
    }
    
    //--------------------------------------------------------------------------
    //-- smoothness prior for diffusion parameters
    //--------------------------------------------------------------------------
    
    inline REAL smoothnessSegs(Particle *p)
    {
            REAL energy = 0;
            if (segment_smooth > 0)
            {
                REAL dx = (sigma_s*voxsize_w); 
                REAL dy = (sigma_s*voxsize_h); 
                REAL dz = (sigma_s*voxsize_d); 
 
                 pVector R = p->R;
 
                 pVector Rs = R; 
                 energy += Particle_smoothnessInteraction(Rs,p);
                 Rs = R; Rs.x += dx;
                 energy += Particle_smoothnessInteraction(Rs,p);
                Rs = R; Rs.x -= dx;
                energy += Particle_smoothnessInteraction(Rs,p);
                Rs = R; Rs.y += dy;
                energy += Particle_smoothnessInteraction(Rs,p);
                Rs = R; Rs.y -= dy;
                energy += Particle_smoothnessInteraction(Rs,p);
                Rs = R; Rs.z += dz;
                energy += Particle_smoothnessInteraction(Rs,p);
                Rs = R; Rs.z -= dz;
                energy += Particle_smoothnessInteraction(Rs,p);
            }
            return -energy;
    }
    
    
  	inline REAL Particle_smoothnessInteraction(pVector &R,Particle *p)
    {
        REAL intEnergy = 0;
        
        int cnt1;
        if (SpatProb(R)>0)
        {
             Particle **P1 = pcontainer->getCell(R,cnt1);
             for (int k = 0; k < cnt1; k++)
             {
                Particle *p1 = P1[k];      
                if (p1->ID != p->ID)
                {
                    REAL dot = fabs(p1->N*p->N);
                    REAL wu = powf(dot,segment_smooth_range)*segment_smooth;
                    intEnergy += ( (p->Di-p1->Di)*(p->Di-p1->Di)+ (p->Da-p1->Da)*(p->Da-p1->Da) + (p->Dp-p1->Dp)*(p->Dp-p1->Dp) )*wu;                                                                        
                }
             }
        }
        return intEnergy;
    }    

    
    
   








};


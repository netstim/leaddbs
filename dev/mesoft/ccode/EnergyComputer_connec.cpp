
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
		
        
        
   
        numbshells = int(info.numbvals);
        lmax =  int(info.lmax/2+1);
        nip = sinterp->nverts;
	//	fprintf(stderr,"Data directions (sizeLUT * dirs) : %i * %i\n",nalpha,nip);
                  
        
        
        
		particle_weight = info.particle_weight;
        
		ex_strength = info.data_lambda;///num_qspace_samples;
        particle_length_sq = (info.particle_len*info.particle_len);
		in_strength = info.conpot_lambda;
        kappa = info.kappa_connection;
        bound_attract = info.bound_attract;
		
		curv_hard = info.curv_hardthres;

		sigma_s = info.particle_width;

        w_vf = int((REAL)w/sigma_s);
        h_vf = int((REAL)h/sigma_s);
        d_vf = int((REAL)d/sigma_s);
        
              
        vf_smooth =  info.voxelSmooth_vf;
        fiber_smooth = info.fiberSmooth;
        

        eigen_energy = info.chempot_particle;                    // 0.1
		eigencon_energy = info.chempot_connection*eigen_energy;  // 0.1

        
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
	inline REAL computeExternalEnergy(Particle *tp, Particle *dp)
	{
        
        REAL m = SpatProb(tp->R);
		if (m == 0)
		{
			return -INFINITY;
		}
        
        
        
        
        REAL *alphaSM_self = (REAL*) malloc(sizeof(REAL)*lmax*numbshells); 
        REAL *alphaMM_self = (REAL*) malloc(sizeof(REAL)*lmax); 
        for (int k = 0; k < lmax; k++)
            alphaMM_self[k] = 0;
        for (int k = 0; k < lmax*numbshells; k++)
            alphaSM_self[k] = 0;

        
        REAL *alphaSM = (REAL*) malloc(sizeof(REAL)*lmax*numbshells); 
        REAL *alphaMM = (REAL*) malloc(sizeof(REAL)*lmax); 
        for (int k = 0; k < lmax; k++)
            alphaMM[k] = 0;
        for (int k = 0; k < lmax*numbshells; k++)
            alphaSM[k] = 0;

        
        
        
        
        int cnt1;
        getModelSignalCorrelation_modelfree(tp, alphaSM_self)   ;         
        mminteract_modelfree(1, tp->w*tp->w,alphaMM_self );
        Particle **P1 = pcontainer->getCell(tp->R,cnt1);
        if (cnt1 > 0)
        {
            for (int k = 0; k < cnt1; k++)
            {
                Particle *p1 = P1[k];            
                if (p1 != dp)
                {
                    REAL fac = 2 * tp->w*p1->w;
                    mminteract_modelfree(tp->N*p1->N, fac,alphaMM_self );
                }
            }
            
        }

      
        int at_least_one_other = 0;
        for (int k = 0; k < cnt1; k++)
        {
            Particle *p1 = P1[k];            
            if (p1 != dp)
            {
                at_least_one_other = 1;
                getModelSignalCorrelation_modelfree(p1, alphaSM)   ;         
                for (int j = k; j < cnt1; j++)
                {
                    Particle *p2 = P1[j];
                    if (p2 != dp)
                    {
                        REAL fac = ((k!=j)?2:1) * p1->w*p2->w;
                        mminteract_modelfree(p2->N*p1->N, fac,alphaMM );
                    }
                }            
            }
        }			
            
        
        int show = 0;
        if (mtrand.rand() > 0.99)
        {
            show = 1;
        show = 0;
            //fprintf(stderr,"show %f \n",INFINITY);
        }
        
        
            
        int cigar = 0;
        for (int j = 0; j < numbshells;j++)
        {
            int k = 1;

            if (alphaSM[k+j*lmax]> 0 || alphaSM_self[k+j*lmax] > 0)
            {
                cigar = 1;
                break;
            }
            
            
        }
         dbgflag = 0;
        if (cigar == 1)
        {
            
            free(alphaSM);
            free(alphaMM);
            free(alphaSM_self);
            free(alphaMM_self);
            dbgflag = 0;
            
            return -INFINITY;
        }
        
        
        
        REAL m0 = 0;
        REAL m0_with =0;      
            
        REAL energy = 0;
       
        REAL myeps = 0.00000;
        
        int showss =mtrand.frand() >0.9999;
        
        for (int k = 0; k < lmax; k++)
        {
            REAL SM_with = 0;
            REAL SM = 0;
            for (int j = 0; j < numbshells; j++)
            {
                REAL sm_with = (alphaSM[k+j*lmax]+alphaSM_self[k+j*lmax]); 
                REAL sm = alphaSM[k+j*lmax];
                
                SM_with += sm_with*sm_with;
                SM += sm*sm;
            }
            const REAL myeps = 0;
            REAL pp = SM_with / fabs(alphaMM[k]+alphaMM_self[k]+myeps);
            REAL ww = SM / fabs(alphaMM[k]+myeps);
                  
            
            
            
            
//             REAL thres = 0.0001;
//             if (pp > thres)
//                 pp = thres;
//             if (ww > thres)
//                 ww = thres;
//             


            REAL facy = 1;// (lmax+1);
            
//             if (showss) 
//             {
//                 fprintf(stderr,"%i) %f %f %f %f %f %f %f\n",k,energy,facy,pp,ww,SM_with,SM,alphaMM[k]+alphaMM_self[k]);
//             }
//             
            
            if (at_least_one_other == 1)
                energy += (pp-ww) * facy;
            else
                energy += (pp ) * facy;
            
            

            if (show==1)
            {
              if (k == 1 && SM_with / (alphaMM[k]+alphaMM_self[k]+myeps)  > 1)
                {
                
//                fprintf(stderr,"%i %i(%f %f),\n ",k, at_least_one_other,  SM_with / (alphaMM[k]+alphaMM_self[k]) ,   SM / (alphaMM[k])  );
                  fprintf(stderr,"%i %i (%f %f %f %f),\n ",k, at_least_one_other, SM_with / (alphaMM[k]+alphaMM_self[k]+myeps) ,SM_with, alphaMM[k], alphaMM_self[k]);                
                
                  for (int k = 0; k < cnt1; k++)
                    {
                        Particle *p1 = P1[k];            
                        {
                            for (int j = k; j < cnt1; j++)
                            {
                                Particle *p2 = P1[j];
                                {
                                    fprintf(stderr,"%f ",p2->N*p1->N);
                                }
                            }            
                        }
                    }			
                    fprintf(stderr,"\n\n"  );
                
                
                }
            }

        }
        
        
//         if (show==1)
//         {
//            fprintf(stderr,"\n\n");
//         }
//         if (show==1)
//             fprintf(stderr,"%f,\n\n ",energy);
//         
        energy = energy*ex_strength;
        

        energy -= eigen_energy;
            
            
            
        free(alphaSM);
        free(alphaMM);
        free(alphaSM_self);
        free(alphaMM_self);
                    
        
        return energy;
        
	}
    
//     
// 	inline REAL computeExternalEnergy(pVector &R)
// 	{
//         
//         
//         REAL energy = 0;
// 
//         int cnt1;
//         Particle **P1 = pcontainer->getCell(R,cnt1);
//         if (cnt1 > 0)
//         {
//             REAL *alphaSM = (REAL*) malloc(sizeof(REAL)*nalpha); 
//             REAL *alphaMM = (REAL*) malloc(sizeof(REAL)*lmax); 
//             for (int k = 0; k < lmax; k++)
//                 alphaMM[k] = 0;
//             for (int k = 0; k < nalpha; k++)
//                 alphaSM[k] = 0;
//             
//             int numsegs = 0;
//             
//             for (int k = 0; k < cnt1; k++)
//             {
//                 Particle *p1 = P1[k];            
//                 if (p1->active)
//                 {
//                     numsegs++;
//                     getModelSignalCorrelation_modelfree(p1, alphaSM)   ;         
//                     for (int j = k; j < cnt1; j++)
//                     {
//                         Particle *p2 = P1[j];
//                         if (p2->active)
//                         {
//                             REAL fac = ((k!=j)?2:1) * p1->w*p2->w;
//                             mminteract_modelfree(p2->N*p1->N, fac,alphaMM );
//                         }
//                     }            
//                 }
//             }			
//             
//             int show = 0;
//             if (mtrand.rand() > 0.99)
//                 show = 1;
//             show = 0;
//             
//             if (numsegs>0)
//             {
//                 for (int k = 0; k < lmax; k++)
//                 {
//                     REAL SM = 0;
//                     for (int j = 0; j < numbshells; j++)
//                     {
//                         SM += alphaSM[k+j*lmax]*alphaSM[k+j*lmax];
//                     }
//                     energy += SM / alphaMM[k] *(lmax+1);
//                      if (show==1)
//                      {
//                        fprintf(stderr,"%i(%f %f %f),\n ",k,sqrt(SM),alphaMM[k],SM / alphaMM[k]);
//                      }
// 
//                  }
//                 if (show==1)
//                 {
//                  fprintf(stderr,"\n\n");
//                 }
//             }
//             if (show==1)
//                 fprintf(stderr,"%f,\n\n ",energy);
//             
//             energy -= eigen_energy*numsegs;
//             
//             
//             
//             free(alphaSM);
//             free(alphaMM);
//             
//         }
//         
//         
//         return energy*ex_strength;
//         
// 	}
//     
    
    
//     
//     
// 
//     //--------- computes data likelihood for specific particle    
// 	inline REAL computeExternalEnergy(Particle *tp ,Particle *dp)
// 	{
// 
// 		REAL m = SpatProb(tp->R);
// 		if (m == 0)
// 		{
// 			return -INFINITY;
// 		}
// 
// 
//         REAL energy = computeExternalEnergy(tp->R);
//         tp->active = false;
//         energy -= computeExternalEnergy(tp->R);
//         tp->active = true;
//                                               
//         
//         return energy;
//   
// 	}
//     
//     
    
    


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
                    return eigencon_energy*bound_attract; 
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
        
	 //  REAL bendEnergy = (R1-R2).norm_square();
     
        REAL energy;
//         if (fiber_smooth>0)
//         {
//             REAL vfi1=0;
//             REAL vfi2=0;
//             getVF_int(p1->R,vfi1);
//             getVF_int(p2->R,vfi2);
//             REAL dfi = vfi1-vfi2; dfi *= dfi;
//             REAL dDi = p1->Di-p2->Di; dDi *= dDi;
//             REAL dDa = p1->Da-p2->Da; dDa *= dDa;
//             REAL dDp = p1->Dp-p2->Dp; dDp *= dDp;
//             
//             
//             REAL f_smooth = 0;
//             if (restrictions&1 && restrictions&2)
//                 f_smooth = dDi;
//             else if (restrictions&1)
//                 f_smooth = dDi+dDa;
//             else if (restrictions&2)
//                 f_smooth = dDi+dDp;
//             else 
//                 f_smooth = dDi+dDa+dDp;
//                 
//             
//     		energy = (eigencon_energy-bendEnergy - fiber_smooth*f_smooth )*in_strength; 
//         }
//         else
            energy = eigencon_energy-bendEnergy/particle_length_sq*in_strength; 


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


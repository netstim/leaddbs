

#ifndef _PARTICLEGRID
#define _PARTICLEGRID

#include "auxilary_classes.cpp"
#include <list>


struct NeighborTracker  // to run over the neighbors
{
    int cellidx[8];
    int cellidx_c[8];

    int cellcnt;
    int pcnt;
    
    int gnum;

};

///////////////////////////////////////////////////////////////////////////////////
// we have three gridarrays: two for finding neighbors to make connections (grid 0,1)
// and one grid at the finest granularity for finding particles that share a voxel
//
// grid 0(1) contains only those particles that have unconnected port 0(1)
///////////////////

#define NUMGRIDS 3


template <class T>
class ParticleGrid
{

	//////////////// Container
public:
	T *particles; // particles in linear array
	int pcnt; // actual number of particles
	int concnt; // number of connections
	int celloverflows;

	T **ID_2_address;

    list<int> freeslots; // keeps track of empty slots in the particle list
    
    bool needs_reallocation; // if there was not enough memory to create a particle
    
    REAL len;
    
public:
    
    

	int  capacity; // maximal number of particles
	int increase_step;

	/////////////////// Grid
    
    class GridArray
    {
    
    public:
        
        T **grid;   // the grid

        // grid size	
        int  nx;
        int  ny;
        int  nz;

        // scaling factor for grid
        REAL mulx;
        REAL muly;
        REAL mulz;


        int  csize; // particle capacity of single cell in grid
        int *occnt; // occupation count of grid cells
        int  gridsize; // total number of cells
        
        REAL cellsize_x; // size of cells in mm
        REAL cellsize_y; // size of cells in mm
        REAL cellsize_z; // size of cells in mm
        
        
        
        GridArray()
        {
            nx = 0; ny = 0; nz = 0; csize = 0;
            gridsize = 0;
            grid = (T**) 0;
            occnt = (int*) 0;
                        
        }
        
        
        int init(int _nx,int _ny,int _nz,int cellcapacity, REAL cellsize_x, REAL cellsize_y, REAL cellsize_z)
        {

            ////// involvin the grid
            nx = _nx; ny = _ny; nz = _nz; csize = cellcapacity;
            gridsize = nx*ny*nz;
            
    
            grid = (T**) malloc(sizeof(T*)*nx*ny*nz*csize);
            occnt = (int*) malloc(sizeof(int)*nx*ny*nz);

            if (grid == 0 || occnt == 0)
            {
                fprintf(stderr,"error: Out of Memory\n");
                return -1;
            }

            for (int i = 0;i < gridsize;i++)
                occnt[i] = 0;

            
            
            
            
            
            this->cellsize_x = cellsize_x;
            this->cellsize_y = cellsize_y;
            this->cellsize_z = cellsize_z;
            mulx = 1/cellsize_x;
            muly = 1/cellsize_y;
            mulz = 1/cellsize_z;  
            
            return 0;
        }

        void freeMem()
        {
            if (grid != 0)
                free(grid);
            if (occnt != 0)
                free(occnt);        
        }
        
        int getIndex(pVector R)
        {
            int xint = int(R.x*mulx);
            if (xint < 0) { fprintf(stderr,"error: out of grid\n");
                return -1;}
            if (xint >= nx)  {  fprintf(stderr,"error: out of grid\n"); 
                return -1;}
            int yint = int(R.y*muly);
            if (yint < 0)  { fprintf(stderr,"error: out of grid\n"); 
                return -1;}
            if (yint >= ny)  { fprintf(stderr,"error: out of grid\n"); 
                return -1;}
            int zint = int(R.z*mulz);
            if (zint < 0) { fprintf(stderr,"error: out of grid\n"); 
                return -1;}
            if (zint >= nz)  { fprintf(stderr,"error: out of grid\n"); 
                return -1;}

            return xint + nx*(yint + ny*zint);
        }
        
    } g1[NUMGRIDS];


public:


	ParticleGrid()
	{

		//// involving the container
		capacity = 0;
		particles = 0;
		ID_2_address = 0;
		pcnt = 0;		
		concnt = 0;
		celloverflows = 0;		
        
        needs_reallocation = false;
		increase_step = 100000;

	}
	

    //-------------- does all the memory allocs
	int allocate(int _capacity,
			  int _nx, int _ny, int _nz, REAL cellsize, int cellcapacity,
			  int _nx2, int _ny2, int _nz2, REAL *cellsize2, int cellcapacity2, REAL len) 
	{

        this->len = len;

		//// involving the container
		capacity = _capacity;
		particles = (T*) mxMalloc(sizeof(T)*capacity);
		ID_2_address = (T**) mxMalloc(sizeof(T*)*capacity);

		if (particles == 0 || ID_2_address == 0)
		{
			fprintf(stderr,"error: Out of Memory\n");
			capacity = 0;
			return -1;
		}
		else
		{
			fprintf(stderr,"Allocated Memory for %i particles \n",capacity);		
		}

		pcnt = 0;	
        int i;
		for (i = 0;i < capacity;i++)
		{
			ID_2_address[i] = &(particles[i]);   // initliaze pointer in LUT
			particles[i].ID = i;		     // initialize unique IDs
		}

		////// involvin the grid
        if (g1[0].init(_nx,_ny,_nz,cellcapacity,cellsize,cellsize,cellsize) == -1)
            return -1;
        if (g1[1].init(_nx,_ny,_nz,cellcapacity,cellsize,cellsize,cellsize) == -1)
            return -1;
        if (g1[2].init(_nx2,_ny2,_nz2,cellcapacity2,cellsize2[0],cellsize2[1],cellsize2[2]) == -1)
            return -1;

        return 1;
	}


    //-------------- if we need more memory we increase the particle list size
	int reallocate()
	{
		int  new_capacity = capacity + increase_step;
        
        T* new_particles = (T*) mxMalloc(sizeof(T)*new_capacity);
 		T** new_ID_2_address = (T**) mxMalloc(sizeof(T*)*new_capacity);
        memcpy(new_particles,particles,sizeof(T)*capacity);
        memcpy(new_ID_2_address,ID_2_address,sizeof(T*)*capacity);        
		fprintf(stderr,"   now %i particles are allocated \n",new_capacity);
        
        int  i;
        for (i = 0; i < capacity; i++)
		{
			new_ID_2_address[i] = new_ID_2_address[i] - particles + new_particles;   // shift address
		}
        
		for (i = capacity; i < new_capacity; i++)
		{
			new_particles[i].ID = i;		     // initialize unique IDs
			new_ID_2_address[i] = &(new_particles[i]) ;   // initliaze pointer in LUT
		}
        
        for (i = 0; i < pcnt; i++)
        {            
            Particle *p = &(particles[i]);
            if (p->active)
            {
                Particle *new_p = &(new_particles[i]);
                for (int k = 0; k<NUMGRIDS;k++)
                {
                    if (p->gridindex[k] != -1)
                        g1[k].grid[p->gridindex[k]] = new_p;
                }
            }
            
        }
        
        mxFree(particles);
        mxFree(ID_2_address);
        
		particles = new_particles;
		ID_2_address = new_ID_2_address;
		capacity = new_capacity;
        needs_reallocation = false;
		return 1;
	}


	
	~ParticleGrid()
	{
		if (particles != 0)
			mxFree(particles);        
		if (ID_2_address != 0)
			mxFree(ID_2_address);
        for (int k = 0; k < NUMGRIDS; k++)
            g1[k].freeMem();
	}


    //-------------- converts particle ID to index in array
	int ID_2_index(int ID)
	{
		if (ID == -1)
			return -1;
		else
			return (ID_2_address[ID] - particles);

	}


    //------------- return the number of particles
    int getNumParticles()
    {
        int cnt; 
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (PARTMOD)
        #endif
        {
            cnt = pcnt - freeslots.size();
        }
        return cnt;
    }
    
    
    
    //----------------- creates a new particle
	T* newParticle(pVector R, pVector N)
	{		
		// get free place in container;
		if (pcnt >= capacity)
		{
			fprintf(stderr,"capacity overflow , needs reallocation ...\n");
            needs_reallocation = true;
            return 0;
		}
        
        // continous positions in the grids
        pVector Rs[NUMGRIDS]; 
        Rs[0] = R+ N*len;
        Rs[1] = R- N*len;
        Rs[2] = R;
        
        // get the corresponding indices
        int idx[NUMGRIDS];
        for (int k = 0; k < NUMGRIDS; k++)
        {
            idx[k] = g1[k].getIndex(Rs[k]);
            if (idx[k] == -1)
                return 0;
        };
        
		if (g1[0].occnt[idx[0]] < g1[0].csize && g1[1].occnt[idx[1]] < g1[1].csize && g1[2].occnt[idx[2]] < g1[2].csize)
		{
            T *p;
            #ifdef PARALLEL_OPENMP
            #pragma omp critical (PARTMOD)
            #endif
            {
                if (freeslots.empty())
                {
                   p = &(particles[pcnt]);               
                   pcnt++;					             
                }
                else
                {
                    p = &(particles[freeslots.back()]);
                    freeslots.pop_back();
                }
                p->active = true;
                p->R = R;	
                p->mID = -1;
                p->pID = -1;	
                
                for (int k = 0; k < NUMGRIDS; k++)
                {
                    insertIntoGrid(k, p,idx[k]);
                }
            }
			return p;
		} 
		else
		{
			celloverflows++;
			fprintf(stderr,"error: cell overflow \n");
			return 0;
		}
	}



    //----------------- moves a  particle    
	inline bool tryUpdateGrid(int k)
	{		
        
		T* p = &(particles[k]);
		
		/////// find new grid cells
        
        pVector Rs[NUMGRIDS]; 
        Rs[0] = p->R+ p->N*len;
        Rs[1] = p->R- p->N*len;
        Rs[2] = p->R;
        
        int idx[NUMGRIDS];
        for (int j = 0; j < NUMGRIDS; j++)
        {
            idx[j] = g1[j].getIndex(Rs[j]);
            if (idx[j] == -1)
                 return false;
        }
        int cellidx[NUMGRIDS];
        for (int j = 0; j < NUMGRIDS; j++)
        {
            if (p->gridindex[j] != -1)
            {
                cellidx[j] = p->gridindex[j]/g1[j].csize;
                if (idx[j] != cellidx[j]) // cell has changed 
                {
                    if (!(g1[j].occnt[idx[j]] < g1[j].csize))
                        return false;
                }            
            }
        }
        for (int j = 0; j < NUMGRIDS; j++)
        {
            #ifdef PARALLEL_OPENMP
            #pragma omp critical (PARTMOD)
            #endif
            {
                if (p->gridindex[j] != -1)
                {
                    // remove from old position in grid;
                     removeFromGrid(j,p);

                    // insert at new position in grid
                     insertIntoGrid(j, p,idx[j]);
                }
            }	
		}
        return true;
	}

    

    //------------- insert p in jth grid    
    inline void insertIntoGrid(int j, T *p)
    {
        if (j==0)
        {
            pVector Rs = p->R+ p->N*len;
            int idx = g1[j].getIndex(Rs);
            insertIntoGrid(j, p,idx);
        }
        else if (j == 1)
        {
            pVector Rs = p->R- p->N*len;
            int idx = g1[j].getIndex(Rs);
            insertIntoGrid(j, p,idx);
        }
        else if (j == 2)
        {
            pVector Rs = p->R;
            int idx = g1[j].getIndex(Rs);
            insertIntoGrid(j, p,idx);
        }            
    }
    
    //------------- insert p in jth grid at cell idx
    inline void insertIntoGrid(int j, T *p,int idx)
    {
        p->gridindex[j] = idx*g1[j].csize + g1[j].occnt[idx];
        g1[j].grid[p->gridindex[j]] = p;
        g1[j].occnt[idx]++;
    }
    
    
    
    
    //-------------- removes a particle
    inline void remove (T* p)
    {
        int pnum = ID_2_index(p->ID);
        if (pnum != -1)
            remove(pnum);
	}
    
    
    //-------------- removes a particle
	inline void remove(int k)
	{		
		T* p = &(particles[k]);
        #ifdef PARALLEL_OPENMP
        #pragma omp critical (PARTMOD)
        #endif
        {
            // remove pending connections
            if (p->mID != -1)
                destroyConnection(p,-1);
            if (p->pID != -1)
                destroyConnection(p,+1);
    
            for (int j = 0; j < NUMGRIDS; j++)
            {
                 removeFromGrid(j,p);
            }

            p->active = false;
            freeslots.push_back(k);
        }

	}
    
    //------------- remove p with  p->gridindex[j] from j-th grid
    inline void removeFromGrid(int j,T *p)
    { 
        int grdindex = p->gridindex[j];
        int cellidx = grdindex/g1[j].csize;
        g1[j].grid[grdindex] = g1[j].grid[cellidx*g1[j].csize + g1[j].occnt[cellidx]-1];
        g1[j].grid[grdindex]->gridindex[j] = grdindex;
        g1[j].occnt[cellidx]--;
        p->gridindex[j] = -1;
    }
    
    
    //------------- defrags the particle list
    void defrag()
    {      
     
       int newpcnt = pcnt-freeslots.size();
       int k;
       for (k = pcnt-1; k >= 0; k--)
       {
            if (particles[k].active)
            {
                
                int free = -1;
                while (!freeslots.empty())
                {
                    free = freeslots.back();
                    freeslots.pop_back();
                    if (free > k)
                    {
                        free = -1;
                        continue;
                    }
                    else
                        break;                    
                }
                
                if (free == -1)               
                    break;
                
                Particle tmp = particles[free];
                particles[free] = particles[k];
                particles[k] = tmp;
                for (int j = 0; j < NUMGRIDS; j++)         
                {
                    if (particles[free].gridindex[j] != -1)
                        g1[j].grid[particles[free].gridindex[j]] = &(particles[free]);
                }
                ID_2_address[particles[free].ID] = &(particles[free]);
                ID_2_address[particles[k].ID] = &(particles[k]);
                                
            }
           
       }
       freeslots.clear();
       pcnt = newpcnt;
       
    }
    
    //------------- get cell iterator (only needed for grid 2, the 'voxel'-grid)
    T**  getCell(pVector R, int &cnt)
    {         
        int x = floor(R.x*g1[2].mulx); 
        int y = floor(R.y*g1[2].muly); 
        int z = floor(R.z*g1[2].mulz); 
         
        int idx = (x + g1[2].nx*(y+z*g1[2].ny))    ;  
        
        cnt = g1[2].occnt[idx];
        return &(g1[2].grid[g1[2].csize*idx]);
    }
    T**  getCell(int x,int y,int z, int &cnt)
    {         
         
        int idx = (x + g1[2].nx*(y+z*g1[2].ny))    ;  
        
        cnt = g1[2].occnt[idx];
        return &(g1[2].grid[g1[2].csize*idx]);
    }
     
    
    GridArray *getVoxelGridArray()
    {
            return &(g1[2]);
    }
        
    //------------- get neighborhood iterator (only needed for grid 0,1, the 'conn'-grids)
	inline void computeNeighbors(pVector &R,NeighborTracker &nbtrack, int ep)
	{
        int idx;
        if (ep == 1)
            idx = 0;
        else
            idx = 1;
        
		REAL xfrac = R.x*g1[idx].mulx;
		REAL yfrac = R.y*g1[idx].muly;
		REAL zfrac = R.z*g1[idx].mulz;
		int xint = int(floor(xfrac));
		int yint = int(floor(yfrac));
		int zint = int(floor(zfrac));

		int dx = -1;
		if (xfrac-xint > 0.5) dx = 1;
		if (xint <= 0) { xint = 0; dx = 1; }
		if (xint >= g1[idx].nx-1) { xint = g1[idx].nx-1; dx = -1; }

		int dy = -1;
		if (yfrac-yint > 0.5) dy = 1;
		if (yint <= 0) {yint = 0; dy = 1; }
		if (yint >= g1[idx].ny-1) {yint = g1[idx].ny-1; dy = -1;}

		int dz = -1;
		if (zfrac-zint > 0.5) dz = 1;
		if (zint <= 0) {zint = 0; dz = 1; }
		if (zint >= g1[idx].nz-1) {zint = g1[idx].nz-1; dz = -1;}


		nbtrack.cellidx[0] = xint + g1[idx].nx*(yint+zint*g1[idx].ny);
		nbtrack.cellidx[1] = nbtrack.cellidx[0] + dx;
		nbtrack.cellidx[2] = nbtrack.cellidx[1] + dy*g1[idx].nx;
		nbtrack.cellidx[3] = nbtrack.cellidx[2] - dx;		
		nbtrack.cellidx[4] = nbtrack.cellidx[0] + dz*g1[idx].nx*g1[idx].ny;
		nbtrack.cellidx[5] = nbtrack.cellidx[4] + dx;
		nbtrack.cellidx[6] = nbtrack.cellidx[5] + dy*g1[idx].nx;
		nbtrack.cellidx[7] = nbtrack.cellidx[6] - dx;


		nbtrack.cellidx_c[0] = g1[idx].csize*nbtrack.cellidx[0];
		nbtrack.cellidx_c[1] = g1[idx].csize*nbtrack.cellidx[1];
		nbtrack.cellidx_c[2] = g1[idx].csize*nbtrack.cellidx[2];
		nbtrack.cellidx_c[3] = g1[idx].csize*nbtrack.cellidx[3];
		nbtrack.cellidx_c[4] = g1[idx].csize*nbtrack.cellidx[4];
		nbtrack.cellidx_c[5] = g1[idx].csize*nbtrack.cellidx[5];
		nbtrack.cellidx_c[6] = g1[idx].csize*nbtrack.cellidx[6];
		nbtrack.cellidx_c[7] = g1[idx].csize*nbtrack.cellidx[7];



		nbtrack.cellcnt = 0;
		nbtrack.pcnt = 0;
        nbtrack.gnum = idx;

	}
    
    //-------- to iterate through the neighborhood
	inline int getNeighborCount(NeighborTracker &nbtrack)
    {
        int idx = nbtrack.gnum;
        
        int cnt = 0;
        for (int k = 0; k < 8; k++)            
            cnt += g1[idx].occnt[nbtrack.cellidx[k]];
        return cnt;
    }    
	inline T *getNextNeighbor(NeighborTracker &nbtrack)
	{
		
		if (nbtrack.pcnt < g1[nbtrack.gnum].occnt[nbtrack.cellidx[nbtrack.cellcnt]])
		{
			
			return g1[nbtrack.gnum].grid[nbtrack.cellidx_c[nbtrack.cellcnt] + (nbtrack.pcnt++)];

		}
		else
		{
					
			for(;;)
			{
				nbtrack.cellcnt++;
				if (nbtrack.cellcnt >= 8)
					return 0;
				if (g1[nbtrack.gnum].occnt[nbtrack.cellidx[nbtrack.cellcnt]] > 0)
					break;
			}
		
			nbtrack.pcnt = 1;
			return g1[nbtrack.gnum].grid[nbtrack.cellidx_c[nbtrack.cellcnt]];
		}
    }


    //-------- create a connection between two particle at two specific ports
	inline void createConnection(T *P1,int ep1, T *P2, int ep2)
	{
		if (ep1 == -1) 
        {
			P1->mID = P2->ID;            
            removeFromGrid(1,P1);  
            
        }
		else
        {
			P1->pID = P2->ID;
           removeFromGrid(0,P1);            
            
        }
		if (ep2 == -1) 
        {
			P2->mID = P1->ID;
            removeFromGrid(1,P2);
        }
		else
        {
			P2->pID = P1->ID;
            removeFromGrid(0,P2);
        }

		concnt++;
	}

    //-------- destroys a connection between two particle at two specific ports
	inline void destroyConnection(T *P1,int ep1, T *P2, int ep2)
	{
		if (ep1 == -1) 
        {
			P1->mID = -1;
            insertIntoGrid(1,P1);
        }
		else
        {
			P1->pID = -1;
            insertIntoGrid(0,P1);
        }
        
		if (ep2 == -1) 
        {
			P2->mID = -1;
            insertIntoGrid(1,P2);
        }
		else
        {
			P2->pID = -1;
            insertIntoGrid(0,P2);
        }
		concnt--;
	}

	inline void destroyConnection(T *P1,int ep1)
	{

		T *P2 = 0;
		int ep2;
		if (ep1 == 1)
		{
			P2 = ID_2_address[P1->pID];
			P1->pID = -1;
            insertIntoGrid(0,P1);
		}
		else		
		{
			P2 = ID_2_address[P1->mID];
			P1->mID = -1;
            insertIntoGrid(1,P1);
		}
		if (P2->mID == P1->ID)
		{
			P2->mID = -1;
            insertIntoGrid(1,P2);
		}
		else		
		{
			P2->pID = -1;
            insertIntoGrid(0,P2);
		}
		concnt--;

	}
    
    
    
};



#endif

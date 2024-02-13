"""
Created on 03.02.2017

@author: Christian Schmidt cschmidt18057 (at) gmail.com
Copyright 2017
This file is part of FanPy
"""

import numpy as np
import copy

class Axon():
    """
    Axon(params)
    
    Contains the axon parameters for the axon model based on the myelinated axon
    model by C C McIntyre:
    
    McIntyre, C C and Richardson, A G and Grill, W M. Modelling the excitability of
    mammalian nerve fibers: influence of afterpotentials on the recovery cycle.
    J Neurophysiol, vol. 87, pp. 995-1006, Feb. 2002
    
    Parameters for 2 um and 3 um fiber diameter are taken from:
    Sotiropoulos, S N  and Steinmetz, P N. Assessing the direct effects of deep
    brain stimulation using embedded axon models. J Neural Eng, vol. 4,
    pp. 107-19, Jun 2007.
    
    Parameters
    ----------
    params : dict 
        'centered' : (True, False)
            If True, the axon reference nodes are centered to the center compartment
        'diameter' : float
            Diameter of the axonal fiber. Supported diameters are available by AXONPARAMETERS key.
       
    """  
    
    # condg for fiber diameter 2.0 um and 3.0 um is determined using the polyfit from the
    # data 5.7 um and larger. Resulting fitting function: 0.0001x^2+0.0139x+0.5235
    AXONPARAMETERS = {
            2.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 60, 'deltax': 200.0, 'ranvier_length': 1.0,  
                'para1_length': 3.0, 'para2_length': 10.0, 'axon_diameter': 1.6,
                'node_diameter': 1.4, 'para1_diameter': 1.4, 'para2_diameter': 1.6,
                'condg': 0.552, 'lamellas': 30
                },
            3.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 60, 'deltax': 278.0, 'ranvier_length': 1.0,  
                'para1_length': 3.0, 'para2_length': 17.0, 'axon_diameter': 2.1,
                'node_diameter': 1.52, 'para1_diameter': 1.52, 'para2_diameter': 2.1,
                'condg': 0.567, 'lamellas': 43
                },
            #(fiberD==5.7) {g=0.605 axonD=3.4 nodeD=1.9 paraD1=1.9 paraD2=3.4 deltax=500 paralength2=35 nl=80}
            5.7: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 500.0, 'ranvier_length': 1.0, 
                'para1_length': 3.0, 'para2_length': 35.0, 'axon_diameter': 3.4,
                'node_diameter': 1.9, 'para1_diameter': 1.9, 'para2_diameter': 3.4,
                'condg': 0.605, 'lamellas': 80
                },
            #(fiberD==7.3) {g=0.630 axonD=4.6 nodeD=2.4 paraD1=2.4 paraD2=4.6 deltax=750 paralength2=38 nl=100}
            7.3: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 750.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 38.0, 'axon_diameter': 4.6,
                'node_diameter': 2.4, 'para1_diameter': 2.4, 'para2_diameter': 4.6,
                'condg': 0.630, 'lamellas': 100
                },
            #(fiberD==8.7) {g=0.661 axonD=5.8 nodeD=2.8 paraD1=2.8 paraD2=5.8 deltax=1000 paralength2=40 nl=110}
            8.7: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1000.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 40.0, 'axon_diameter': 5.8,
                'node_diameter': 2.8, 'para1_diameter': 2.8, 'para2_diameter': 5.8,
                'condg': 0.661, 'lamellas': 110
                },
            #(fiberD==10.0) {g=0.690 axonD=6.9 nodeD=3.3 paraD1=3.3 paraD2=6.9 deltax=1150 paralength2=46 nl=120}
            10.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1150.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 46.0, 'axon_diameter': 6.9,
                'node_diameter': 3.3, 'para1_diameter': 3.3, 'para2_diameter': 6.9,
                'condg': 0.690, 'lamellas': 120
                },
            #(fiberD==11.5) {g=0.700 axonD=8.1 nodeD=3.7 paraD1=3.7 paraD2=8.1 deltax=1250 paralength2=50 nl=130}
            11.5: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1250.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 50.0, 'axon_diameter': 8.1,
                'node_diameter': 3.7, 'para1_diameter': 3.7, 'para2_diameter': 8.1,
                'condg': 0.700, 'lamellas': 130
                },
            #(fiberD==12.8) {g=0.719 axonD=9.2 nodeD=4.2 paraD1=4.2 paraD2=9.2 deltax=1350 paralength2=54 nl=135}
            12.8: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1350.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 54.0, 'axon_diameter': 9.2,
                'node_diameter': 4.2, 'para1_diameter': 4.2, 'para2_diameter': 9.2,
                'condg': 0.719, 'lamellas': 135
                },
            #(fiberD==14.0) {g=0.739 axonD=10.4 nodeD=4.7 paraD1=4.7 paraD2=10.4 deltax=1400 paralength2=56 nl=140}
            14.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1400.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 56.0, 'axon_diameter': 10.4,
                'node_diameter': 4.7, 'para1_diameter': 4.7, 'para2_diameter': 10.4,
                'condg': 0.739, 'lamellas': 140
                },
            #(fiberD==15.0) {g=0.767 axonD=11.5 nodeD=5.0 paraD1=5.0 paraD2=11.5 deltax=1450 paralength2=58 nl=145}
            15.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1450.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 58.0, 'axon_diameter': 11.5,
                'node_diameter': 5.0, 'para1_diameter': 5.0, 'para2_diameter': 11.5,
                'condg': 0.767, 'lamellas': 145
                },
            #(fiberD==16.0) {g=0.791 axonD=12.7 nodeD=5.5 paraD1=5.5 paraD2=12.7 deltax=1500 paralength2=60 nl=150}
            16.0: {'ranvier_nodes': 21, 'para1_nodes': 40, 'para2_nodes': 40,
                'inter_nodes': 120, 'deltax': 1500.0, 'ranvier_length': 1.0,
                'para1_length': 3.0, 'para2_length': 60.0, 'axon_diameter': 12.7,
                'node_diameter': 5.5, 'para1_diameter': 5.5, 'para2_diameter': 12.7,
                'condg': 0.791, 'lamellas': 150
                },
            }
    
    def __init__(self, params):
        self.__params = params
        if params is None:
            self.__params = dict()
        if not 'centered' in self.__params:
            self.__params['centered'] = True
        if not 'diameter' in params:
            self.__params['diameter'] = 5.7
        
        self.__axonparams = self.__create_axon_parameters(self.__params['diameter'])            
        self.__ref_nodes = self.__create_ref_nodes()
        self.__create_nodes()
        
    def get_axonparams(self):
        """
        Returns the axonparameter as given by paper mentioned in the class documentation.
        
        Returns
        -------
        dict
            'total_nodes'    int
                Total number of compartments
            'ranvier_nodes'  int
                Number of node of Ranvier compartments
            'para1_nodes'    int
                Number of first paranodal compartments
            'para2_nodes'    int
                Number of second paranodal compartments
            'inter_nodes'    int
                Number of internodal compartments
            'ranvier_length' float
                Length of the Node of Ranvier
            'para1_length'   float
                Length of the first paranodal compartments
            'para2_length'   float
                Length of the second paranodal compartments
            'inter_length'   float
                Length of the internodal compartments
            'deltax'         float
                Length from a node of Ranvier to the next
            'fiberD'         float
                Fiber diameter
            'axon_diameter': float
                Diameter of the axon
            'node_diameter': float
                Diameter of the nodes of Ranvier
            'para1_diameter': float
                Diameter of the first paranodal compartments
            'para2_diameter': float,
                Diameter of the second paranodal compartments
            'condg':          float
                Axial conductivity
            'lamellas':       int
                Number of lamellas in the myelin sheath
            
        """ 
        return self.__axonparams

    def __create_ref_nodes(self):
        """ Topology: Node,MYSA,FLUT,STIN(6x),FLUT,MYSA;Node,... : 11 in total
        Determine the location of the nodes along the axon at which
        the extracellular potential will be applied.
        """
        n = self.__axonparams['total_nodes']
        l_ranvier = self.__axonparams['ranvier_length']
        l_para1 = self.__axonparams['para1_length']
        l_para2 = self.__axonparams['para2_length']
        l_inter = self.__axonparams['inter_length']
     
        nranvier = self.__axonparams['ranvier_nodes']
        nstin = int(float(self.__axonparams['inter_nodes'])/(nranvier-1))
        nmysa = int(float(self.__axonparams['para1_nodes'])/(nranvier-1))
        nflut = int(float(self.__axonparams['para2_nodes'])/(nranvier-1))
        ncomp = nstin+nmysa+nflut+1      

        ref_nodes = np.zeros(n,)
        for i in range(0, n):
            j = i % ncomp
            if j == 0 or j == 1:  # mysa <-> ranvier node
                if i == 0:  # first of all nodes
                    ref_nodes[i] = l_ranvier/2.
                else:
                    ref_nodes[i] = ref_nodes[i-1]+l_para1/2.+l_ranvier/2.
            elif (j > 1 and j < (int(nmysa/2)+1)) or (j > (int(nmysa/2)+nflut+nstin) and j <(nmysa+nflut+nstin)): # mysa <-> mysa node
                ref_nodes[i] = ref_nodes[i-1]+l_para1
            elif j == (int(nmysa/2)+1) or j == (int(nmysa/2)+nflut+nstin):  # mysa <-> flut node
                ref_nodes[i] = ref_nodes[i-1]+l_para1/2.+l_para2/2.
            elif (j > (int(nmysa/2)+1) and j < (int(nmysa/2)+int(nflut/2)+1)) or (j > (int(nmysa/2)+int(nflut/2)+nstin) and j < (int(nmysa/2)+nflut+nstin)): # flut <-> flut node
                ref_nodes[i] = ref_nodes[i-1]+l_para2
            elif j == (int(nmysa/2)+int(nflut/2)+1) or j == (int(nmysa/2)+int(nflut/2)+nstin):  # flut <-> stin node
                ref_nodes[i] = ref_nodes[i-1]+l_para2/2.+l_inter/2.
            else:  # stin <-> stin node
                ref_nodes[i] = ref_nodes[i-1]+l_inter

        return ref_nodes * 1e-6  # um -> m

    def __create_nodes(self):
        if self.__params['centered']:
            ref_nodes = self.__ref_nodes - self.__ref_nodes[0]/2 - \
                                                       self.__ref_nodes[-1]/2
        else:
            ref_nodes = self.__ref_nodes
        # standard position(along x axis in [x,y,z])
        nodes = np.zeros((len(ref_nodes),3))
        nodes[:,0] = ref_nodes
        nodes[:,1] = np.zeros(len(ref_nodes),)
        nodes[:,2] = np.zeros(len(ref_nodes),)
        self.__nodes = nodes
        return self.__nodes
    
    def get_ref_nodes(self):
        """  
        Reference locations of the center of the compartments of the axon
    
        Returns
        -------
        numpy array
            1D array of the reference locations
    
        """
        return self.__ref_nodes

    def get_nodes(self):
        """
        Locations of the center of the compartments of the axon.
        The axon nodes are aligned along the x-axis (normal to yz-plane).
    
        Returns
        -------
        numpy array
            3D array of the reference locations in cartesian coordinates
    
        """   
        return self.__nodes
    
    def __create_axon_parameters(self, diameter):
        if diameter not in self.AXONPARAMETERS:
            raise Exception ("Couldn't find axon with diameter "+str(diameter)+" in axon parameter lookup map.")
        
        axon_parameters = copy.deepcopy(self.AXONPARAMETERS[diameter])
        
        nranvier = axon_parameters['ranvier_nodes']
        nstin = int(float(axon_parameters['inter_nodes'])/(nranvier-1))
        nmysa = int(float(axon_parameters['para1_nodes'])/(nranvier-1))
        nflut = int(float(axon_parameters['para2_nodes'])/(nranvier-1))
        
        axon_parameters['fiberD']=diameter
        axon_parameters['total_nodes'] = axon_parameters['ranvier_nodes'] + \
            axon_parameters['para1_nodes']+axon_parameters['para2_nodes'] + \
            axon_parameters['inter_nodes']
        axon_parameters['inter_length'] = (axon_parameters['deltax'] -
                                       axon_parameters['ranvier_length'] -
                                       nmysa*axon_parameters['para1_length'] -
                                       nflut*axon_parameters['para2_length'])/float(nstin)
        return axon_parameters
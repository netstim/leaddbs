"""
Created on 30.01.2017

@author: Christian Schmidt cschmidt18057 (at) gmail.com
Copyright 2017
This file is part of FanPy
"""

import numpy as np
from enum import Enum


class DielectricProperties:
    """
    DielectricProperties(tissue)
    
    Dielectric properties of biological tissue based on a 4-term ColeCole model
    with parameters determined in:
    Gabriel S, Lau R W, Gabriel C, 1996, The dielectric properties of
    biological tissue: III. Parametric models for the dielectric spectrum of
    tissues, Phys Med Biol, vol. 41, pp. 2271-2293    
    
    Parameters
    ----------
    tissue : ETissue 
        name of the biological tissue

    """

    FREQUENCIES = 'frequencies'
    CONDUCTIVITY = 'conductivity'
    PERMITTIVITY = 'permittivity'

    def __init__(self, tissue):
        self.tissue = tissue
        self.params = ColeColeParams(tissue)
        self.poles = 4
        
    def set_poles(self, poles):
        """
        Number of poles for the ColeCole model.
        
        Parameters
        ----------
        poles : int 
            1 to 4 poles, beginning with the pole with lowest frequency
            default = 1

        """
        if poles > 0 and poles < 5:
            self.poles = poles
        else:
            print('Number of relaxation poles has to be in [1,4]')

    def get_dielectrics(self, frequencies):
        """
        Dielectric tissue properties for the given frequency (frequencies)
        
        Parameters
        ----------
        frequencies : float or numpy.array 
            Frequency or frequencies at which the dielectric should be determined
            
        Returns
        -------
        dict
            FREQUENCIES: float or numpy.array
                frequencies to the corresponding dielectric properties
            CONDUCTIVITY: float or numpy.array
                conductivity values
            PERMITTIVITY: float or numpy.array
                relative permittivity values
    
        """
        isVector = isinstance(frequencies, (list, tuple, np.ndarray))
        freq = frequencies
        if isVector is False:
            freq = np.array([frequencies])
        
        idxs_zero = np.where(frequencies == 0)
        idxs_nonzero = np.where(frequencies != 0)
        e0 = 8.854e-12
        su = np.zeros(len(freq),)
        # order of the parameter vectors is flipped to be sorted from
        # low frequency pole to high frequency pole
        for i in range(3, (4-self.poles)-1, -1):
            su = su + \
                self.params.de[i]/(1 +
                                   (2j*np.pi*freq*self.params.t[i]) **
                                   (1-self.params.a[i]))

        # treat DC (f=0) components
        # self.params.de is flipped to be sorted from low frequency pole to
        # high frequency pole
        perm_dc = np.sum(self.params.de[::-1][0:self.poles])+self.params.einf
        cond_dc = self.params.sig
        
        # generate complex permittivity vector
        permc = self.params.einf+ \
            su[idxs_nonzero]+self.params.sig/(2j*np.pi*freq[idxs_nonzero]*e0)
        
        # generate permittivity vector
        perm = np.zeros(len(freq),)
        perm[idxs_nonzero] = permc.real
        perm[idxs_zero] = perm_dc
        
        # generate conductivity vector
        cond = np.zeros(len(freq),)
        cond[idxs_nonzero] = -permc.imag*2*np.pi*freq[idxs_nonzero]*e0
        cond[idxs_zero] = cond_dc
        
        if isVector is False:
            cond = cond[0]
            perm = perm[0]
            freq = freq[0]
        reslt=[cond,perm]
        #return reslt
        return cond,perm
        #return {
        #        self.FREQUENCIES: freq,
        #        self.CONDUCTIVITY: cond,
        #        self.PERMITTIVITY: perm
        #        }


class ColeColeParams:
    """
    ColeColeParams(tissue)
    
    Parameters of the ColeCole model for the given tissue.
    
    Parameters
    ----------
    tissue : ETissue 
        name of the biological tissue

    """

    def __init__(self, tissue):

        if tissue == 3:
            self.name = 'Brain (Grey Matter)'
            self.einf = 4.000
            self.sig = 0.020
            self.de = np.array([45.00, 400, 2.00e+5, 4.50e+7])
            self.a = np.array([0.100, 0.150, 0.220, 0.000])
            self.t = np.array([7.958e-12, 15.915e-9, 106.103e-6, 5.305e-3])
        elif tissue == 2:
            self.name = 'Brain (White Matter)'
            self.einf = 4.000
            self.sig = 0.020
            self.de = np.array([32.00, 100, 4.00e+4, 3.50e+7])
            self.a = np.array([0.100, 0.100, 0.300, 0.020])
            self.t = np.array([7.958e-12, 7.958e-9, 53.052e-6, 7.958e-3])
        elif tissue == 1:
            self.name = 'Cerebro Spinal Fluid'
            self.einf = 4.000
            self.sig = 2.000
            self.de = np.array([65.00, 40, 0.00e+0, 0.00e+0])
            self.a = np.array([0.100, 0.000, 0.000, 0.000])
            self.t = np.array([7.958e-12, 1.592e-9, 159.155e-6, 15.915e-3])

        else:
            raise NotImplementedError('Unknown tissue type ', tissue)

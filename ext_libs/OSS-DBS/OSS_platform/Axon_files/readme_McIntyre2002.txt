/*--------------------------------------------------------------------
Based on MRGaxon.hoc from
NEURON Yale ModelDB: Spinal Motor Neuron (McInytre et al 2002
https://senselab.med.yale.edu/modeldb/showModel.cshtml?model=3810&file=/MRGaxon/MRGaxon.hoc#tabs-1

2/02
Cameron C. McIntyre
SIMULATION OF PNS MYELINATED AXON

This model is described in detail in:

McIntyre CC, Richardson AG, and Grill WM. Modeling the excitability of
mammalian nerve fibers: influence of afterpotentials on the recovery
cycle. Journal of Neurophysiology 87:995-1006, 2002.

This model can not be used with NEURON v5.1 as errors in the
extracellular mechanism of v5.1 exist related to xc. The original
stimulations were run on v4.3.1. NEURON v5.2 has corrected the 
limitations in v5.1 and can be used to run this model.

Modifications:
01/17
Christian Schmidt
Apply external extracellular potential to the sections of the model.
----------------------------------------------------------------------*/
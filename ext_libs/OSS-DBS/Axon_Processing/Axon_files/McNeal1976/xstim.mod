: $Id$

COMMENT
This mechanism is intended to be used in conjunction 
with the extracellular mechanism.  Pointers specified 
at the hoc level must be used to connect the 
extracellular mechanism's e_extracellular 
to this mechanism's ex.

xstim does two useful things:

1. Serves as a target for Vector.play() to facilitate 
extracellular stimulation.  Assumes that one has initialized 
a Vector to hold the time sequence of the stimulus current.
This Vector is to be played into the GLOBAL variable is 
(GLOBAL so only one Vector.play() needs to be executed), 
which is multiplied by the RANGE variable rx ("transfer 
resistance between the stimulus electrode and the local 
node").  This product, called ex in this mechanism, is the 
extracellular potential at the local node, i.e. is used to 
drive local e_extracellular.

2. Allows local storage of xyz coordinates interpolated from 
the pt3d data.  These coordinates are used by hoc code that 
computes the transfer resistance that couples the membrane 
to extracellular stimulating and recording electrodes.

This mechanism uses the BEFORE BREAKPOINT block 
to ensure that the stimulus potential is computed 
prior to the solution step.

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = is*rx*(1e6)
}

ENDCOMMENT

NEURON {
  SUFFIX xstim
  POINTER ex
  GLOBAL is
  RANGE rx
  RANGE x, y, z
}

PARAMETER {
  : default transfer resistance between stim electrodes and axon
  rx = 1 (megohm) : mV/nA
  x = 0 (1) : spatial coords
  y = 0 (1)
  z = 0 (1)
}

ASSIGNED {
  v (millivolts)
  ex (millivolts)
  is (milliamp)
}

INITIAL {
  ex = is*rx*(1e6)
}

BEFORE BREAKPOINT { : before each cy' = f(y,t) setup
  ex = is*rx*(1e6)
}


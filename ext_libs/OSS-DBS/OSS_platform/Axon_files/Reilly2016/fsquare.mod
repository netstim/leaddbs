COMMENT
square.mod
Generates a square wave.
User specifies
del  time at which first cycle starts
dp   duration of a phase (duration of a half cycle)
num  number of cycles
amp1 level for first half cycle
amp2 level for second half cycle
: bal  nonzero forces amp2 = -amp1
20150417 NTC
ENDCOMMENT

NEURON {
  POINT_PROCESS Fsquare
:  RANGE del, dp, num, amp1, amp2, bal
  RANGE del, dp, num, amp1, amp2
  POINTER x
}

PARAMETER {
  del = 0 (ms) <0, 1e9> : time at which first cycle starts
  dp = 0 (ms) <0, 1e9> : phase duration (half cycle duration)
  num = 0 (1) : how many cycles
  amp1 = 0 (1) : level for first half cycle
  amp2 = 0 (1) : level for second half cycle
:  bal = 0 (1) : nonzero forces amp2 = -amp1
}

ASSIGNED {
  x (1)
  on (1)
  tally (1) : how many more cycles are to be generated
}

UNITSOFF
FUNCTION nonneg(x) {
  nonneg = x
  if (x<0) {
    nonneg = 0
  } else {
    nonneg = x
  }
}

INITIAL {
  on = 0
  x = 0
  : force these to be nonnegative
  del = nonneg(del)
  dp = nonneg(dp)
  num = nonneg(num)

  : do nothing if num == 0 or dp == 0
  if (num*dp>0) {
    tally = num
    net_send(del,1) : to start first phase
  }

:  if (bal!=0) {
:    amp2 = -amp1
:  }
}
UNITSON

NET_RECEIVE (w) {
  : respond only to self-events
  if (tally>0) { : generate output until all cycles have been completed
    if (flag == 1) { : start a new cycle
      if (on == 0) {
        on = 1 : signal that it's "on"
      }
      x = amp1 : enter phase 1
      : prepare for phase 2
      net_send(dp, 2)
    }
    if (flag == 2) {
      x = amp2 : enter phase 2
      tally = tally - 1
      if (tally>0) {
        net_send(dp, 1) : prepare for next cycle
      } else {
        net_send(dp, 3) : end of waveform
      }
    }
  }
  if (flag == 3) { : no more cycles to generate
    on = 0 : signal that it's "off"
    x = 0
  }
}


: $Id: fzap.mod,v 1.5 2015/05/28 02:51:15 ted Exp ted $

COMMENT
fzap.mod

A bogus point process that contains the variable x, 
which oscillates starting at t = del >= 0.
The frequency f of the oscillation increases linearly with time
from f0 at t == del to f1 at t == del + dur, 
where both del and dur are > 0.

fzap uses the event delivery system to ensure compatibility with adaptive integration.

=================
NOTES AND CAVEATS
=================

1.  If x were a RANGE variable, an assignment statement would 
have to be inserted into proc advance() in order for the 
value of x to be used by other mechanisms--e.g.
proc advance() {
  is_xtra = Fzap[0].x
  fadvance()
}
However, that would be incompatible with adaptive integration.
To eliminate the need for such an assignment statement, x is a 
POINTER.  This preserves compatibility with adaptive integration.

2.  On every fadvance, the statements that evaluate Fzap's x 
should be executed before the statements in any client mechanism 
that relies on the value of Fzap's x.  To that end, the value of 
x is computed in a BEFORE BREAKPOINT block, which will take care
of any client mechanism that uses Fzap's x in a BREAKPOINT block.

However, some client mechanisms may have their own 
BEFORE BREAKPOINT blocks that need the value of Fzap's x.  
xtra is such a mechanism.  In this situation, care is required 
to ensure that the statements in Fzap's BEFORE BREAKPOINT block
are executed first.  This can be done by compiling the mod file 
that defines Fzap _before_ the client mechanism's mod file.

There are two ways to make this happen:
A.  Invoke nrnivmodl with a command line that presents the file 
names in the desired sequence.  UNIX/Linux users may be quite 
comfortable with this.
B.  Choose mod file names so that Fzap's mod file appears before 
the name of any client mod files in an alphabetical listing.
For the example of Fzap and xtra, the file names fzap.mod and 
xtra.mod would be quite suitable.  This is more convenient for 
users of all operating systems, but especially MSWin and OS X, 
whose users are accustomed to compiling all mod files in a 
directory with mknrndll or "drag and drop," respectively.

12/11/2008 NTC
ENDCOMMENT

NEURON {
  POINT_PROCESS Fzap
  RANGE del, dur, f0, f1, amp, f
  POINTER x
}

UNITS {
  PI = (pi) (1)
}

PARAMETER {
  del (ms)
  dur (ms)
  f0 (1/s)  : frequency is in Hz
  f1 (1/s)
  amp (1)
}

ASSIGNED {
  f (1/s)
  x (1)
  on (1)
}

INITIAL {
  f = 0
  x = 0
  on = 0

  if (del<0) { del=0 }
  if (dur<0) { dur=0 }
  if (f0<=0) { f0=0 (1/s) }
  if (f1<=0) { f1=0 (1/s) }

  : do nothing if dur == 0
  if (dur>0) {
    net_send(del, 1)  : to turn it on and start frequency ramp
  }
}

COMMENT
The angular velocity in radians/sec is w = 2*PI*f, 
where f is the instantaneous frequency in Hz.

Assume for the moment that the frequency ramp starts at t = 0.
f = f0 + (f1 - f0)*t/dur

Then the angular displacement is
theta = 2*PI * ( f0*t + (f1 - f0)*(t^2)/(2*dur) ) 
      = 2*PI * t * (f0 + (f1 - f0)*t/(2*dur))
But the ramp starts at t = del, so just substitute t-del for every occurrence of t
in the formula for theta.
ENDCOMMENT

BEFORE BREAKPOINT {
  if (on==0) {
    f = 0
    x = 0
  } else {
    f = f0 + (f1 - f0)*(t-del)/dur
    x = amp * sin( 2*PI * (t-del) * (f0 + (f1 - f0)*(t-del)/(2*dur)) * (0.001) )
  }
}

NET_RECEIVE (w) {
  : respond only to self-events with flag > 0
  if (flag == 1) {
    if (on==0) {
      on = 1  : turn it on
      net_send(dur, 1)  : to stop frequency ramp, freezing frequency at f1
    } else {
      on = 0  : turn it off
    }
  }
}

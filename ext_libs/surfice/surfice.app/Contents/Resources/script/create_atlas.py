import gl
gl.resetdefaults()
gl.azimuthelevation(240, 15)
gl.	azimuthelevation(240, 15);
gl.meshload('AICHAhr.lh.mz3');
#format: atlashide(layer, (regions))
# gl.atlashide(0, (1, 2, 3, 12, 22, 44));
n = gl.atlasmaxindex(0)
msk = ()
for i in range(n):
  msk = msk + (i,)
  gl.atlashide(0,msk)
  gl.wait(20);
gl.atlashide(0,())

import gl
gl.resetdefaults()
kframesperrotation = 180
gl.meshload('BrainMesh_ICBM152_smoothed.lh.mz3')
gl.meshcurv()
gl.shadername('hidecurves')
gl.overlayload('CIT168.mz3')
gl.azimuthelevation(250, 35)
gl.meshhemisphere(-1)
for i in range(1, kframesperrotation * 5):
  if ((i % kframesperrotation) == 0):
    rot = (i / kframesperrotation)
    if rot == 1:
        gl.shadername('outline')
    elif rot == 2:
        gl.shadername('wire')
    elif rot == 3:
        gl.shadername('squares')
        gl.shaderambientocclusion(0)
    else :
        gl.shadername('xGaps')
  gl.azimuth( int(round(360.0/kframesperrotation)))
  gl.wait(20)
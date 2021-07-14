import gl
gl.resetdefaults()
kframesperrotation = 180
gl.meshload('BrainMesh_ICBM152_smoothed.lh.mz3')
gl.meshcurv()
gl.overlayload('CIT168.mz3')
gl.azimuthelevation(250, 35)
gl.meshhemisphere(-1)
for i in range(0, kframesperrotation * 7):
  rot = (i / kframesperrotation)
  if ((i % kframesperrotation) == 0):
    shader = 'Default'
    if rot == 0:
        shader = 'hidecurves'
    elif rot == 1:
        shader = 'outline'
    elif rot == 2:
        shader = 'wire'
    elif rot == 3:
        shader = 'squares'
    elif rot == 4:
        shader = 'xGaps'
    else :
        shader = 'MatCap'
    gl.shadername(shader)
    print('shader: ' + shader)
  azimuth = round((i * 360.0/kframesperrotation) % 360)
  xrayObj= 1.0
  if (rot == 5):
    xrayObj = abs(180 - azimuth) / 180
  xrayOver= 0.0
  if (rot == 6):
    xrayOver = 1.0 - abs(180 - azimuth) / 180
  gl. shaderxray(xrayObj, xrayOver);
  gl.azimuthelevation(int(azimuth), 35);
  gl.wait(20)
  #if (azimuth == 270):
  #  gl.savebmp(str(rot)+'.png')
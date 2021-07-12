import gl
gl.resetdefaults()
gl.azimuthelevation(240, 15)
gl.meshload('AICHAhr.lh.mz3')
#gl.atlashide(0,1, 2, 3)
gl.atlasstatmap('AICHAhr.lh.mz3','surficetemp.mz3',(4,5,6),(7,3,4))
gl.meshload('lh.pial');
gl.overlayload('surficetemp.mz3');
gl.shaderxray(1.0, 0.3);



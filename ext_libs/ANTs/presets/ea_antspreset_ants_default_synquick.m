function apref = ea_antspreset_ants_default_synquick(cmd)

if ischar(cmd) % 'query' mode, return the name of the preset
    apref = 'ANTs Default: SyN Quick';
    return
end

apref.rigid.gradientstep = '0.1';
apref.rigid.metric = 'MI';
apref.rigid.metricparams = '1,32,Regular,0.25';
apref.rigid.convergence='[1000x500x250x0,1e-6,10]'; % Rigid convergence params
apref.rigid.shrinkfactors='12x8x4x2'; % Rigid shrink factors
apref.rigid.smoothingsigmas='4x3x2x1vox'; % Rigid Smoothing sigmas

apref.affine.gradientstep = '0.1';
apref.affine.metric = 'MI';
apref.affine.metricparams = '1,32,Regular,0.25';
apref.affine.convergence='[1000x500x250x0,1e-6,10]'; % Affine convergence params
apref.affine.shrinkfactors='12x8x4x2'; % Affine shrink factors
apref.affine.smoothingsigmas='4x3x2x1vox'; % Affine Smoothing sigmas

apref.syn.gradientstep = '0.1,3,0';
apref.syn.metric = 'MI';
apref.syn.metricparams = '1,32';
apref.syn.convergence='[100x100x70x50x0,1e-6,10]'; % SyN convergence params
apref.syn.shrinkfactors='10x6x4x2x1'; % SyN shrink factors
apref.syn.smoothingsigmas='5x3x2x1x0vox'; % SyN Smoothing sigmas

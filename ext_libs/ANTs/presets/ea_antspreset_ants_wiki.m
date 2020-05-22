function apref = ea_antspreset_ants_wiki(cmd)
% from https://github.com/ANTsX/ANTs/wiki/Anatomy-of-an-antsRegistration-call

if ischar(cmd) % 'query' mode, return the name of the preset
    apref = 'ANTs Anatomy of an antsRegistration call';
    return
end

apref.rigid.gradientstep = '0.1';
apref.rigid.metric = 'MI';
apref.rigid.metricparams = '1,32,Regular,0.25';
apref.rigid.convergence='[1000x500x250x100,1e-6,10]'; % Rigid convergence params
apref.rigid.shrinkfactors='8x4x2x1'; % Rigid shrink factors
apref.rigid.smoothingsigmas='3x2x1x0vox'; % Rigid Smoothing sigmas

apref.affine.gradientstep = '0.1';
apref.affine.metric = 'MI';
apref.affine.metricparams = '1,32,Regular,0.25';
apref.affine.convergence='[1000x500x250x100,1e-6,10]'; % Affine convergence params
apref.affine.shrinkfactors='8x4x2x1'; % Affine shrink factors
apref.affine.smoothingsigmas='3x2x1x0vox'; % Affine Smoothing sigmas

apref.syn.gradientstep = '0.1,3,0';
apref.syn.metric = 'CC';
apref.syn.metricparams = '1,4';
apref.syn.convergence='[100x70x50x20,1e-6,10]'; % SyN convergence params
apref.syn.shrinkfactors='8x4x2x1'; % SyN shrink factors
apref.syn.smoothingsigmas='3x2x1x0vox'; % SyN Smoothing sigmas

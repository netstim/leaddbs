function apref=ea_antspreset_SyN_fine(list)
if nargin % return name of preset
    apref='SyN Fine';
    return
end

prefs=ea_prefs;
switch prefs.machine.normsettings.ants_metric
    case 'Mutual Information'
        apref.metric='MI';
        apref.metricsuffix=',32,Regular,0.75';
    case 'ANTS Neighborhood Cross Correlation'
        apref.metric='CC';
        apref.metricsuffix=',4';
    case 'Global Correlation'
        apref.metric='GC';
        apref.metricsuffix=',15,Random,0.05';
end

switch prefs.machine.normsettings.ants_strategy
    case 'SyN';
        apref.antsmode='SyN';
        apref.antsmode_suffix='[0.1,3,0]';
    case 'BSplineSyN'
        apref.antsmode='BSplineSyN';
        apref.antsmode_suffix='[0.1,26,0,3]'; % as in example script in Tustison 2013
    case 'Exponential'
        apref.antsmode='Exponential';
        apref.antsmode_suffix='[0.1,26,0,3]'; % need to find proper values.
    case 'BSplineExponential'
        apref.antsmode='BSplineExponential';
        apref.antsmode_suffix='[0.1,26,0,3]'; % need to find proper values.
end

apref.small_large_dissociation=256; % Value from which imgsize is treated as "large" or "small", see below.
% Convergence
apref.convergence.rigid.large='[1000x500x250x100,1e-6,10]'; % Rigid convergence params for large volumes
apref.convergence.affine.large='[1000x500x250x100,1e-6,10]'; % Affine convergence params for large volumes
apref.convergence.syn.large='[100x100x70x50x20,1e-6,10]'; % SyN convergence params for large volumes

apref.convergence.rigid.small='[1000x500x250x100,1e-6,10]'; % Rigid convergence params for small volumes
apref.convergence.affine.small='[1000x500x250x100,1e-6,10]'; % Affine convergence params for small volumes
apref.convergence.syn.small='[100x70x50x20,1e-6,10]'; % SyN convergence params for small volumes

apref.convergence.scrf='[20x5,1e-6,10]'; % SyN subcortical focus stage convergence params

% Affine Convergence
apref.shrinkfactors.rigid.large='12x8x4x2'; % Rigid shrink factors for large volumes
apref.shrinkfactors.affine.large='12x8x4x2'; % Affine shrink factors for large volumes
apref.shrinkfactors.syn.large='10x6x4x2x1'; % SyN shrink factors for large volumes

apref.shrinkfactors.rigid.small='8x4x2x1'; % Rigid shrink factors for small volumes
apref.shrinkfactors.affine.small='8x4x2x1'; % Affine shrink factors for small volumes
apref.shrinkfactors.syn.small='8x4x2x1'; % SyN shrink factors for small volumes

apref.shrinkfactors.scrf='2x1'; % SyN subcortical focus stage shrink factors

% Smoothing Sigmas
apref.smoothingsigmas.rigid.large='4x3x2x1vox'; % Rigid Smoothing sigmas for large volumes
apref.smoothingsigmas.affine.large='4x3x2x1vox'; % Affine Smoothing sigmas for large volumes
apref.smoothingsigmas.syn.large='5x3x2x1x0vox'; % SyN Smoothing sigmas for large volumes

apref.smoothingsigmas.rigid.small='3x2x1x0vox'; % Rigid Smoothing sigmas for small volumes
apref.smoothingsigmas.affine.small='3x2x1x0vox'; % Affine Smoothing sigmas for small volumes
apref.smoothingsigmas.syn.small='3x2x1x0vox'; % SyN Smoothing sigmas for small volumes

apref.smoothingsigmas.scrf='1x0vox'; % SyN subcortical focus stage smoothing sigmas

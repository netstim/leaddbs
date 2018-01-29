function apref = ea_antspreset_legacy_default(cmd)

if ischar(cmd) % 'query' mode, return the name of the preset
    apref = 'Legacy: "Default"';
    return
end

normsettings = cmd; % cmd is the 'normsettings' struct
switch normsettings.ants_metric
    case 'Mutual Information'
        apref.metric='MI';
        apref.metricsuffix=',32,Regular,0.25';
    case 'ANTS Neighborhood Cross Correlation'
        apref.metric='CC';
        apref.metricsuffix=',4';
    case 'Global Correlation'
        apref.metric='GC';
        apref.metricsuffix=',15,Random,0.05';
end

switch normsettings.ants_strategy
    case 'SyN'
        apref.antsmode='SyN';
        apref.antsmode_suffix='[0.1,3,0]';
    case 'BSplineSyN'
        apref.antsmode='BSplineSyN';
        apref.antsmode_suffix='[0.1,26,0,3]'; % as in example script in Tustison 2013
end

% Convergence
apref.convergence.rigid='[500x250x120x0,1e-6,10]'; % Rigid convergence params
apref.convergence.affine='[500x250x120x0,1e-6,10]'; % Affine convergence params
apref.convergence.syn='[200x50x25x0,1e-6,10]'; % SyN convergence params
apref.convergence.scrf='[200x50x25x0,1e-6,10]'; % SyN subcortical focus stage convergence params

% Affine Convergence
apref.shrinkfactors.rigid='12x8x4x1'; % Rigid shrink factors
apref.shrinkfactors.affine='12x8x4x1'; % Affine shrink factors
apref.shrinkfactors.syn='6x4x2x1'; % SyN shrink factors
apref.shrinkfactors.scrf='4x2x2x1'; % SyN subcortical focus stage shrink factors

% Smoothing Sigmas
apref.smoothingsigmas.rigid='4x3x2x1vox'; % Rigid Smoothing sigmas
apref.smoothingsigmas.affine='4x3x2x1vox'; % Affine Smoothing sigmas
apref.smoothingsigmas.syn='3x2x1x1vox'; % SyN Smoothing sigmas
apref.smoothingsigmas.scrf='2x1x0x0vox'; % SyN subcortical focus stage smoothing sigmas

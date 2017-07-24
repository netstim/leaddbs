function apref=ea_antspreset_effective_lowvariance(list)
if nargin % return name of preset
    apref='Effective (Low Variance)';
    return
end

prefs=ea_prefs;
switch prefs.machine.normsettings.ants_metric
    case 'Mutual Information'
        apref.metric='MI';
        apref.metricsuffix=',32,Random,0.25';
    case 'ANTS Neighborhood Cross Correlation'
        apref.metric='CC';
        apref.metricsuffix=',4';
    case 'Global Correlation'
        apref.metric='GC';
        apref.metricsuffix=',15,Random,0.05';
end

switch prefs.machine.normsettings.ants_strategy
    case 'SyN'
        apref.antsmode='SyN';
        apref.antsmode_suffix='[0.3,4,3]';
    case 'BSplineSyN'
        apref.antsmode='BSplineSyN';
        apref.antsmode_suffix='[0.1,26,2,4]'; % as in example script in Tustison 2013
end

% Convergence
apref.convergence.rigid='[1000x500x250x0,1e-6,10]'; % Rigid convergence params
apref.convergence.affine='[1000x500x250x0,1e-6,10]'; % Affine convergence params
apref.convergence.syn='[1000x500x250x0,1e-6,7]'; % SyN convergence params
apref.convergence.scrf='[200x50x10x0,1e-6,7]'; % SyN subcortical focus stage convergence params

% Shrink Factors (please keep in mind these refer to fixed space which is
% typically 0.5 mm resolution in Lead-DBS whereas usually input images are
% resliced to 0.7 mm. Thus should consider not going down to shrink factor
% of 1 at all.
apref.shrinkfactors.rigid='12x8x4x1'; % Rigid shrink factors
apref.shrinkfactors.affine='8x4x2x1'; % Affine shrink factors
apref.shrinkfactors.syn='8x4x4x1'; % SyN shrink factors
apref.shrinkfactors.scrf='4x4x2x1'; % SyN subcortical focus stage shrink factors

% Smoothing Sigmas
apref.smoothingsigmas.rigid='4x3x2x1vox'; % Rigid Smoothing sigmas
apref.smoothingsigmas.affine='4x3x2x1vox'; % Affine Smoothing sigmas
apref.smoothingsigmas.syn='4x3x1x1vox'; % SyN Smoothing sigmas
apref.smoothingsigmas.scrf='2x2x1x1vox'; % SyN subcortical focus stage smoothing sigmas

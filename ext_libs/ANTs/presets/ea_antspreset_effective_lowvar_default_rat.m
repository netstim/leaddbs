function apref = ea_antspreset_effective_lowvar_default(cmd)

if ischar(cmd) % 'query' mode, return the name of the preset
    apref = 'Effective: Low Variance, Default';
    return
end

normsettings = cmd; % cmd is the 'normsettings' struct
switch normsettings.ants_metric
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

switch normsettings.ants_strategy
    case 'SyN'
        apref.antsmode='SyN';
        apref.antsmode_suffix='[0.3,4,3]';
        
        apref.convergence.syn='[1000x500x500x0,1e-7,7]'; % SyN convergence params
        apref.convergence.scrf='[200x50x10x0,1e-6,7]'; % SyN subcortical focus stage convergence params
        apref.smoothingsigmas.syn='4x3x1x1vox'; % SyN Smoothing sigmas
        apref.smoothingsigmas.scrf='2x2x1x1vox'; % SyN subcortical focus stage smoothing sigmas

        apref.shrinkfactors.syn='16x6x6x1'; % SyN shrink factors
        apref.shrinkfactors.scrf='6x6x4x1'; % SyN subcortical focus stage shrink factors
       
        % Shrink Factors (please keep in mind these refer to fixed space which is
        % typically 0.5 mm resolution in Lead-DBS whereas usually input images are
        % resliced to 0.7 mm. Thus should consider not going down to shrink factor
        % of 1 at all.
        apref.shrinkfactors.rigid='12x8x6x1'; % Rigid shrink factors
        apref.shrinkfactors.affine='12x8x6x1'; % Affine shrink factors

        
    case 'BSplineSyN'
        apref.antsmode='BSplineSyN';
        apref.antsmode_suffix='[0.1,26,0,3]';
        
        apref.convergence.syn='[100x70x50x20,1e-6,7]'; % SyN convergence params
        apref.convergence.scrf='[100x70x50x20,1e-6,7]'; % SyN subcortical focus stage convergence params
        
        apref.smoothingsigmas.syn='3x2x1x0vox'; % SyN Smoothing sigmas
        apref.smoothingsigmas.scrf='2x2x1x0vox'; % SyN subcortical focus stage smoothing sigmas
        
        apref.shrinkfactors.syn='6x4x2x1'; % SyN shrink factors
        apref.shrinkfactors.scrf='6x4x2x1'; % SyN subcortical focus stage shrink factors

        % Shrink Factors (please keep in mind these refer to fixed space which is
        % typically 0.5 mm resolution in Lead-DBS whereas usually input images are
        % resliced to 0.7 mm. Thus should consider not going down to shrink factor
        % of 1 at all.
        apref.shrinkfactors.rigid='12x8x4x1'; % Rigid shrink factors
        apref.shrinkfactors.affine='12x8x4x1'; % Affine shrink factors

        
end

% Convergence
apref.convergence.rigid='[1000x500x250x0,1e-7,10]'; % Rigid convergence params
apref.convergence.affine='[1000x500x250x0,1e-7,10]'; % Affine convergence params


% Smoothing Sigmas
apref.smoothingsigmas.rigid='5x4x3x1vox'; % Rigid Smoothing sigmas
apref.smoothingsigmas.affine='5x4x3x1vox'; % Affine Smoothing sigmas

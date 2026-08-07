function [Rpred, Rpred_p, Ihat] = ea_unified_perm_sweetspot_predictfromvals(obj, vals, patsel, corrType, efieldTIn)
% Replicates the sweetspotmapping E-Fields prediction logic (see
% obj.basepredictionon) that crossval()/ea_compute_unified_model.m use, but
% takes the voxelwise map(s) directly rather than recomputing them, so it
% can be reused cheaply across many permuted maps (ea_unified_perm_sweetspot_predict.m
% calls this once per permutation, up to Nperm times).
%
% vals    - {side} cell, one Nvoxels x 1 column vector per side (e.g. a
%           reprojected empirical or permuted training R-map).
% patsel  - patient indices to predict for.
% corrType - correlation type for the final Rpred (matches the training
%            permtest's corrType, not necessarily obj's own default).
% efieldTIn (optional) - precomputed obj.results.sweetspotmapping.efield{side}(patsel,:)'
%            per side, since that transpose/slice is independent of vals
%            and otherwise gets recomputed on every one of the Nperm calls
%            a caller like ea_unified_perm_sweetspot_predict.m makes.

if ~exist('efieldTIn', 'var') || isempty(efieldTIn)
    efieldTIn = {};
end

Ihat = nan(length(patsel), numel(vals));
for side = 1:numel(vals)
    v = vals{side};
    if numel(efieldTIn) >= side && ~isempty(efieldTIn{side})
        efieldT = efieldTIn{side};
    else
        efieldT = obj.results.sweetspotmapping.efield{side}(patsel,:)';
    end

    vmasked = v;
    if ~obj.posvisible
        vmasked(vmasked>0) = nan;
    end
    if ~obj.negvisible
        vmasked(vmasked<0) = nan;
    end

    switch lower(obj.basepredictionon)
        case 'profile of scores: spearman'
            Ihat(:,side) = atanh(ea_corr(vmasked, efieldT, 'spearman'));
        case 'profile of scores: pearson'
            Ihat(:,side) = atanh(ea_corr(vmasked, efieldT, 'pearson'));
        case 'profile of scores: bend'
            Ihat(:,side) = atanh(ea_corr(vmasked, efieldT, 'bend'));
        case 'mean of scores'
            Ihat(:,side) = ea_nanmean(vmasked.*efieldT,1);
        case 'sum of scores'
            Ihat(:,side) = ea_nansum(vmasked.*efieldT,1);
        case 'peak of scores'
            Ihat(:,side) = ea_discfibers_getpeak(v.*efieldT, obj.posvisible, obj.negvisible, 'peak');
        case 'peak 5% of scores'
            Ihat(:,side) = ea_discfibers_getpeak(v.*efieldT, obj.posvisible, obj.negvisible, 'peak5');
    end
end
Ihat = ea_nanmean(Ihat,2);
[Rpred, Rpred_p] = corr(obj.responsevar(patsel), Ihat, 'type', corrType, 'rows', 'pairwise');

function [vals,pvals,gvalOut,gpatselOut,thisvalsOut,nanidxOut] = ea_unified_perm_sweetspot_calcstats(obj,patsel,Iperm,skipsigthresh,gvalIn,gpatselIn,thisvalsIn,nanidxIn)
% Voxelwise Correlations statistic for sweetspotmapping, restricted to the
% same scope permutation testing is gated to (obj.statsettings.statfamily
% == 'Correlations', stimulationmodel in {'Electric Field','Sigmoid Field'}
% -- 'VTA' is not in those stat-test files' own .compatibility list, so it's
% excluded here too). A dedicated function rather than a change to
% ea_unifiedmapping_calcstats.m, which is shared by fiberfiltering/
% networkmapping and doesn't return p-values or support the caching
% arguments below -- and reusing it unmodified would actually be SLOWER for
% permtest's use case, since it has no caching mechanism and would redo the
% full (expensive) coverage-masking on every one of Nperm calls.
%
% Deliberately kept in numeric lockstep with ea_unifiedmapping_calcstats.m's
% 'sweetspotmapping' branch (same coverage-masking, same voxel-x-patient
% orientation), so the empirical (unpermuted) map this produces matches
% what obj.draw() already shows for the real map. Because this duplicates
% that logic rather than sharing it, ea_unified_perm_sweetspot.m runs
% ea_unified_perm_sweetspot_consistency_check.m once per permtest call to
% catch the two implementations drifting apart if the shared function is
% ever changed without mirroring the change here.
%
% patsel          - patient indices (obj.patientselection if not supplied).
% Iperm           - substituted in place of obj.responsevar when given (the
%                    permutation mechanism); real responsevar used otherwise.
% skipsigthresh   - if true, vals are returned unmasked regardless of
%                    obj.showsignificantonly (permutation testing wants the
%                    raw, unmasked map; thresholding happens post-hoc).
% gvalIn/gpatselIn, thisvalsIn/nanidxIn - caching triplet: each is
%                    independent of I/Iperm once the patient selection and
%                    coverage mask are fixed, so a caller running this
%                    Nperm times (permtest's parfor loop) can compute them
%                    once and skip re-materializing full-size arrays on
%                    every subsequent call.
%
% Returns per-side ({side}, no group dimension -- sweetspotmapping doesn't
% support groups in ea_unifiedmapping) cells: vals/pvals (Nvoxels x 1 each),
% plus the three cache outputs to feed back in on the next call.

if ~exist('Iperm','var') || isempty(Iperm)
    I = obj.responsevar;
else
    I = Iperm;
end

if ~exist('skipsigthresh','var') || isempty(skipsigthresh)
    skipsigthresh = false;
end
if ~exist('gvalIn','var'); gvalIn = {}; end
if ~exist('gpatselIn','var'); gpatselIn = []; end
if ~exist('thisvalsIn','var'); thisvalsIn = {}; end
if ~exist('nanidxIn','var'); nanidxIn = {}; end

if ~exist('patsel','var') || isempty(patsel)
    patsel = obj.patientselection;
end

% obj.results.sweetspotmapping.efield{side} is [nVAT x nVox]; transpose to
% [nVox x nVAT] to match ea_unifiedmapping_calcstats.m's convention (rows =
% voxels, columns = patients/VAT instances) before coverage-masking.
val = obj.results.sweetspotmapping.efield';
for i = 1:length(val)
    val{i} = val{i}';
end

if size(I,2)==1
    I = [I,I];
end
if obj.mirrorsides
    I = [I;I];
end
if ~isempty(obj.covars)
    for i = 1:length(obj.covars)
        if obj.mirrorsides
            covars_i = [obj.covars{i}; obj.covars{i}];
        else
            covars_i = obj.covars{i};
        end
        for side = 1:2
            if size(covars_i,2)==1
                I(:,side) = ea_resid(covars_i,I(:,side));
            elseif size(covars_i,2)==2
                I(:,side) = ea_resid(covars_i(:,side),I(:,side));
            end
        end
    end
end

corrType = obj.statsettings.stattest;
if strcmpi(corrType, 'Percentage Bend')
    corrType = 'bend'; % ea_corr's own naming for this test
end

if ~isempty(gpatselIn)
    gpatsel = gpatselIn; % reuse precomputed selection (independent of I/Iperm)
else
    gpatsel = patsel;
    if obj.mirrorsides
        gpatsel = [gpatsel, gpatsel+length(obj.allpatients)];
    end
end
gpatselOut = gpatsel;

gval = val;
vals = cell(1,numel(gval));
pvals = cell(1,numel(gval));
gvalOut = cell(1,numel(gval));
thisvalsOut = cell(1,numel(gval));
nanidxOut = cell(1,numel(gval));

for side = 1:numel(gval)
    if numel(gvalIn)>=side && ~isempty(gvalIn{side})
        gval{side} = gvalIn{side}; % reuse precomputed coverage-masked efield data
    else
        switch obj.statsettings.stimulationmodel
            case 'Electric Field'
                if isfield(obj.statsettings,'nanthreshold') && ~isempty(obj.statsettings.nanthreshold)
                    gval{side}(gval{side} <= obj.statsettings.nanthreshold) = nan;
                end
                Nmap = ea_nansum((gval{side}(:,gpatsel) > obj.statsettings.efieldthreshold_spot), 2);
                gval{side}(Nmap < round((obj.statsettings.connthreshold/100)*length(gpatsel)), gpatsel) = nan;
            case 'Sigmoid Field'
                gval{side}(:,gpatsel) = ea_SigmoidFromEfield(gval{side}(:,gpatsel));
                Nmap = ea_nansum(gval{side}(:,gpatsel) > obj.statsettings.efieldthreshold_spot, 2);
                gval{side}(Nmap < round((obj.statsettings.connthreshold/100)*length(gpatsel)), gpatsel) = nan;
            otherwise
                error('ea_unified_perm_sweetspot_calcstats:stimulationmodel', ...
                    'stimulationmodel ''%s'' is not compatible with Correlations (only ''Electric Field''/''Sigmoid Field'' are).', ...
                    obj.statsettings.stimulationmodel);
        end
    end
    gvalOut{side} = gval{side};

    if numel(thisvalsIn)>=side && ~isempty(thisvalsIn{side})
        thisvals = thisvalsIn{side}; % reuse precomputed, patient-sliced & NaN-filtered input
        nanidx = nanidxIn{side};
    else
        thisvals = gval{side}(:,gpatsel)'; % -> [nSelPatients x nVox], ea_corr's expected orientation
        % Matches ea_unifiedmapping_calcstats.m's voxel-inclusion criterion
        % exactly: nonempty = sum(gval{side}(:,gpatsel),2,'omitnan')>0. Not
        % just "not all-NaN" (isnan(ea_nansum(...)), the old TSProjectSweetspot
        % idiom this was originally ported from) -- a voxel whose covered
        % values are all exactly 0 (e.g. every contributing patient has zero
        % field magnitude there) has a defined, non-NaN sum, but that sum is
        % 0, which the live code correctly treats as "no real signal" and
        % excludes. Using the old, looser criterion let such a voxel through
        % here while ea_unifiedmapping_calcstats.m correctly left it NaN --
        % caught by ea_unified_perm_sweetspot_consistency_check.m.
        nanidx = ~(sum(thisvals, 1, 'omitnan') > 0);
        thisvals = thisvals(:, ~nanidx);
    end
    thisvalsOut{side} = thisvals;
    nanidxOut{side} = nanidx;

    [R,p] = ea_corr(thisvals, I(gpatsel,side), corrType);

    vals{side} = nan(numel(nanidx),1);
    pvals{side} = nan(numel(nanidx),1);
    pvals{side}(~nanidx) = p;

    if obj.showsignificantonly && ~skipsigthresh
        R = ea_unified_perm_corrsignan(R, p, obj);
    end
    vals{side}(~nanidx) = R;
end
end

function vals = ea_unified_perm_corrsignan(vals, ps, obj)
% Same FDR/Bonferroni multiple-comparison masking as the shared explorer
% path (ea_corrsignan in ea_unifiedmapping_calcstats.m), duplicated here
% rather than shared since that function isn't exported.
nnanidx = ~isnan(vals);
numtests = sum(nnanidx);

switch lower(obj.multcompstrategy)
    case 'fdr'
        pnnan = ps(nnanidx);
        [psort,idx] = sort(pnnan); %#ok<ASGLU>
        pranks = zeros(1,length(psort));
        for rank = 1:length(pranks)
            pranks(idx(rank)) = rank;
        end
        pnnan = pnnan.*numtests;
        if ~isequal(size(pranks),size(pnnan))
            pranks = pranks';
        end
        pnnan = pnnan./pranks;
        ps(nnanidx) = pnnan;
    case 'bonferroni'
        ps(nnanidx) = ps(nnanidx).*numtests;
end
ps(~nnanidx) = 1;
vals(ps>obj.alphalevel) = nan;
end

function tfcestat = ea_unified_perm_palm_tfce(X, datatype, geom, E, H, conn, deltah, area)
% CITE IF USED: this is a port of PALM (Permutation Analysis of Linear
% Models)'s palm_tfce.m -- https://github.com/andersonwinkler/PALM -- by
% Anderson M. Winkler (FMRIB / University of Oxford, Sep/2013). If you
% report results computed with this function, cite:
%   Winkler AM, Ridgway GR, Douaud G, Nichols TE, Smith SM (2016). Faster
%   permutation inference in brain imaging. NeuroImage 141:502-516. (PALM)
%   Smith SM, Nichols TE (2009). Threshold-free cluster enhancement:
%   addressing problems of smoothing, threshold dependence and
%   localisation in cluster inference. NeuroImage 44:83-98. (the TFCE
%   method itself, which palm_tfce.m implements)
%
% Computes the TFCE (Threshold-Free Cluster Enhancement) statistic for a
% volume, or for vertexwise/facewise surface data.
%
% Adapted from palm_tfce.m for this codebase:
% - PALM's opts/plm/y struct plumbing (whole-pipeline machinery for an
%   entire permutation analysis run, most of it unrelated to TFCE itself)
%   is replaced with plain arguments -- everything actually used by
%   palm_tfce.m's TFCE math is exposed here directly, nothing else.
% - Both of palm_tfce.m's branches (volume via bwconncomp; vertex/facewise
%   via palm_dpxlabel, ported alongside as ea_unified_perm_palm_dpxlabel.m) are
%   kept, even though only the volume branch has a caller in this codebase
%   today (ea_unified_perm_sweetspot_tfce.m) -- this is deliberate, so the
%   vertex/facewise path is already here, ported and untouched, whenever a
%   surface-based tool (e.g. network mapping) needs TFCE later.
% - PALM's separate mask+injection step (mask = S.data; D(mask) = X) is
%   replaced with NaN-based masking: this codebase's statistic maps
%   (Remp/Rperm) are already full-grid-length vectors with NaN marking
%   uncovered positions, and NaN >= h is always false in MATLAB -- so a NaN
%   position can never join a cluster at any threshold, reproducing PALM's
%   masking behavior without a separate injection step. NaN is restored at
%   those positions in the output (their accumulated value would otherwise
%   be a misleading 0, not "not applicable").
%
% Usage:
%   tfcestat = ea_unified_perm_palm_tfce(X, datatype, geom, E, H, conn, deltah, area)
%
% X        - statistic map, Nx1 or 1xN, NaN = excluded/uncovered position.
% datatype - 'volume' | 'vertex' | 'facewise'.
% geom     - for 'volume': the 3D grid size, e.g. space{side}.dim(1:3).
%            for 'vertex'/'facewise': a sparse adjacency matrix (see PALM's
%            palm_adjacency.m for how one would be built from a surface
%            mesh -- not ported here, out of scope until a caller with
%            actual surface data exists).
% E, H     - TFCE extent/height exponents. Default 0.5 / 2 (PALM/FSL default).
% conn     - volume connectivity for bwconncomp (6, 18, or 26). Default 6
%            (face-adjacency, PALM's default). Unused for vertex/facewise
%            (connectivity there comes from the adjacency matrix itself).
% deltah   - threshold step size. 0 (default) means "auto": dh = max(X)/100,
%            matching PALM's opts.tfce.deltah==0 sentinel/'auto' convention.
% area     - vertex/facewise only: per-vertex/per-face area weights, same
%            length as X. Default: all-ones (equivalent to plain vertex/
%            face counting) -- matches PALM's plm.Yarea when real surface
%            area weights aren't available. Unused for 'volume'.
%
% Returns tfcestat, same shape as X, NaN preserved at excluded positions.

if ~exist('E','var') || isempty(E)
    E = 0.5;
end
if ~exist('H','var') || isempty(H)
    H = 2;
end
if ~exist('conn','var') || isempty(conn)
    conn = 6;
end
if ~exist('deltah','var') || isempty(deltah)
    deltah = 0; % 0 = 'auto', matching PALM's opts.tfce.deltah==0 sentinel
end

Xshape = size(X);
X = X(:);
nanmask = isnan(X);
maxX = max(X(~nanmask));

switch lower(datatype)
    case 'volume'
        dim = geom;
        dim = dim(:)';
        D = nan(dim);
        D(:) = X;

        tfcestat3d = zeros(dim);
        if deltah == 0
            dh = maxX/100;
            for h = dh:dh:maxX
                CC = bwconncomp(D>=h, conn);
                integ = cellfun(@numel,CC.PixelIdxList).^E * h^H;
                for c = 1:CC.NumObjects
                    tfcestat3d(CC.PixelIdxList{c}) = ...
                        tfcestat3d(CC.PixelIdxList{c}) + integ(c);
                end
            end
        else
            dh = deltah;
            h  = dh;
            CC = bwconncomp(D>=h, conn);
            while CC.NumObjects
                integ = cellfun(@numel,CC.PixelIdxList).^E * h^H;
                for c = 1:CC.NumObjects
                    tfcestat3d(CC.PixelIdxList{c}) = ...
                        tfcestat3d(CC.PixelIdxList{c}) + integ(c);
                end
                h  = h + dh;
                CC = bwconncomp(D>=h, conn);
            end
        end
        tfcestat = tfcestat3d(:) * dh;

    case {'vertex','facewise'}
        adj = geom;
        if ~exist('area','var') || isempty(area)
            area = ones(size(X));
        else
            area = area(:);
        end
        D = X; % NaN positions never satisfy D>=h -- same masking-for-free logic as the volume branch above
        tfcestat = zeros(size(D));

        if deltah == 0
            dh = maxX/100;
            for h = dh:dh:maxX
                dpxl = ea_unified_perm_palm_dpxlabel(D>=h, adj);
                U = unique(dpxl(dpxl>0))';
                for u = 1:numel(U)
                    idx = dpxl==U(u);
                    tfcestat(idx) = tfcestat(idx) + sum(area(idx)).^E * h^H;
                end
            end
        else
            dh = deltah;
            h  = dh;
            dpxl = ea_unified_perm_palm_dpxlabel(D>=h, adj);
            U = unique(dpxl(dpxl>0))';
            while numel(U)
                for u = 1:numel(U)
                    idx = dpxl==U(u);
                    tfcestat(idx) = tfcestat(idx) + sum(area(idx)).^E * h^H;
                end
                h    = h + dh;
                dpxl = ea_unified_perm_palm_dpxlabel(D>=h, adj);
                U    = unique(dpxl(dpxl>0))';
            end
        end
        tfcestat = tfcestat(:) * dh;

    otherwise
        error('ea_unified_perm_palm_tfce:datatype', 'datatype must be ''volume'', ''vertex'', or ''facewise''.');
end

tfcestat(nanmask) = nan;
tfcestat = reshape(tfcestat, Xshape);

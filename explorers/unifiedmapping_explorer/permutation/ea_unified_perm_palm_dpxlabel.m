function dpvl = ea_unified_perm_palm_dpxlabel(dpx,adj)
% Ported, near-verbatim, from PALM (Permutation Analysis of Linear Models)'s
% palm_dpxlabel.m -- https://github.com/andersonwinkler/PALM -- used here as
% the connected-component labeling engine for the vertex/facewise branch of
% ea_unified_perm_palm_tfce.m (the volume branch uses MATLAB's bwconncomp
% instead; this function is for surface-mesh data via a sparse adjacency
% matrix). If you use TFCE with vertex/facewise data via this file, cite
% PALM (Winkler et al. 2016, NeuroImage 141:502-516) -- see
% ea_unified_perm_palm_tfce.m's header for the full citation note.
%
% Label a DPV file (vertexwise data).
%
% Usage:
% dpxl = ea_unified_perm_palm_dpxlabel(dpx,adj)
%
% - dpx  : Data per vertex or per face to be labelled. It should be
%          a logical vector, and the values marked as "true" are
%          the ones labelled.
% - adj  : Sparse adjacency matrix (see palm_adjacency for details).
% - dpvl : Labelled data.
%
% References:
% * Dulmage AL, Mendelsohn NS. Coverings of bipartite graphs.
%   Journal Canadien de Mathematiques. 1958;10(0):517-534.
%   http://dx.doi.org/10.4153/CJM-1958-052-0
% * Pothen A, Fan CJ. Computing the block triangular form of a
%   sparse matrix. ACM Transactions on Mathematical Software.
%   1990;16(4):303-324. http://dx.doi.org/10.1145/98267.98287
%
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Feb/2012 (first version)
% Jul/2015 (this version)
% http://brainder.org

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PALM -- Permutation Analysis of Linear Models
% Copyright (C) 2015 Anderson M. Winkler
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Remove from the adjacency matrix the entries that won't be labelled.
adj(~ dpx,:) = [];
adj(:,~ dpx) = [];

% Partition the adjacency with the Dulmage-Mendelsohn algorithm.
% See the MATLAB help for details.
[p,~,r]   = dmperm(adj);
C         = zeros(size(adj,1),1);
for j = 1:numel(r)-1
    C(p(r(j):r(j+1)-1)) = j;
end
dpvl      = zeros(size(dpx));
dpvl(dpx) = C;

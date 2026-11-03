function O = ea_permtest_cv_scores(cvb, PermIdx)
% Score matrix for ea_permtest_cv_predict.m: columns 1:Nvoters = the
% models' scores, Nvoters+1 = combined target (NaN if none),
% Nvoters+1+v = weight variable of model v (NaN if its weight is fixed).
%
% cvb     - from ea_permtest_cv_adapter_disctract.m.
% PermIdx - optional Npatients x B x (Nvoters + 1 + Nvoters): for every
%           version and column, the patient whose value each patient
%           receives (see ea_permtest_derive_shuffles.m). Without it, the
%           real scores are returned (B = 1).
%
% O - Npatients x (2 * Nvoters + 1) x B.

V = numel(cvb.voterLabels);
N = cvb.Npatients;
S = nan(N, 2 * V + 1);
S(:, 1:V) = cvb.outcomes;
if ~isempty(cvb.extra.target)
    S(:, V + 1) = cvb.extra.target;
end
S(:, V + 1 + (1:V)) = cvb.extra.weightvars;
if ~exist('PermIdx', 'var') || isempty(PermIdx)
    O = S;
    return
end
B = size(PermIdx, 2);
O = nan(N, 2 * V + 1, B);
for k = 1:size(S, 2)
    col = S(:, k);
    O(:, k, :) = reshape(col(PermIdx(:, :, k)), N, 1, B);
end

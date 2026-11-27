function [ranks, ties, adj] = ea_permtest_midranks(X)
% Ranks along each row of X, tied values get the mean of their ranks
% ("midranks", as tiedrank / corr's Spearman use).
%
% X     - M x N, no NaN.
% ranks - M x N (multiples of 0.5).
% ties  - M x 1 logical, true when a row contains tied values.
% adj   - M x 1 tie adjustment as corr's Spearman uses it:
%         sum over groups of tied values of (t^3 - t) / 2.

[M, N] = size(X);
[S, order] = sort(X, 2);
isStart = [true(M, 1), diff(S, 1, 2) ~= 0];
groupOf = cumsum(isStart, 2) + (0:M-1)' * N; % tie group, unique over rows
pos = repmat(1:N, M, 1);
sumPos = accumarray(groupOf(:), pos(:), [M * N, 1]);
count = accumarray(groupOf(:), 1, [M * N, 1]);
ranks = zeros(M, N);
ranks(sub2ind([M, N], repmat((1:M)', 1, N), order)) = sumPos(groupOf) ./ count(groupOf);
ties = ~all(isStart, 2);
if nargout > 2
    rowOfGroup = floor(((1:M * N)' - 1) / N) + 1;
    adj = accumarray(rowOfGroup, (count.^3 - count) / 2, [M, 1]);
end

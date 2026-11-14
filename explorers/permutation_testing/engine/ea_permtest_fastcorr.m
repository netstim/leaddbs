function [R, P] = ea_permtest_fastcorr(prep, Y)
% R and p of the explorer's correlation test for many outcomes at once,
% using the fiber ranks / values prepared by ea_permtest_fastcorr_prepare.m
% (see there for how corr's results are reproduced).
%
% prep - from ea_permtest_fastcorr_prepare (prep.ok must be true).
% Y    - N x B outcomes (e.g. B shuffles), rows as passed to the test; NaN
%        exactly where the real outcome is NaN.
%
% R, P - Nfibers x B.

if ~prep.ok
    ea_error(['Fast correlation path not available: ', prep.reason]);
end
if any(~isnan(Y(~prep.rows, :)), 'all') || any(isnan(Y(prep.rows, :)), 'all')
    ea_error('The shuffled outcomes have missing scores in other places than the real outcome.');
end
Y = Y(prep.rows, :);
n = prep.n;
switch prep.type
    case 'pearson'
        Y = Y - mean(Y, 1);
        Y = Y ./ vecnorm(Y, 2, 1);
        R = prep.Zx * Y;
        R(abs(R) > 1) = sign(R(abs(R) > 1));
        t = R .* sqrt((n - 2) ./ (1 - R.^2));
        P = 2 * tcdf(-abs(t), n - 2);
    case 'spearman'
        [Yr, ~, yadj] = ea_permtest_midranks(Y');
        Yr = Yr';
        yadj = yadj';
        % sums of squared rank differences, exact (multiples of 0.25)
        D = prep.sx2 + sum(Yr.^2, 1) - 2 * (prep.Xr * Yr);
        n3const = (n+1)*n*(n-1) ./ 3;
        R = nan(size(D));
        P = nan(size(D));
        tb = prep.tieBranch;
        if any(tb)
            % corr's ties branch (corrSpearman>vectorizedcorr / pvalSpearman)
            Dt = D(tb, :);
            xadj = prep.xadj(tb);
            meanD = (n3const - (xadj + yadj)./3) ./ 2;
            stdD = sqrt((n3const./2 - xadj./3) .* (n3const./2 - yadj./3) ./ (n-1));
            n3const2 = (n+1)*n*(n-1)/2;
            stdD((xadj == n3const2) | (yadj == n3const2)) = 0;
            coef = (meanD - Dt) ./ (sqrt(n-1)*stdD);
            Dt(stdD == 0) = NaN;
            r = (meanD - Dt) ./ (sqrt(n-1)*stdD);
            t = Inf*sign(r);
            ok = (abs(r) < 1);
            t(ok) = r(ok) .* sqrt((n-2)./(1-r(ok).^2));
            P(tb, :) = 2*tcdf(-abs(t), n-2);
            k = find(abs(coef) > 1);
            coef(k) = coef(k)./abs(coef(k)); % preserves NaNs
            R(tb, :) = coef;
        end
        if any(~tb)
            % corr's no-ties branch: Pearson of ranks, p from AS 89
            Y0 = Yr - (n + 1) / 2;
            dy = vecnorm(Y0, 2, 1);
            coef = (prep.X0(~tb, :) * Y0) ./ prep.dx(~tb) ./ dy;
            k = find(abs(coef) > 1);
            coef(k) = coef(k)./abs(coef(k));
            meanD = n3const./2;
            stdD = n3const ./ (2*sqrt(n-1));
            Dn = meanD - round(coef*(sqrt(n-1)*stdD));
            n3c = (n^3 - n)/3;
            p = 2 * local_as89(max(Dn, n3c - Dn), n, n3c);
            p(p > 1) = 1;
            R(~tb, :) = coef;
            P(~tb, :) = p;
        end
end


function p = local_as89(D, n, n3const)
% Upper tail probability of Spearman's D without ties: Algorithm AS 89
% (Best & Roberts 1975, Appl Stat 24:377-379), Edgeworth series
% approximation for n >= 10, as used by corr.
c = [.2274 .2531 .1745 .0758 .1033 .3932 .0879 .0151 .0072 .0831 .0131 .00046];
x = (2 * (D - 1) ./ n3const - 1) * sqrt(n - 1);
y = x .* x;
u = x .* (c(1) + (c(2) + c(3)/n)/n + ...
    y .* (-c(4) + (c(5) + c(6)/n)/n - ...
    y .* (c(7) + c(8)/n - y .* (c(9) - c(10)/n + ...
    y .* (c(11) - c(12) * y)/n))/n))/n;
p = u ./ exp(.5 * y) + 0.5 * erfc(x ./ sqrt(2));
p(p > 1) = 1;
p(p < 0) = 0;

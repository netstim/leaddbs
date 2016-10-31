function [Iab,Pab,Pa,Pb] = mutualInformation(a,b,normalize,nbins)
%
% function [Iab,Pab,Pa,Pb] = mutualInformation(a,b,[nbins])
%
% Computes the mutual information between two vectors. Uses hist2 to
% compute the joint histogram (which ignores the tails of the two marginal
% distributions. Mutual information is:
%    I(a,b) = H(a) + H(b) - H(a,b)
% where
%    H(a) = - sum P(a) log[P(a)]
%
% Normalized mutual information is:
%    [H(a) + H(b)] / H(a,b)
%
% Default nbins: sqrt(length(v1)/10)
% Default normalize: 0
% 
% djh, 3/2005

% Default normalize
if ~exist('normalize','var')
    normalize = 0;
end

% Default nbins
if ~exist('nbins','var')
    nbins = round(sqrt(length(a)/10));
end

% Joint histogram
abHist = hist2(a,b,nbins);

% Marginal histograms
aHist = sum(abHist,1);
bHist = sum(abHist,2);

% Probabilities
N = sum(aHist);
Pa = aHist/N;
Pb = bHist/N;
Pab = abHist/N;

% Disable divide by 0 and log of 0 warnings
warning('off');
Ha = (Pa .* log(Pa));
id = isfinite(Ha);
Ha = - sum(Ha(id));

Hb = (Pb .* log(Pb));
id = isfinite(Hb);
Hb = - sum(Hb(id));

Hab = (Pab .* log(Pab));
id = isfinite(Hab);
Hab = - sum(Hab(id));
warning('on');

if normalize
    Iab = (Ha + Hb) / Hab;
else
    Iab = Ha + Hb - Hab;
end

return

% Test
a = randn(1e5,1);
b = randn(1e5,1);
norm = 1;
n=100;
mutualInformation(a,a,norm,n)
mutualInformation(b,b,norm,n)
mutualInformation(a,b,norm,n)

[Iab,Pab,Pa,Pb] = mutualInformation(a,b);
figure(1)
showIm(Pab)
figure(2)
showIm(Pb*Pa)

[djh,djhHdr] = tfiReadAnalyze('/Users/david/data/anatomy/djh/dh_mprage1');
[inplane,inplaneHdr] = tfiReadAnalyze('/Users/david/data/MLR4-data/djh-retino080603/Inplane/Anatomy/dh080603retino+05+InPlane_T1_24sl.img');


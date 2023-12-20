function [r] = ea_circcorr(y,X,exponent)


[~,score] = pca(X);
if ~exist('exponent','var')
    exponent=1.5;
end

res=ea_nanzscore(ea_resid(score(:,1:5),y));
weights=abs(res).^exponent;

cnt=1; chunks=40; % a chunk of 40 seems to be a good speed tradeoff for this operation.
%ea_dispercent(0,'Iterating voxels');
nsz=size(X,2);
for chunk=1:chunks:nsz
    if (cnt+chunks-1)>nsz
        from=cnt; to=nsz;
    else
        from=cnt; to=cnt+chunks-1;
    end

    R=weightedcorrs([y,X(:,from:to)],weights);
    r(from:to)=R(2:end,1);
    cnt=cnt+chunks;
    %ea_dispercent(chunk/nsz);
end
%ea_dispercent(1,'end');




function R = weightedcorrs(Y, w)
%
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%                                                                                               %
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 23 July 2008                                                      %
%            Updated: 6 June 2012                                                 %
%                                                                                               %
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%
% % ======================================================================
%

% % Check input
% ctrl = isvector(w) & isreal(w) & ~any(isnan(w)) & ~any(isinf(w)) & all(w > 0);
% if ctrl
%   w = w(:) / sum(w);                                                          % w is column vector
% else
%   error('Check w: it needs be a vector of real positive numbers with no infinite or nan values!')
% end
% ctrl = isreal(Y) & ~any(isnan(Y)) & ~any(isinf(Y)) & (size(size(Y), 2) == 2);
% if ~ctrl
%   error('Check Y: it needs be a 2D matrix of real numbers with no infinite or nan values!')
% end
% ctrl = length(w) == size(Y, 1);
% if ~ctrl
%   error('size(Y, 1) has to be equal to length(w)!')
% end

% clean from nan/inf:
idx=all(isfinite(Y),2);
Y=Y(idx,:);
w=w(idx);

[T, N] = size(Y);                                                             % T: number of observations; N: number of variables
temp = Y - repmat(w' * Y, T, 1);                                              % Remove mean (which is, also, weighted)
temp = temp' * (temp .* repmat(w, 1, N));                                     % Covariance Matrix (which is weighted)
temp = 0.5 * (temp + temp');                                                  % Must be exactly symmetric
R = diag(temp);                                                               % Variances
R = temp ./ sqrt(R * R');                                                     % Matrix of Weighted Correlation Coefficients

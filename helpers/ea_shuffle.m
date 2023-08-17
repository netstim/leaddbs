function [xPerm, permInd] = ea_shuffle(X, N, Sel, rngseed)
% Generate N random permutations of input vector X
% Can also shuffle only the [Sel]ected elements in X
% Return the shuffled X and permutation indices as row vectors

if ~exist('Sel', 'var')
    Sel = 1:length(X);
end

if ~exist('rngseed', 'var')
    rngseed = 'default';
end

rng(rngseed);

if isvector(X)
    if iscolumn(X)
        X = X';
    end
    
    xPerm = repmat(X, N, 1);
    permInd = repmat(1:length(X), N, 1);
    
    for i=1:N
        permInd(i, Sel) = Sel(randperm(length(Sel)));
        xPerm(i, :) = xPerm(i, permInd(i, :));
    end

elseif ismatrix(X)
    xPerm(1,:,:) = X; 
    xPerm = repmat(xPerm,N,1,1); 
    permInd = repmat(1:length(X), N, 1);

    for i=1:N
        permInd(i, Sel) = Sel(randperm(length(Sel)));
        xPerm(i, :, :) = xPerm(i, permInd(i, :), :); % permute lines, keep relationship between the different variables within a line
    end
end
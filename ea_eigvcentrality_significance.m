function [ixes]=ea_eigvcentrality_significance(XYZV)




% clear NAN vars:
XYZV(isnan(XYZV(:,4)),:)=[];


% setup randomization:
ea_dispercent(0,'Permutation based statistics on eigenvector centrality');
maxiter=10000;
for iter=1:maxiter
    dat=[XYZV,XYZV(:,4),XYZV(:,4)];
    if iter>1
        dat(:,4:end)=dat(randperm(length(dat)),4:end);
    end
    dat=zscore(dat);
    % calculate proximity matrix P:
    P=squareform(pdist(dat)); % distance matrix
    
    P=1./P; % convert to proximity matrix
    % calculate centrality on matrix P (for values that lie close to each
    % other both in spatial and value domains).
    P(logical(eye(size(P,1))))=1;
    v = eigenvector_centrality_und(P);    
    csize(iter,:)=v;
    ea_dispercent(iter/maxiter);
end
ea_dispercent(1,'end');

nullmodel=csize(2,:);
nullmodel=nullmodel(:);

% determine cutoff for p-value.
nullmodel=sort(nullmodel);

mincsizeval=nullmodel(ceil((1-0.1)*numel(nullmodel)));

ixes=csize(1,:)>mincsizeval;

disp([num2str(sum(ixes)),' values identified above threshold at p>0.05.']);
%end



% GM=fitgmdist(XYZV,1);
% P=posterior(GM,XYZV);
% idx=cluster(GM,XYZV);

%mdl = fitglm(XYZV(:,1:end-1),XYZV(:,end))


function   v = eigenvector_centrality_und(CIJ)
%EIGENVECTOR_CENTRALITY_UND      Spectral measure of centrality
%
%   v = eigenvector_centrality_und(CIJ)
%
%   Eigenector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector 
%   corresponding to the largest eigenvalue of the adjacency matrix.
%
%   Inputs:     CIJ,        binary/weighted undirected adjacency matrix.
%
%   Outputs:      v,        eigenvector associated with the largest
%                           eigenvalue of the adjacency matrix CIJ.
%
%   Reference: Newman, MEJ (2002). The mathematics of networks.
%
%   Xi-Nian Zuo, Chinese Academy of Sciences, 2010
%   Rick Betzel, Indiana University, 2012

n = length(CIJ) ;
if n < 1000
    [V,~] = eig(CIJ) ;
    ec = abs(V(:,n)) ;
else
    [V, ~] = eigs(sparse(CIJ)) ;
    ec = abs(V(:,1)) ;
end
v = reshape(ec, length(ec), 1);
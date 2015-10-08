function [ixes]=ea_eigvcentrality_significance(XYZV)




% clear NAN vars:
XYZV(isnan(XYZV(:,4)),:)=[];


% setup randomization:
ea_dispercent(0,'Permutation based statistics on strength centrality');
maxiter=5001;
csize=zeros(maxiter,size(XYZV,1));
for iter=1:maxiter
    dat=XYZV;
    if iter>1
        dat(:,4)=dat(randperm(length(dat)),4);
    end
    
Peuc=squareform(pdist(dat(:,1:3))); % distance matrix
Pval=squareform(pdist(dat(:,4))); % distance matrix

% exp
Peuc(logical(eye(size(Peuc,1))))=1./exp(Peuc(logical(eye(size(Peuc,1)))));
Pval(logical(eye(size(Peuc,1))))=1./exp(Pval(logical(eye(size(Peuc,1)))));
[X,Y]=meshgrid(dat(:,4),dat(:,4));
vals=mean(cat(3,X,Y),3);
Pval=vals.^Pval;

P=Peuc.^Pval;

P(logical(eye(size(P,1))))=0;
P=sum(P);

    csize(iter,:)=P;
    m(iter)=max(P);
    
    ea_dispercent(iter/maxiter);
end
ea_dispercent(1,'end');



%% maximum testing:

nullmodel=csize(2:end,:);
nullmodel=max(nullmodel,[],2);
nullmodel=sort(nullmodel,'descend');
mincsizeval=nullmodel(floor((0.05)*numel(nullmodel)));
realvals=csize(1,:);
ixes=realvals>mincsizeval;
disp([num2str(sum(ixes)),' values identified above threshold at p>0.05 in conservative testing.']);



%% pairwise testing:



nullmodel=csize(2:end,:);
nullmodel=sort(nullmodel,2,'descend'); % to be able to compare every 1st, 2nd, 3rd, etc. value with each other

nullmodel=sort(nullmodel,1,'descend'); % to be able to set up a threshold..

[realvals,ids]=sort(csize(1,:),'descend');
signodecnt=0;
keyboard
for node=1:size(nullmodel,2)
    thisnodesnullmodel=nullmodel(:,node); % all 1st, 2nd, 3rd, etc. values..
    mincsizeval=thisnodesnullmodel(floor((0.05)*numel(thisnodesnullmodel)));
    ps=realvals(node)>thisnodesnullmodel;
    p(node)=sum(~ps)/length(ps);
    if realvals(node)>mincsizeval
       signodecnt=signodecnt+1;
    else
        break
    end
end

FDR=mafdr(


ids=ids(1:signodecnt);
ixes=zeros(1,size(XYZV,1));
ixes(ids)=1;
ixes=logical(ixes);
disp([num2str(signodecnt),' values identified above threshold at p>0.05.']);
%keyboard
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

function [str] = strengths_und(CIJ)
%STRENGTHS_UND        Strength
%
%   str = strengths_und(CIJ);
%
%   Node strength is the sum of weights of links connected to the node.
%
%   Input:      CIJ,    undirected weighted connection matrix
%
%   Output:     str,    node strength
%
%
%   Olaf Sporns, Indiana University, 2002/2006/2008

% compute strengths
str = nansum(CIJ);        % strength



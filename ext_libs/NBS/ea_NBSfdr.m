function [n_cnt,con_mat,pval]=ea_NBSfdr(varargin)
%NBSfdr Simes procedure to determine edges surviving a false discovery
%rate (FDR) threshold. A p-value for each edge is determined using 
%permuted data that has been supplied as part of the structure STATS. 
%
%   [N_CNT,CON_MAT,PVAL]=NBSfdr(STATS) operates on the network and
%   associated statistical data in the structure STATS. 
%
%   [...]=NBSfdr(STATS,H) attempts to write out progress to uiwaitbar 
%   with handle H.  
%
%   [...]=NBSfdr(STATS,H,GLM) GLM is mandatory if the permuted data has 
%   not been precomputed, which is the case when STATS.test_stat is empty. 
%   In this situation, NBSglm is repeatedly called to compute a test
%   statistic for each permutation. Slower than precomputation, but saves 
%   memory. 
%
%   A STATS structure contains the following fields:
%       STATS.thresh:     Primary test statistic threshold [redundant here,
%                         only used for NBS]    
%       STATS.alpha:      False discovery rate threshold (user specified)
%       STATS.N:          Number of nodes in network
%       STATS.test_stat:  K+1 x J array of test statistics. The first 
%                         row is the oberved test statistic. The remaining 
%                         rows are samples of the test statistic under the 
%                         null hypothesis. Each column corresponds to a 
%                         seperate edge. K is number of permutations. J is
%                         the number of edges. The test statistics can be
%                         computed with NBSglm. Columns are mapped to
%                         edges such that column i=1:J corresponds to the
%                         edge with index ind_upper(i), where ind_upper are
%                         the indexes of the upper trianguler elements. 
%                         ind_upper = find(triu(ones(N,N),1)); 
%       STATS.size        {'extent'} | 'intensity' 
%                         Measure used to assess size of a network 
%                         component [redundant here, only used for NBS] 
%                          
%   Outputs:
%       N_CNT:            Set to 1 if at least one edge survives the false
%                         discovery rate threshold, otherwise set to 0
%       CON_MAT:          Cell array containing an N x N 
%                         upper-triangular, adjacency matrix specifying 
%                         edges surviving false discovery rate
%                         threshold
%       PVAL:             set to STATS.alpha
%   
%   Remarks:
%       If no edges survive the false discovery rate, CON_MAT and
%       PVAL are returned empty and N_CNT = 0.  
%
%       STATS.test_stat is empty if the number of permutations is too large
%       to precompute. See Limit parameter in NBSrun for details. In this 
%       situation, NBSglm is repeatedly called (J times) to compute test 
%       statistics for each permutation. 
%
%   azalesky@unimelb.edu.au

STATS=varargin{1}; 
if nargin==2
    H=varargin{2}; 
elseif nargin==3
    H=varargin{2};
    GLM=varargin{3}; 
end

%Number of nodes
N=STATS.N; 

%Number of edges
J=N*(N-1)/2;


if isempty(STATS.test_stat)
    %Compute test statistics on the fly
    test_stat=zeros(2,N*(N-1)/2);
    %Get desired number of permutations
    K=GLM.perms;
    %Set to 1, since NBSglm will be called separately for each permutation
    GLM.perms=1; 
else
    %Number of permutations
    K=size(STATS.test_stat,1)-1;
end

%Index of upper triangular elements of connectivity matrix
ind_upper=zeros(J,1); 
ind_upper=find(triu(ones(N,N),1)); 

pvals=zeros(1,J);
for i=1:K
    if isempty(STATS.test_stat)
        %Determine p-values for each edge by repeatedly calling NBSglm for 
        %each permutation.
        test_stat=ea_NBSglm(GLM);  
        %Added to v1.1.2 hanged < to <=
        pvals=pvals+(test_stat(1,:)<=test_stat(2,:));
    else
        %Determine p-values for each edge using precomputed test statistics
        %Added v1.1.2 changed < to <=
        pvals=pvals+(STATS.test_stat(1,:)<=STATS.test_stat(i+1,:));
    end
    try uiwaitbar(H,i/K); catch; end
end
pvals=pvals/K; 

%Sort p-values
ind_srt=zeros(1,J); 
[pvals,ind_srt]=sort(pvals);

%Simes procedure
pvals=(pvals<=(1:J)/J*STATS.alpha);
ind=find(pvals); 

if ~isempty(ind)
    %Maximum index
    ind=ind(length(ind));
    n_cnt=1;
    con_mat{1}=spalloc(N,N,ind);
    con_mat{1}(ind_upper(ind_srt(1:ind)))=1; 
    con_mat{n_cnt}=triu(con_mat{n_cnt},1);
    %Set pval to the alpha significance threshold 
    pval=STATS.alpha; 
else
    con_mat=[];
    n_cnt=0; 
    pval=[];
end



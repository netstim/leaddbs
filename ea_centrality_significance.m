function [ixes]=ea_centrality_significance(XYZV)




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
    
Peuc=squareform(ea_pdist(dat(:,1:3))); % distance matrix
Pval=squareform(ea_pdist(dat(:,4))); % distance matrix

% exp
zexp=2; % 1=zscore, 2=exp
switch zexp
    case 1
        Peuc(~logical(eye(size(Peuc,1))))=zscore(Peuc(~logical(eye(size(Peuc,1)))));
        Peuc(logical(eye(size(Peuc,1))))=1;
        Pval(~logical(eye(size(Pval,1))))=zscore(Pval(~logical(eye(size(Pval,1)))));
        Pval(logical(eye(size(Pval,1))))=1;
        [X,Y]=meshgrid(dat(:,4),dat(:,4));
        vals=mean(cat(3,X,Y),3);
        vals(:)=zscore(vals(:));
    case 2
        Peuc(~logical(eye(size(Peuc,1))))=1./exp(Peuc(~logical(eye(size(Peuc,1)))));
        Peuc(logical(eye(size(Peuc,1))))=1;
        Pval(~logical(eye(size(Pval,1))))=1./exp(Pval(~logical(eye(size(Pval,1)))));
        Pval(logical(eye(size(Pval,1))))=1;
        [X,Y]=meshgrid(dat(:,4),dat(:,4));
        vals=mean(cat(3,X,Y),3);
end
Pval=vals.*Pval;

P=Peuc.*Pval;

P(logical(eye(size(P,1))))=0;
P=sum(P); % degree centrality

    csize(iter,:)=P;%.*dat(:,4)';
    %m(iter)=max(P);
    
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
disp([num2str(sum(ixes)),' values identified above threshold at p>0.05 in conservative maximum testing.']);



%% pairwise testing:

nullmodel=csize(2:end,:);
nullmodel=sort(nullmodel,2,'descend'); % to be able to compare every 1st, 2nd, 3rd, etc. value with each other

nullmodel=sort(nullmodel,1,'descend'); % to be able to set up a threshold..

[realvals,ids]=sort(csize(1,:),'descend');
signodecnt=0;

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

ids=ea_fdr_bh(p,0.05); % calculate fdr after Benjamini & Hochberg

ixes=zeros(1,size(XYZV,1));
try
ixes(ids)=1;
end
ixes=logical(ixes);

disp([num2str(sum(ixes)),' values identified above threshold at p>0.05 in pairwise testing.']);
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





function [h crit_p adj_p]=ea_fdr_bh(pvals,q,method,report)
% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery 
%            rate (FDR) of a family of hypothesis tests. FDR is the expected
%            proportion of rejected hypotheses that are mistakenly rejected 
%            (i.e., the null hypothesis is actually true for those tests). 
%            FDR is a somewhat less conservative/more powerful method for 
%            correcting for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
% Usage:
%  >> [h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are 
%             significant (i.e., their null hypotheses are rejected).  If 
%             no p-values are significant, crit_p=0.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted 
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
% Example:
%   [dummy p_null]=ttest(randn(12,15)); %15 tests where the null hypothesis
%                                       %is true
%   [dummy p_effect]=ttest(randn(12,5)+1); %5 tests where the null
%                                          %hypothesis is false
%   [h crit_p adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%
%
% For a review on false discovery rate control and other contemporary
% techniques for correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis 
% of event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values
% 5/14/2013- D.H.J. Poot, Erasmus MC, improved run-time complexity

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        error('Some p-values are greater than 1.');
    end
end

if nargin<2,
    q=.05;
end

if nargin<3,
    method='pdep';
end

if nargin<4,
    report='no';
end

s=size(pvals);
if (length(s)>2) || s(1)>1,
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[dummy, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=(1:m)*q/m;
    wtd_p=m*p_sorted./(1:m);
    
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./(1:m));
    thresh=(1:m)*q/denom;
    wtd_p=denom*p_sorted./[1:m];
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

if nargout>2,
    %compute adjusted p-values
    adj_p=zeros(1,m)*NaN;
    [wtd_p_sorted, wtd_p_sindex] = sort( wtd_p );
    nextfill = 1;
    for k = 1 : m
        if wtd_p_sindex(k)>=nextfill
            adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
            nextfill = wtd_p_sindex(k)+1;
            if nextfill>m
                break;
            end;
        end;
    end;
    adj_p=reshape(adj_p(unsort_ids),s);
end

rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
end

if strcmpi(report,'yes'),
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep'),
        fprintf('FDR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end


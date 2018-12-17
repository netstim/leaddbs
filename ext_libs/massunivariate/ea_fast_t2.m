% fast_t2 - Computes a 2D matrix of independent sample t-tests relatively 
%           quickly.
%
% Usage:
% [p_values, t_scores, mn_dif, stder]=fast_t2(dataA,dataB,tail,verblevel);
%
% Inputs:
%  dataA   - 2D matrix of data for Group A
%  dataB   - 2D matrix of data for Group B
%
% Optional Inputs:
%  tail        - [1 | 0 | -1] If tail=1, the alternative hypothesis is that
%                the mean of Group A is greater than that of Group B.  If tail=0,
%                the alternative hypothesis is that the mean of the groups are different
%                than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%                is that the mean of Group A is less than that of Group B.
%                {default: 0}
%  verblevel   - An integer specifiying the amount of information you want
%                this function to provide about what it is doing during runtime.
%                 Options are:
%                    0 - quiet, only show errors, warnings, and EEGLAB reports
%                    1 - stuff anyone should probably know
%                    2 - stuff you should know the first time you start working
%                        with a data set {default value}
%                    3 - stuff that might help you debug (show all
%                        reports)
%
% Outputs:
%  p_values    - p-value of difference between groups at each time point 
%                and electrode (no correction for multiple comparisons)
%  t_scores    - t-score of difference between groups at each time point 
%                and electrode
%  mn_dif      - mean voltage difference between groups at each time point 
%                and electrode (Group A-Group B)
%  stder       - standard error of the mean voltage difference at each time
%                point and electrode
%
% Work based on:
% David Groppe
% May, 2010
% Kutaslab, San Diego

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%


function [p_values, t_scores, mn_dif, stder]=ea_fast_t2(dataA,dataB,tail,verblevel)

if nargin<2,
    error('You need to provide data.');
end

if nargin<3,
    tail=0; %default two-tailed test
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<4,
   verblevel=2; 
end

[n_chanA, n_subsA]=size(dataA);
[n_chanB, n_subsB]=size(dataB);
if n_chanA~=n_chanB,
    error('dataA and dataB have a different number of channels. They need to be the same.')
end


df=n_subsA+n_subsB-2;
if verblevel~=0,
    fprintf('fast_t2: Number of channels: %d\n',n_chanA);
    fprintf('fast_t2: Number of participants in Group A: %d\n',n_subsA);
    fprintf('fast_t2: Number of participants in Group B: %d\n',n_subsB);
    fprintf('t-score degrees of freedom: %d\n',df);
end

smA=sum(dataA,2);
mnA=smA/n_subsA;
ssA=sum(dataA.^2,2)-(smA.^2)/n_subsA;

smB=sum(dataB,2);
mnB=smB/n_subsB;
ssB=sum(dataB.^2,2)-(smB.^2)/n_subsB;

mult_fact=(n_subsA+n_subsB)/(n_subsA*n_subsB);
pooled_var=(ssA+ssB)/df;
stder=sqrt(pooled_var*mult_fact);

mn_dif=mnA-mnB;
t_scores=mn_dif./stder;

if tail<0,
    %lower tailed test
    p_values=tcdf(t_scores,df);
elseif tail>0,
    %upper tailed test
    p_values=1-tcdf(t_scores,df);
else
    %two tailed test
    p_values=tcdf(t_scores,df);
    ids=find(p_values>.5); %t-scores above zero
    p_values(ids)=1-p_values(ids); 
    p_values=p_values*2; %double for two tailed test
end
  


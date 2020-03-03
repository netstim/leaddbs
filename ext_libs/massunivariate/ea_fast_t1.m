% fast_t1 - Computes a 2D matrix of one sample/repeated measures t-tests 
%           relatively quickly.
%
% Usage:
% [p_values, t_scores, mn, stder]=fast_t1(data,tail,verblevel);
%
% Inputs:
%  data   - 2D matrix of data
%
% Optional Inputs:
%  tail        - [1 | 0 | -1] If tail=1, the alternative hypothesis is that the
%                mean of data is greater than 0 (upper tailed test).  If tail=0,
%                the alternative hypothesis is that the mean of data is different
%                than 0 (two tailed test).  If tail=-1, the alternative hypothesis
%                is that the mean of the data is less than 0 (lower tailed test).
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
%  p_values    - p-value at each time point and electrode (no correction
%                for multiple comparisons)
%  t_scores    - t-score at each time point and electrode
%  mn          - mean voltage at each time point and electrode
%  stder       - standard error of the mean voltage at each time point and 
%                electrode
%
% Work based on:
% David Groppe
% May, 2010
% Kutaslab, San Diego

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%


function [p_values, t_scores, mn, stder]=ea_fast_t1(data,tail,verblevel)

if nargin<1,
    error('You need to provide data.');
end

if nargin<2,
    tail=0; %default two-tailed test
elseif (tail~=0) && (tail~=1) && (tail~=-1),
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<3,
   verblevel=2; 
end

[n_chan, n_subs]=size(data);
df=n_subs-1;
if n_subs<2,
    error('You need data from at least two observations (e.g., participants) to perform a hypothesis test.')
end

if verblevel~=0,
    fprintf('fast_t1: Number of channels: %d\n',n_chan);
    fprintf('fast_t1: Number of participants: %d\n',n_subs);
    fprintf('t-score degrees of freedom: %d\n',df);
end

sm=sum(data,2);
mn=sm/n_subs;
sm_sqrs=sum(data.^2,2)-(sm.^2)/n_subs;
stder=sqrt(sm_sqrs/(n_subs*df));
t_scores=mn./stder;
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
  

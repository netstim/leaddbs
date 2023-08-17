%function [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY,n_perm,tail,alpha_level,stat,reports,seed_state)
%
% mult_comp_perm_corr-Permutation test based on Pearson's linear correlation
% coefficient (r) or Spearman's rank correlation coefficient (rho).  This
% function can perform the test on one or more pairs of variables.  When
% applying the test to multiple variables, the "max statistic" method is
% used for adjusting the p-values of each variable for multiple comparisons
% (Groppe, Urbach, & Kutas, 2011).  Like Bonferroni correction, this method
% adjusts p-values in a way that controls the family-wise error rate.
% However, the permutation method will be more powerful than Bonferroni
% correction when different variables in the test are correlated.
%
% Required Inputs:
%  dataX   - A vector or 2D matrix of data (Observation x Variable)
%  dataY   - A vector or 2D matrix of data (Observation x Variable).  Each
%            element of dataX should be paired with an element of dataY.
%
% Optional Inputs:
%  n_perm      - Number of permutations used to estimate the distribution of
%                the null hypothesis. If the number of observations is less
%                than or equal to 7, all possible permutations are used and
%                this optional input has no effect.  If the number of
%                observations is greater than 7, n_perm specifies the
%                number of random permutations computed. Manly (1997)
%                suggests using at least 1000 permutations for an alpha level
%                of 0.05 and at least 5000 permutations for an alpha level of
%                0.01. {default: 5000}
%  tail        - [1, 0, or -1] If tail=1, the alternative hypothesis is
%                that the correlation between dataX and dataY is positive
%                (upper tailed test). If tail=0, the alternative hypothesis
%                is that the correlation between dataX and dataY is
%                non-zero (two-tailed test). If tail=-1, the alternative
%                hypothesis is that the correlation between dataX and dataY
%                is negative (lower tailed test). {default: 0}
%  alpha_level - Desired family-wise alpha level. Note, because of the finite
%                number of possible permutations, the exact desired family-wise
%                alpha may not be possible. Thus, the closest approximation
%                is used and output as est_alpha. {default: .05}
%  stat        - 'linear', 'rank' or 'rankit'.  If 'linear', then Pearson's
%                r will be the correlation coefficient used.  If 'rank',
%                Spearman's rho will be used. If 'rankit', Bliss's rankit
%                correlation coefficient will be used. {default: 'linear'}
%  reports     - [0 or 1] If 0, function proceeds with no command line
%                reports. Otherwise, function reports what it is doing to
%                the command line. {default: 0}
%  seed_state  - The initial state of the random number generating stream
%                (see MATLAB documentation for "randstream"). If you pass
%                a value from a previous run of this function, it should
%                reproduce the exact same values. Note, this input only has
%                an effect if you have more than 12 observations.
%
% Outputs:
%  pval       - p-value (adjusted for multiple comparisons) of each
%               variable
%  corr_obs   - correlation coefficient of the observations for each variable
%  crit_corr  - Lower and upper critical correlation coefficients for given
%               alpha level. Correlations that exceed these values
%               significantly deviate from the null hypothesis.  For upper
%               tailed tests, the lower critical coefficient is NaN. The
%               opposite is true of lower tailed tests.
%  est_alpha  - The estimated family-wise alpha level of the test.  With
%               permutation tests, a finite number of p-values are possible.
%               This function tries to use an alpha level that is as close
%               as possible to the desired alpha level.  However, if the
%               sample size is small, a very limited number of p-values are
%               possible and the desired family-wise alpha level may be
%               impossible to approximately achieve.
%  seed_state - The initial state of the random number generating stream
%               (see MATLAB documentation for "randstream") used to
%               generate the permutations. You can use this to reproduce
%               the output of this function.  If the number of observations
%               is less than or equal to 7, seed_state='exact' since
%               all possible permutations are used in lieu of random
%               permutations.
%
% Note:
% -Unlike a parametric test (e.g., an ANOVA), a discrete set of p-values
% are possible (at most the number of possible permutations).  Since the
% number of possible permutations grows rapidly with the number of
% observations, this is only an issue for small sample sizes (e.g., 6
% independent observations).  When you have such a small sample size, the
% limited number of p-values may make the test overly conservative (e.g.,
% you might be forced to use an alpha level of .0286 since it is the biggest
% possible alpha level less than .05).
%
% -The null hypothesis of the permutation test is that every possible order
% of the observations in dataX and dataY are equally likely.  It is a
% non-parametric test because it makes no assumptions about the
% distributions of dataX and dataY.
%
%
% Example:
% >> dataX=randn(16,5); %5 variables, 16 observations
% >> dataY=randn(16,5);
% >> dataY(:,1)=.5*dataX(:,1)+dataY(:,1); %first variable of dataX and dataY are positively correlated
% >> dataY(:,3)=-dataX(:,3)+dataY(:,3); %third variable of dataX and dataY are negatively correlated
% >> [pval, corr_obs, crit_corr, est_alpha, seed_state]=mult_comp_perm_corr(dataX,dataY);
% >> disp(pval); %adjusted p-values
%
%
% References:
% Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis
% of event-related brain potentials/fields I: A critical tutorial review.
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x
%
% Manly, B.F.J. (1997) Randomization, Bootstrap, and Monte Carlo Methods in
% Biology. 2nd ed. Chapman and Hall, London.
%
%
% Author:
% David Groppe
% Jan, 2012
% Laboratory for Multimodal Human Brain Mapping
% Cushing Neuroscience Institute
% Manhasset, NY
%

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 3/17/2013-Randomization is now compatible with Matlab v13. Thanks to
% Aaron Newman for the fix.
%
% 8/14/2015-If setting seed state fails, keep going. This makes the
% function compatible with Matlab v14.

% Notes:
% Script test_corr_perm.m found this function to have accurate false postitive
% rates

function [corr_obs, pval, crit_corr, est_alpha, seed_state] = ea_permcorr(dataX,dataY,stat,n_perm,tail,alpha_level,rows,reports,seed_state)

if nargin<2
    error('You need to provide two sets of data.');
end

if nargin<3
    stat='rank';
end

if ~ismember(lower(stat), {'linear', 'pearson', 'rank', 'spearman', 'rankit'})
    error('Argument ''stat'' needs to be ''linear'', ''pearson'', ''rank'', ''spearman'' or ''rankit''');
end

if nargin<4
    n_perm=5000;
end

if nargin<5
    tail=0;
elseif (tail~=0) && (tail~=1) && (tail~=-1)
    error('Argument ''tail'' needs to be 0,1, or -1.');
end

if nargin<6
    alpha_level=0.05;
elseif (alpha_level>=1) || (alpha_level<=0)
   error('Argument ''alpha_level'' needs to be a number between 0 and 1.');
end

if nargin<7
    rows = 'complete';
end

if nargin<8
    reports = 0;
end

% Get random # generator state
defaultStream=RandStream.getGlobalStream;
if nargin<9 || isempty(seed_state)
    % Store state of random number generator
    seed_state=defaultStream.State;
else
    defaultStream.State=seed_state; % Reset random number generator to saved state
end

% Only use rows without NaN values
if strcmpi(rows, 'complete')
    ind = ~any(isnan(dataX),2) & ~any(isnan(dataY),2);
    if ~any(ind)
        ea_cprintf('CmdWinErrors', 'X or Y contains only NaNs!\n');
        corr_obs = nan;
        pval = nan;
        crit_corr = [nan nan];
        est_alpha = nan;
        return;
    elseif ~all(ind)
        ea_cprintf('CmdWinWarnings', 'X or Y is missing values, using only rows without NaNs.\n');
        dataX = dataX(ind);
        dataY = dataY(ind);
    end
end

[n_obsX, n_varX]=size(dataX);
if n_obsX<2
    error('You need data from at least two observations in data set X to perform a hypothesis test.')
end

[n_obsY, n_varY]=size(dataY);
if n_obsY<2
    error('You need data from at least two observations in data set Y to perform a hypothesis test.')
end

if n_obsY~=n_obsX
   error('You need to have the same number of observations in data set X and Y');
end

if n_varY~=n_varX
   error('You need to have the same number of variables in data set X and Y');
end

if n_obsX<=5
    n_psbl_prms=factorial(n_obsX);
    if reports
        watchit(sprintf(['Due to the very limited number of observations,' ...
            ' the total number of possible permutations is small.\nThus only a limited number of p-values (at most %d) are possible and the test might be overly conservative.'], ...
            n_psbl_prms));
    end
end

if reports
    fprintf('mult_comp_perm_corr: Number of variables: %d\n',n_varX);
    fprintf('mult_comp_perm_corr: Number of observations: %d\n',n_obsX);
end

if ismember(lower(stat), {'rank', 'spearman'})
    % Convert values to ranks
    dataX = ea_vec2rank(dataX);
    dataY = ea_vec2rank(dataY);
elseif strcmpi(stat, 'rankit') % https://doi.org/10.1177/0013164414557639
    % Convert values to ranks
    dataX = ea_vec2rank(dataX);
    dataY = ea_vec2rank(dataY);
    % Normalizing transformation
    dataX = norminv((dataX-0.5)/n_obsX);
    dataY = norminv((dataY-0.5)/n_obsX);
end

%% Set up permutation test
if n_obsX<=7
    n_perm=factorial(n_obsX); % Total number of possible permutations
    exact=1;
    seed_state='exact';
    if reports
        fprintf('Due to the limited number of observations, all possible permutations of the data will be computed instead of random permutations.\n');
    end
else
    exact=0;
end

if reports
    fprintf('Executing permutation test with %d permutations...\n',n_perm);
    fprintf('Permutations completed: ');
end

%% Compute permutations

% Constant factors for computing correlations
ssX=sum(dataX.^2)-(sum(dataX).^2)/n_obsX;
ssY=sum(dataY.^2)-(sum(dataY).^2)/n_obsX;
sum_dataX=sum(dataX);
sum_dataY=sum(dataY);
sqrt_ssYssX=sqrt(ssY.*ssX);
constant=sum_dataX.*sum_dataY/n_obsX;
% ssXY=sum(dataX.*dataY)-sum(dataX).*sum(dataY)/n_obsX;
ssXY=sum(dataX.*dataY)-constant;

%% Computes t-scores of observations at all variables and time points
corr_obs=ssXY./sqrt_ssYssX;

%%
if exact
    % Use all possible permutations
    all_perms=perms(1:n_obsX);
    mx_corr=zeros(1,n_perm);
    for perm=1:n_perm
        if ~rem(perm,100)
            if reports
                if ~rem(perm-100,1000)
                    fprintf('%d',perm);
                else
                    fprintf(', %d',perm);
                end
                if ~rem(perm,1000)
                    fprintf('\n');
                end
            end
        end
        % Randomly set order of dataX
        ssXY=sum(dataX(all_perms(perm,:),:).*dataY)-constant;
        r_perm=ssXY./sqrt_ssYssX;

        % Get most extreme correlation
        if tail==1
            mx_corr(perm)=max(r_perm);
        elseif tail==-1
            mx_corr(perm)=min(r_perm);
        else
            [mx_corr(perm) id]=max(abs(r_perm));
            mx_corr(perm)=mx_corr(perm)*sign(r_perm(id)); % Add sign back
        end
    end
else
    % Use random permutations
    mx_corr=zeros(1,n_perm);
    for perm=1:n_perm
        if ~rem(perm,100)
            if reports
                if ~rem(perm-100,1000)
                    fprintf('%d',perm);
                else
                    fprintf(', %d',perm);
                end
                if ~rem(perm,1000)
                    fprintf('\n');
                end
            end
        end
        % Randomly set order of dataX
        orderX=randperm(n_obsX);
        ssXY=sum(dataX(orderX,:).*dataY)-constant;
        r_perm=ssXY./sqrt_ssYssX;

        % Get most extreme correlation
        if tail==1
            mx_corr(perm)=max(r_perm);
        elseif tail==-1
            mx_corr(perm)=min(r_perm);
        else
            [mx_corr(perm), id]=max(abs(r_perm));
            mx_corr(perm)=mx_corr(perm)*sign(r_perm(id)); % Add sign back
        end
    end
end

% End permutations completed line
if reports && rem(perm,1000)
    fprintf('\n');
end

%% Compute p-values
pval=zeros(1,n_varX);
for t=1:n_varX
    if tail==0
        pval(t)= ( sum( abs(mx_corr) > abs(corr_obs(t)) ) + 1 ) / (n_perm + 1);
    elseif tail==1
        pval(t)= ( sum( mx_corr > corr_obs(t) ) + 1 ) / (n_perm + 1);
    elseif tail==-1
        pval(t)= ( sum( mx_corr < corr_obs(t) ) + 1 ) / (n_perm + 1);
    end
end

%% Compute critical t-scores for specified alpha level
if tail==0
    %two-tailed
    crit_corr(1)=prctile(mx_corr,100*alpha_level/2);
    crit_corr(2)=prctile(mx_corr,100-100*alpha_level/2);
    est_alpha=mean(mx_corr>=crit_corr(2))+mean(mx_corr<=crit_corr(1));
elseif tail==1
    %upper tailed
    crit_corr(1)=NaN;
    crit_corr(2)=prctile(mx_corr,100-100*alpha_level);
    est_alpha=mean(mx_corr>=crit_corr(2));
else
    %tail=-1, lower tailed
    crit_corr(1)=prctile(mx_corr,alpha_level*100);
    est_alpha=mean(mx_corr<=crit_corr(1));
    crit_corr(2)=NaN;
end
if reports
    fprintf('Desired family-wise alpha level: %f\n',alpha_level);
    fprintf('Estimated actual family-wise alpha level for returned values of crit_corr: %f\n',est_alpha);
end


function watchit(msg)
%function watchit(msg)
%
% Displays a warning message on the Matlab command line.  Used by
% several Mass Univariate ERP Toolbox functions.
%

disp(' ');
disp('****************** Warning ******************');
disp(msg);
disp(' ');

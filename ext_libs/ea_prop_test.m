function [h,p, chi2stat,df] = ea_prop_test(X , N, correct)

% [h,p, chi2stat,df] = prop_test(X , N, correct)
% 
% A simple Chi-square test to compare two proportions
% It is a 2 sided test with alpha=0.05
%
% Input:
% X = vector with number of success for each sample (e.g. [20 22])
% N = vector of total counts for each sample (e.g. [48 29])
% correct = true/false : Yates continuity correction for small samples?
%
% Output:
% h = hypothesis (H1/H0)
% p = p value
% chi2stat= Chi-square value
% df = degrees of freedom (always equal to 1: 2 samples)
%
% Needs chi2cdf from the Statistics toolbox
% Inspired by prop.test() in "R" but much more basic
%
% Example: [h,p,chi]=prop_test([20 22],[48 29], true)
% The above example tests if 20/48 differs from 22/29 using Yate's correction

if (length(X)~= 2)||(length(X)~=length(N))
    disp('Error: bad vector length')
elseif (X(1)>N(1))|| (X(2)>N(2))
    disp('Error: bad counts (X>N)')
else  
    df=1; % 2 samples

    % Observed data
    n1 = X(1);
    n2 = X(2);
    N1 = N(1);
    N2 = N(2);

    if ((n1/N1)-(n2/N2))>0
        multiplyby=1;
    else
        multiplyby=-1;
    end

    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    
    if correct == false
        % Standard Chi-square test
        chi2stat = sum((observed-expected).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    else
        % Yates continuity correction        
        chi2stat = sum((abs(observed - expected) - 0.5).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    end
    chi2stat=chi2stat*multiplyby;
    h=0;
    if p<0.05
        h=1;
    end
end
end
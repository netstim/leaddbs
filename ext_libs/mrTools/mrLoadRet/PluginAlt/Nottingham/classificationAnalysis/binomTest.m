function pout=binomTest(s,n,p,Sided)
%function pout=binomTest(s,n,p,Sided)
%
% Performs a binomial test of the number of successes given a total number 
% of outcomes and a probability of success. Can be one or two-sided.
%
% Inputs:
%       s-      (Scalar or Array) The observed numebr of successful outcomes
%       n-      (Scalar or Array) The total number of outcomes (successful or not)
%       p-      (Scalar or Array) The proposed probability of a successful outcome
%       Sided-  (String) can be 'Two','Greater' or 'Lesser'. A value of
%               'Two' will perform a two-sided test, that the actual number
%               of success is different from the expected number of
%               successes in any direction. 'Greater' or 'Lesser' will
%               perform a 1-sided test to examine if the observed number of 
%               successes are either significantly greter than or less than
%               (respectively) the expected number of successes.
%
% Outputs:
%       pout-   The probability of observing the resulting value of s or a
%               value more extreme (the precise meaning of which depends on
%               the value of Sided) given n total outcomes with a 
%               probability of success of p.               
%
%       s, n and p can be scalars or arrays of the same size. The
%       dimensions and size of pout will match that of these inputs.
%
% For example, the signtest is a special case of this where the value of p
% is equal to 0.5 (and a 'success' is dfeined by whether or not a given
% sample is of a particular sign.), but the binomial test and this code is 
% more general allowing the value of p to be any value between 0 and 1.
%
% References:
% http://en.wikipedia.org/wiki/Binomial_test
%
% Adapted from myBinomTest.m
% by Matthew Nelson July 21st, 2009
% matthew.j.nelson.vumail@gmail.com

if nargin<4 || isempty(Sided);    Sided='Two';      end
if nargin<3 || isempty(p);      p=0.5;      end
    
switch lower(Sided)
    case 'two'
        [s,n,p]= EqArrayAndScalars(s,n,p);

        E=p.*n;

        GreaterInds=s>=E;
        pout=repmat(0,size(GreaterInds));
        dE=pout;

        %note that matlab's binocdf(s,n,p) gives the prob. of getting up to AND INCLUDING s # of successes...
        %Calc pout for GreaterInds first
        pout(GreaterInds)=1-binocdf( s(GreaterInds)-1,n(GreaterInds),p(GreaterInds));  %start with the prob of getting >= s # of successes

        %now figure the difference from the expected value, and figure the prob of getting lower than that difference from teh expected value # of successes
        dE(GreaterInds)=s(GreaterInds)-E(GreaterInds);
        pout(GreaterInds)=pout(GreaterInds)+ binocdf(floor(E(GreaterInds)-dE(GreaterInds)),n(GreaterInds),p(GreaterInds));    %the binonmial is a discrete dist. ... so it's value over non-integer args has no menaing... this flooring of E-dE actually doesn't affect the outcome (the result is teh same if the floor was removed) but it's included here as a reminder of the discrete nature of the binomial

        %Calc pout for LesserInds second
        pout(~GreaterInds)=binocdf(s(~GreaterInds),n(~GreaterInds),p(~GreaterInds));  %start with the prob of getting <= s # of successes

        %now figure the difference from the expected value, and figure the prob of getting greater than that difference from teh expected value # of successes
        dE(~GreaterInds)=E(~GreaterInds)-s(~GreaterInds);
        pout(~GreaterInds)=pout(~GreaterInds) + 1-binocdf(ceil(E(~GreaterInds)+dE(~GreaterInds))-1,n(~GreaterInds),p(~GreaterInds));    %Here teh ceiling is needed b/c of teh -1 following it, so that integer and non-integer vals of E+dE will bothe give teh correct value with teh same line of code
    case 'greater'  %one-sided
        pout=1-binocdf(s-1,n,p);  %just report the prob of getting >= s # of successes
    case 'lesser'   %one-sided
        pout=binocdf(s,n,p);  %just report the prob of getting <= s # of successes
    otherwise
        error(['In myBinomTest, Sided variable is: ' Sided '. Unkown sided value.'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout=EqArrayAndScalars(varargin)
%function varargout=EqArrayAndScalars(varargin)
%
% This will compare a collection of inputs that must be either scalars or 
% arrays of the same size. If there is at least one array input, all scalar
% inputs will be replicated to be the array of that same size. If there are
% two or more array inputs that have different sizes, this will return an
% error.
%
% created by Matthew Nelson on April 13th, 2010
% matthew.j.nelson.vumail@gmail.com                        

d=repmat(0,nargin,1);
for ia=1:nargin    
    d(ia)=ndims(varargin{ia});
end
maxnd=max(d);

s=repmat(1,nargin,maxnd);
    
for ia=1:nargin
    s(ia,1:d(ia))=size(varargin{ia});
end
maxs=max(s);

for ia=1:nargin
    if ~all(s(ia,:)==maxs)
        if ~all(s(ia,:)==1)
            error(['Varargin{' num2str(ia) '} needs to be a scalar or equal to the array size of other array inputs.'])
        else
            varargout{ia}=repmat(varargin{ia},maxs);
        end
    else
        varargout{ia}=varargin{ia};
    end
end

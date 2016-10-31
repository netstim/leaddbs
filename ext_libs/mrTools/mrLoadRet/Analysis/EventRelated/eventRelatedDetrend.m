% mydetrend.m
%
%      usage: mydetrend(data,<demean>,<dispfig>)
%         by: justin gardner
%       date: 07/05/05
%    purpose: detrend a vector. if data is
%             a matrix then the detrending goes
%             along the columns.
%       e.g.: mydetrend(rand(1,51)+(0:1:50)*0.1+3,1)
%
function retval = eventRelatedDetrend(data,demean,dispfig)

% check command line arguments
if (nargin == 1)
  demean = 0;
  dispfig = 0;
elseif (nargin == 2)
  dispfig = 0;
elseif (nargin ~= 3)
  help mydetrend;
  return
end

% make into column vector
if (size(data,1)) == 1
  data = data';
end

% get length
n = size(data,1);

% make regression matrix
A = [(0:1/(n-1):1)' ones(n,1)];

% get the slope and offsets
regcoef = ((A'*A)^-1)*A'*data;

% take out the slope
if demean
  retval = data-A*regcoef;
else
  retval = data-A(:,1)*regcoef(1,:);

end

% plot it
if dispfig
  mlrSmartfig('mydetrend');
  for i = 1:size(data,2)
    plot(data,'g');
    hold on
    plot(A*regcoef,'k-');
    plot(retval,'r');
    legend('raw data','regression line','detrended data');
  end
end

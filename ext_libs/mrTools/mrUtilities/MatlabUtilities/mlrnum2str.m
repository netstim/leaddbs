% mynum2str.m
%
%        $Id:$ 
%      usage: mynum2str(num,<sigfigs=2>,<doFixBadChars=false>,<tabs=false>,'compact=true')
%         by: justin gardner
%       date: 09/07/09
%    purpose: num2str that allows setting # of significant figures and also doesn't make so many spaces
%             but still will align numbers from line to line.
%
%             For example, if you set sigfigs=-1 then the program will figure out how many sigfigs are needed
%             to display all numbers and make them align across lines. For example, the following
%             line will produce:
%       e.g.: mynum2str([12.1 -0.001 10.1;-1.3 -12.4 30.01],'sigfigs=-1','compact',false)
%              12.100  -0.001  10.100
%              -1.300 -12.400  30.010
%
%             set sigfigs to the number of significant digits you want to show (numbers are rounded to
%               show the appropriate number of sigfigs).
%             set tabs to true if you want to have each number followed by a tab
%             set compact to false if you want numbers to align from row to row (like in the above example)
%                 the default gives the most compact string possible
%
function s = mynum2str(num,varargin)

% check arguments
if nargin == 0
  help mynum2str
  return
end

% check for empty num, in which case return empty
s = '';
if isempty(num),return,end

% evaluate arguments
sigfigs = [];
doFixBadChars = [];
tabs = [];
compact = [];
getArgs(varargin,{'sigfigs=2','doFixBadChars',false,'tabs',false,'compact',true});

% automatic sigfig
sigFigsEachNum = [];

if sigfigs == -1
  % get how many sigfigs each number needs
  sigFigsEachNum = getSigFigs(num);
  % and get the maximum needed sigfigs
  sigfigs = max(sigFigsEachNum(:));
end

% need sigFigsEachNum if we are doing a compact display
if compact && isempty(sigFigsEachNum)
  sigFigsEachNum = getSigFigs(num,sigfigs);
end

% intialize return
s = '';

% check for non 1d or 2d array
if length(size(num))>2
  disp(sprintf('(mynum2str) Can not display %i dimensional matrix',length(size(s))));
  return
elseif length(size(num)) == 2
  num = num';
  % compute maxnumdigits on the left side of decimal excluding inf and nan
  maxnumdigits = length(sprintf('%i',floor(nanmax(abs(num(~isinf(num(:))))))));
  % if there is an inf or nan, account for that
  if any(isnan(num(:))) | any(isinf(num(:)))
    % if sigfigs is zero then Inf and Nan cannot line up with decimal, so shift over by a space
    if sigfigs == 0
      maxnumdigits = max(maxnumdigits,3);
      if any(num(isinf(num(:)))<0)
	maxnumdigits = max(maxnumdigits,4);
      end
    else
      maxnumdigits = max(maxnumdigits,2);
      if any(num(isinf(num(:)))<0)
	maxnumdigits = max(maxnumdigits,4);
      end
    end
  end
end

% make the string
for j = 1:size(num,2)
  for i = 1:size(num,1)
    % create the formatting string
    if compact
      % create string
      formatString = sprintf('%%s%%0.0%if ',sigFigsEachNum(j,i));
    else
      % figure out how many spaces to add
      if num(i,j) < 0,addspace = '';else addspace = ' ';end
      % get the number of digits to the left of decimal point
      if isnan(num(i,j)) | isinf(num(i,j))
	% for nan or inf, then set spaces so that the Inf or NaN string lines up with decimal
	if sigfigs > 0
	  numDigitsToLeftOfDecimal = 2;
	else
	  numDigitsToLeftOfDecimal = 3;
	end
      else
	% for all other numbers count how many places it takes to represent 
	numDigitsToLeftOfDecimal = length(sprintf('%i',floor(abs(num(i,j)))));
      end
      addspace = [addspace repmat(' ',1,maxnumdigits-numDigitsToLeftOfDecimal)];
      % create string
      formatString = sprintf('%%s%s%%0.0%if ',addspace,sigfigs);
      % add spaces after a nan or inf
      if isnan(num(i,j)) | isinf(num(i,j))
	formatString = sprintf('%s%s',formatString,repmat(' ',1,sigfigs));
      end
    end
    % add tabs if called for
    if tabs
      formatString = sprintf('%s\t',formatString(1:end-1));
    end
    % and update the full string
     s = sprintf(formatString,s,round(num(i,j)*(10^sigfigs))/(10^sigfigs));
  end
  % strip off last space
  s = s(1:end-1);
  % add new line character
  if j ~= size(num,2),s = sprintf('%s\n',s);end
end

% fix bad chars
if doFixBadChars
  s = fixBadChars(s,[],{'.','p'});
  s = s(2:end);
end

%%%%%%%%%%%%%%%%%%%%
%    getSigFigs    %
%%%%%%%%%%%%%%%%%%%%
function retval = getSigFigs(num,maxSigFigs)

% if we want at most 0 sigfigs then just return zero for all elements
if (nargin>=2) && (maxSigFigs == 0)
  retval = zeros(size(num));
  return
end

if nargin == 1
  % maximum number of sigfigs
  maxSigFigs = 6;
  minDiff = 1e-10;
else
  % get the minimum difference that is still considered the same number
  % this is a bit of a hack since at some point due to numerical round
  % off you might be smaller than this limit, but still be different numbers
  % but, this should only mean that you will get less digits than asked
  % for only in the compact case
  minDiff = 10^(-2*maxSigFigs);
end

for i = 1:size(num,1)
  for j = 1:size(num,2)
    % for nan and inf the sigfigs needed are 0
    if isnan(num(i,j)) | isinf(num(i,j))
      retval(i,j) = 0;
    else
      % otherwise find out how many sigfigs are needed
      for iSigfig =  1:maxSigFigs
	% check to see if this number is evenly roundable by this
	% many sig digits
	if (abs(round(num(i,j)*(10^iSigfig))/(10^iSigfig) - num(i,j))) < minDiff
	  break
	end
      end
      retval(i,j) = iSigfig;
    end
  end
end


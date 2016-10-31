% readpetable.m
%
%      usage: readpetable(filename)
%         by: justin gardner
%       date: 08/11/11
%    purpose: 
%
function petable = readpetable(filename,varargin)

petable = [];
% check arguments
if nargin < 1
  help readpetable
  return
end

% petable
petablePath = [];verbose = [];
getArgs(varargin,{'petablePath=~/vnmrsys/tablib','verbose=0'});

% look in the following paths for the petabel
possiblePaths = {'.',petablePath,'/usr4/justin/vnmrsys/tablib'};
foundpetable = false;
for iPath  = 1:length(possiblePaths)
  % see if petable exists
  if isfile(fullfile(possiblePaths{iPath},filename))
    foundfile = true;
    petablePath = possiblePaths{iPath};
    continue;
  end
end
if ~foundfile
  disp(sprintf('(readpetabe) Could not find petable %s',filename));
  return
end

% open petable
fPetable = fopen(fullfile(petablePath,filename),'r');

l = fgetl(fPetable);
while ischar(l)
  % see if it is a comment or empty
  if isempty(findstr('/*',l)) && ~isempty(l)
    % see if we have a variable name
    if ~isempty(findstr('=',l))
      varname = textscan(l,'%s =',1);
      varname = varname{1}{1};
      petable.(varname) = [];
    else
      % must be a line of numbers, get them
      varval = textscan(l,'%d');
      petable.(varname)(end+1,:) =  double(varval{1});
    end
  end
  l = fgetl(fPetable);
end

% close petable
fclose(fPetable);

if verbose
  disp(sprintf('(readpetable) Read petable: %s',fullfile(petablePath,filename)));
end

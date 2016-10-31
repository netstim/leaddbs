% readfdf.m
%
%        $Id:$ 
%      usage: [h d] = readfdf(FDFdirname,<verbose>,<headerOnly>);
%         by: justin gardner
%       date: 08/08/11
%    purpose: Reads Varian FDF (Flexible Data Format) Files, based on VNMRJ manual
%             User Programming, Parameters and Data, Page 365
%
function [d h] = readfdf(FDFdirname,verbose,headerOnly)

h = [];
d = [];

% check arguments
if ~any(nargin == [1 2 3])
  help readfdf
  return
end

% set verbose default
if nargin < 2,verbose = 0;headerOnly = false;end
if nargin < 3,headerOnly = false;end

% check FDF dir
if ~isdir(FDFdirname)
  FDFdirname = setext(FDFdirname,'img');
  if ~isdir(FDFdirname)
    disp(sprintf('(readfdf) Could not find FDF dir %s',FDFdirname));
    return
  end
end

% get file listing
FDFdir = dir(fullfile(FDFdirname,'*image*.fdf'));

% check that we got FDF files
if isempty(FDFdir)
  % try to see if the user meant an img directory
  FDFdirname = setext(FDFdirname,'img');
  FDFdir = dir(fullfile(FDFdirname,'*image*.fdf'));
  if isempty(FDFdir)
    disp(sprintf('(readfdf) Could not find any FDF files in %s',FDFdirname));
    return
  end
end

% read the image
[d h] = readfdffile(fullfile(FDFdirname,FDFdir(1).name),verbose,headerOnly);
if headerOnly,return,end

% if data is stored as slices than continue to read each slice in turn
if length(FDFdir)>1
  %  convert h to a cell array
  htemp = h;clear h;h(1) = htemp;
  for i = 2:length(FDFdir)
    [d(:,:,i) h(i)] = readfdffile(fullfile(FDFdirname,FDFdir(i).name),verbose,headerOnly);
  end
end

%%%%%%%%%%%%%%%%%%%%%
%    readfdffile    %
%%%%%%%%%%%%%%%%%%%%%
function [d h] = readfdffile(filename,verbose,headerOnly);

d = [];

% open file
f = fopen(filename,'r');

% check for error opening
if f==-1
  disp(sprintf('(readfdf:readfdffile) Error opening file: %s',filename));
  return
end

% read magic bytes for check
magicbytes = fread(f,2,'char');
if ((char(magicbytes(1))~='#') || (char(magicbytes(2))~='!'))
  disp(sprintf('(readfdf:readfdffile) ERROR: Magic bytes are not #! for file %s',filename));
  fclose(f);
  return
end

% read the idString - this should always be the string /usr/local/fdf/startup
idString = fgetl(f);
if ~strcmp(idString,'/usr/local/fdf/startup')
  disp(sprintf('(readfdf:readfdffile) !!!! id string is %s and not /usr/local/fdf/startup as expected',idString));
  % just warn
end

% some verbose info
if verbose,disp(sprintf('(readfdf:readfdffile) Reading file %s',filename));end

% now read line by line until we hit the end of header character (ASCII NUL)
thisline = fgetl(f);
while ~isempty(thisline) && (length(thisline) ~= 1)
  % get the varname and varvalue
  [varname varvalue] = strread(thisline,'%s%s','delimiter','=;');
  if length(varname)
    varname = varname{1};
    [vartype varname] = strread(varname,'%s%s','delimiter',' ');
    if length(varname), varname = varname{1}; else varvalue = [];end
  end

  % make sure we got a valid value, then decide whether to save or not
  if length(varvalue)
    % strip * and [] 
    if strcmp(varname(1),'*'),varname = varname(2:end);end
    if strcmp(varname(end-1:end),'[]') varname = varname(1:end-2);end
    % get varvalue
    varvalue = varvalue{1};
    % convert " to '
    quoteloc = find(varvalue=='"');
    varvalue(quoteloc) = '''';
    % now evaluate variable
    h.(varname) = eval(varvalue);
    % convert cell arrays if they don't contain strings
    if iscell(h.(varname)) 
      if ~isstr(h.(varname){1})
	h.(varname) = cell2mat(h.(varname));
      end
    end
  end
  thisline = fgetl(f);
end

% check for compressed data
if isfield(h,'compression') && ~isempty(h.compression)
  disp(sprintf('(readfdf) Compressed data reading is not implemented yet!!!'));
  fclose(f);
  return
end

% now make sure we got certain fields
fieldcheck = {'storage','type','bits'};
for i = 1:length(fieldcheck)
  if ~isfield(h,fieldcheck{i})
    disp(sprintf('(readfdf) Could not find field %s in header',fieldcheck{i}));
    fclose(f);
    return
  end
end

% now get the data type
storageType = find(strcmp(h.storage,{'integer','float'}));
if isempty(storageType)
  disp(sprintf('(readfdf) Unrecogonized data type storage=%s',h.storage));
  fclose(f);
  return
end

% set the storageType to the appropriate matlab name for the correct number of bits and int/float
if storageType == 1
  switch h.bits
   case {8,16,32,64}
    storageType = sprintf('int%i',h.bits);
   otherwise
    disp(sprintf('(readfdf) Unrecogonized number of bits %i',h.bits));
    fclose(f);
    return
  end
elseif storageType == 2
  switch h.bits
   case {8,16}
    disp(sprintf('(readfdf) Number of bits for float: %i not supported',h.bits));
    fclose(f);
    return
   case {32}
    storageType = 'float';
   case {64}
    storageType = 'double';
   otherwise
    disp(sprintf('(readfdf) Unrecogonized number of bits %i',h.bits));
    fclose(f);
    return
  end
end

% compute data size
dataSize = prod(h.matrix)*h.bits/8;

% see how much data is left in file
currentPos = ftell(f);
fseek(f,0,1);
bytesInFile = ftell(f)-currentPos;

% read header only
if headerOnly
  fclose(f);
  return
end

if bytesInFile < dataSize
  disp(sprintf('(readfdf) Found %i bytes in file, but need %i\n',bytesInFile,dataSize));
  fclose(f);
  return
end

% seek back from end the number of bytes needed
fseek(f,-dataSize,1);

% now read the data
d = fread(f,prod(h.matrix),sprintf('%s=>double',storageType));

% reshape
d = reshape(d,h.matrix);

% close file
fclose(f);

% readprocpar.m
%
%      usage: readprocpar.m()
%         by: justin gardner
%       date: 07/03/03
%    purpose: reads procpar into matlab variable
%      usage: procpar = readprocpar;
%
function procpar = readprocpar(procdir,verbose)


% argument check
if (nargin == 0)
  procdir = './';
  verbose = 0;
elseif (nargin == 1) 
  if ~strcmp(getext(procdir),'par')
    procdir = setext(procdir,'fid',0);
  end
  verbose = 0;
elseif (nargin ~= 2)
  help readprocpar;
  return
end
if (nargin >= 1)
  if (procdir(length(procdir)) ~= '/')
    procdir(length(procdir)+1) = '/';
  end
end

if (verbose)
  disp('(readprocpar) Reading procpar...');
end
procpar = [];

% open the procpar
fprocpar = fopen([procdir 'procpar'],'r');
if (fprocpar == -1)
  disp(sprintf('(readprocpar) ERROR: Could not open %s',[procdir 'procpar']));
  return
end

% scan through procpar looking for parameters needed for parsing
line = fgets(fprocpar);
while (line ~= -1) 
  % get the first token of the line
  token = strtok(line);
  % if it is a character and not doublequote then we have an array name
  if (~((token(1) >= '0') & (token(1) <= '9')) & (token(1) ~= '"'))
    fieldname = token;
    % get next line
    line = fgets(fprocpar);
    % first number is just length of array
    [len line] = strtok(line);
    len = str2num(len);
    % just an array of numbers
    if ((length(line) > 2) & (line(2) ~= '"'))
      eval(sprintf('procpar.%s = str2num(line);',fieldname));
      % confirm length
      arraylen = eval(sprintf('length(procpar.%s)',fieldname));
      if (arraylen ~= len)
	sprintf('WARNING: %s should be %i elements, but found only %i',fieldname,len,arraylen);
      end
    % character arrays
    else
      for j = 1:len
	% strip leading white space
	while (length(line) & (line(1) == ' '))
	  line = line(2:length(line));
	end
	token = line;
	% make sure we have a valid string token
	if (isempty(token) | (token(1) ~= '"'))
	  sprintf('WARNING: %s missing string element %i/%i',fieldname,j,len);
	  eval(sprintf('procpar.%s{%i} = '''';',fieldname,j));
	  j = len;
	else
	  str = strtok(line,sprintf('"\n'));
	  % convert single to double quotes
	  str(findstr(str,'''')) = '"';
	  eval(sprintf('procpar.%s{%i} = ''%s'';',fieldname,j,str));
	end
	line = fgets(fprocpar);
      end
    end
  end
  line = fgets(fprocpar);
end

fclose(fprocpar);

% these fields will get set from the petable name if available
if ~isfield(procpar,'accfactor'),procpar.accfactor = 1;end
if ~isfield(procpar,'tSesnseAccfactor'),procpar.tSenseAccfactor = 1;end
if ~isfield(procpar,'numshots'),procpar.numshots = 1;end

% for epi images, there are navigator echos, which
% should be subtracted from the number of lines.
% this can be known from the name of the petable
if (isfield(procpar,'petable'))
  token = procpar.petable{1};
  % the petable name should be something like
  % "epi132alt8k". We want the second number
  if (strncmp(token,'epi',3)) 
    j = 1;
    % go past any intial characters
    while((j < length(token)) && (isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % then skip the numbers
    while((j < length(token)) && (~isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % and skip the next characters
    while((j < length(token)) && (isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % find end of numbers
    k = j;
    while((k < length(token)) && (~isempty(strfind('0123456789',token(k))))),k=k+1;,end
    procpar.numshots = str2num(token(j:k-1));
    % then skip some more numbers
    while((j < length(token)) && (~isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % go past any other charcters (should be k)
    while((j < length(token)) && (isempty(strfind('0123456789',token(j))))),j=j+1;,end
    % and if we have anything left then get the number, that is the acceleration factor
    k = j;
    while((k < length(token)) && (~isempty(strfind('0123456789',token(k))))),k=k+1;,end
    if j < length(token)
      % r means that it is a sense petable
      if token(j+1:k) == 'r'
	procpar.accfactor = str2num(token(j:k-1));
      % t means t-sense acc factor	
      elseif token(j+1:k) == 't'
	procpar.tSenseAccfactor = str2num(token(j:k-1));
      end
    end
    % before we always had 1 navecho for epi, but now there is a field procpar.navecho that should
    % be set by Ken's prep program that reads the petable and gets the t3 field.
    if ~isfield(procpar,'navecho'),procpar.navecho = 1;end
  end
end

% the number of shots should now be read from the field nseg rather than from the petable
if isfield(procpar,'nseg') procpar.numshots = procpar.nseg;end

% if navecho is not set by now, then there were no navecho's
if ~isfield(procpar,'navecho'),procpar.navecho = 0;end

% compute the number of echoes
procpar.navechoes = procpar.navecho*procpar.numshots/procpar.accfactor;

if (verbose),disp('(readprocpar) Finished.');,end
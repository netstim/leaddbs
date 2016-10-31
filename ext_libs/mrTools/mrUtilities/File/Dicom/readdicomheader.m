% readdicomheader.m
%
%      usage: readdicomheader(filename)
%         by: justin gardner
%       date: 02/08/05
%
function retval = readdicomheader(filename)

if (nargin ~= 1)
  help readdicomheader;
  return
end

% check if file exists
if (~isfile(filename))
  disp(sprintf('UHOH: Could not find file %s',filename));
  return
end

% open file
fdicom = fopen(filename,'r');

% get each line and process
line = fgets(fdicom);
linenum = 1;
while (line ~= -1)
  % ignore everything that does not start with a number
  if ((length(line) > 4) && ~isempty(str2num(line(1:4))))
    % get group and element number
    [group line] = strtok(line);group = str2num(group);
    [element line] = strtok(line);element = str2num(element);
    [len line] = strtok(line);len = str2num(len);
    % make sure that these are valid
    if (~isempty(group) && ~isempty(element))
      % find file character string
      idloc = strfind(line,'//');
      if (length(idloc) >= 2)
	idstr = line(idloc(1)+2:idloc(2)-1);
	valstr = line(idloc(2)+2:length(line)-1);
	% check to see if the value length matches the string
	if (length(valstr) ~= len)
%	  disp(sprintf('UHOH: line %i value length %i does not match %i',linenum,len,length(valstr)));
	end
	% convert the idstr into a field name by
	% removing white space
	idname = '';
	[groupname idstr] = strtok(idstr);
	while(~isempty(idstr))
	  [idtok idstr] = strtok(idstr);
	  idname = sprintf('%s_%s',idname,idtok);
	end
	if (~isempty(groupname) && (length(idname) > 1))
	  % remove initial underscore
	  idname = idname(2:length(idname));
	  % remove punctuation
	  nonalpha = '!@#$%^&*()+~`''/';
	  for k = 1:length(nonalpha)
	    idname(strfind(idname,nonalpha(k))) = '_';
	  end
	  % if this is a string value, save as string
	  if (length(strfind(valstr,'\')) ~= 0) || (length(str2num(valstr))~=1)
	    eval(sprintf('retval.%s.%s=''%s'';',groupname,idname,valstr));
	  % otherwise save as number
	  else
	    valnum = str2num(valstr);
	    if (floor(valnum*10)/10 == valnum)
	      eval(sprintf('retval.%s.%s=%i;',groupname,idname,valnum));
	    else
	      eval(sprintf('retval.%s.%s=%f;',groupname,idname,valnum));
	    end
	  end
	end
      end
    end
  end
  linenum = linenum+1;
  line = fgets(fdicom);
end
fclose(fdicom);
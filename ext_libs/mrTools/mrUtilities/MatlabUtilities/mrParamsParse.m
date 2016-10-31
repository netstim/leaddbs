% mrParamsParse.m
%
%      usage: mrParamsParse(vars)
%         by: justin gardner
%       date: 03/13/07
%    purpose: called by mrDefaultParamsGUI and mrDefaultParamsReconcile
%             parses the variable information, see wiki for details on
%             format of vairable string
%
%	$Id$

function [vars varinfo] = mrParamsParse(vars)

% check arguments
if ~any(nargin == [1])
  help mrParamsParse
  return
end

%-------------------------------- first parse the argument
varinfo = cell(size(vars));
for i = 1:length(vars)
  % if the variable is just a string, then
  % it got passed in without a default argument
  % so make it into a cell array with the second
  % element set to empty
  if ~iscell(vars{i})
    thisvars = vars{i};
    vars{i} = {};
    vars{i}{1} = thisvars;
    vars{i}{2} = '0';
    varinfo{i}.type = 'numeric';
    % no default argument
  elseif length(vars{i}) == 1
    vars{i}{2} = '0';
    varinfo{i}.type = 'numeric';
  elseif isempty(vars{i}{2})
    vars{i}{2} = '';
    varinfo{i}.type = 'string';
    
    % default arguments have to be strings so they
    % can be put into the text fields properly. here
    % we change them into strings, but remember what
    % type they were
  elseif length(vars{i}) >= 2
    if islogical(vars{i}{2})
      vars{i}{2} = double(vars{i}{2});
    end
    if isnumeric(vars{i}{2})
      % check to see if it is an array
      if isscalar(vars{i}{2})
%         vars{i}{2} = num2str(vars{i}{2});
        varinfo{i}.type = 'numeric';
      else
        varinfo{i}.type = 'array';
      end
    elseif iscell(vars{i}{2})
      varinfo{i}.type = 'popupmenu';
    else
      varinfo{i}.type = 'string';
    end
  end

  % check to see if name is valid
  fixedName = fixBadChars(vars{i}{1});
  if ~strcmp(fixedName,vars{i}{1})
    disp(sprintf('(mrParamsParse) Changing variable name %s to %s',vars{i}{1},fixedName));
    vars{i}{1} = fixedName;
  end
  % set info in varinfo
  varinfo{i}.name = vars{i}{1};
  varinfo{i}.value = vars{i}{2};
  varinfo{i}.description = '';
  
  %--------------------------------------- check for options
  % set defaults
  varinfo{i}.editable = 1;
  varinfo{i}.visible = 1; 
  varinfo{i}.passValue=0; 
  if length(vars{i}) > 2
    %JB: The following loop is the main reason why mrLoadRet is so slow at installing overlays and opening GUIs 
    %(of the mrParams type). mrParamsParse is called by both defaultReconcileParams and mrParamsDialog.
    % I've modified it to make it more efficient 
    %  - by reorganizing the order of the tests
    %  - minimizing calls to string operations
    %  - replacing calls to evalargs by the core operations we need in this case 
    %       (putting a value in a structure field and evaluating a string containing an equal sign)
    % but it's still quite slow...
    skipNext = 0;
    for j = 3:length(vars{i})
      % skip this argument
      if skipNext
        skipNext = 0;
        continue;
      end
      if ~isempty(vars{i}{j})
        equals = strfind(vars{i}{j},'=');
        spaces = strfind(vars{i}{j},' ');
        if isempty(equals)
          if isempty(spaces) && j < (length(vars{i})) %this is a singleword and there is at least one argument after
              varinfo{i}.(vars{i}{j})=vars{i}{j+1}; %so it's the varname/value pair form
              skipNext = 1;
          else % there is at least one space or it's the last argument
            varinfo{i}.description = vars{i}{j}; %this is a description
          end
        elseif isempty(spaces) || spaces(1)>equals(1) %if there is no space or it is after the first equal sign
            %this is a string to evaluate
            value = vars{i}{j}(equals+1:end);
            numericValue = mrStr2num(value);
            if ~isempty(numericValue)
              varinfo{i}.(vars{i}{j}(1:equals-1))= numericValue;
            else
              varinfo{i}.(vars{i}{j}(1:equals-1))=value;
            end
        else %otherwise it's a description
          varinfo{i}.description = vars{i}{j}; %this is a description
        end
      end
    end
    
% % % %     skipNext = 0;
% % % %     for j = 3:length(vars{i})
% % % %       % skip this argument
% % % %       if skipNext
% % % %         skipNext = 0;
% % % %         continue;
% % % %       end
% % % %       % if this looks like a description then save it as a
% % % %       % description, descriptions either have no equal sign
% % % %       % and are not a single word, or have an equal sign but
% % % %       % have spaces before the equal sign
% % % %       if isempty(vars{i}{j}) || ((length(strfind(vars{i}{j},'=')) ~= 1) && (length(strfind(vars{i}{j},' ')) ~= 0)) || ...
% % % %           ~isempty(strfind(vars{i}{j}(1:strfind(vars{i}{j},'=')),' '))
% % % %         varinfo{i}.description = vars{i}{j};
% % % %         % now look for settings that involve the next parameter
% % % %         % i.e. ones that are like 'varname',varvalue. These are
% % % %         % distinugished from comments by the fact that the varname
% % % %         % has no equal sign but is a single word
% % % %       elseif ((length(strfind(vars{i}{j},'=')) ~= 1) && (length(strfind(vars{i}{j},' ')) == 0)) && (j < (length(vars{i})))
% % % %         % we are going to call evalargs but we want the variables
% % % %         % set as a part of gParams (also do it quietly)--> that is
% % % %         % evalargs, will do the parsing of the variable=value strings
% % % %         varargin{1} = 'gVerbose = 0';
% % % %         varargin{2} = sprintf('varinfo{i}.%s',vars{i}{j});
% % % %         varargin{3} = vars{i}{j+1};
% % % %         % set the argument
% % % %         eval(evalargs(varargin));
% % % %         skipNext = 1;
% % % %       else
% % % %         % we are going to call evalargs but we want the variables
% % % %         % set as a part of gParams (also do it quietly)--> that is
% % % %         % evalargs, will do the parsing of the variable=value strings
% % % %         setparam{1} = 'gVerbose = 0';
% % % %         setparam{2} = sprintf('varinfo{i}.%s',vars{i}{j});
% % % %         % set the argument
% % % %         eval(evalargs(setparam));
% % % % 
% % % %       end
% % % %     end
  end
  % Pushbuttons are defined so that you can press them and they set the value
  % of the field (e.g. you can have the button call the rand function and it
  % will set the field to whatever rand returns). You can bypass this behavior
  % by setting passCallbackOutput to 0
  if ~isfield(varinfo{i},'passCallbackOutput')
    if strcmp(varinfo{i}.type,'pushbutton')
      varinfo{i}.passCallbackOutput=1; 
    else
      varinfo{i}.passCallbackOutput=0; 
    end
  end
  %if we need to return the callback output argument, make sure there is at least one 
  if varinfo{i}.passCallbackOutput
    if isfield(varinfo{i},'callback')
      %get the number of output argument returned by the callback function
      if iscell(varinfo{i}.callback) %a callback can be defined either as a cell array whose first cell contains the function handle/name
        nArgout = nargout(varinfo{i}.callback{1});
      else
        nArgout = nargout(varinfo{i}.callback); %or simply as a function handle
      end
      if ~nArgout   %if the callback function returns no output, set this option to 0
        varinfo{i}.passCallbackOutput=0;
      end
    else
      %if there is no callback defined, there is no argument to return
      varinfo{i}.passCallbackOutput=0;
    end
  end
    
  %for popup menus, check the type
  if strcmp(varinfo{i}.type,'popupmenu')       
    varinfo{i}.popuptype = 'string';
    if ~isempty(vars{i}{2})
      if iscell(vars{i}{2}{1}) % if it is a cell (contingent variable, check its first member
        if ischar(vars{i}{2}{1}{1})
          varinfo{i}.popuptype = 'string';
        else
          varinfo{i}.popuptype = 'numeric';
        end
        % see if the default argument is a string
      elseif  ~ischar(vars{i}{2}{1})
        varinfo{i}.popuptype = 'numeric';
      % otherwise numeric
      end
    end
  end

  % make sure type is in lower case
  varinfo{i}.type = lower(varinfo{i}.type);
  % check for minmax violation
  if strcmp(varinfo{i}.type,'numeric') && isfield(varinfo{i},'minmax')
    if vars{i}{2} < varinfo{i}.minmax(1)
      vars{i}{2} = varinfo{i}.minmax(1);
    elseif vars{i}{2} > varinfo{i}.minmax(2)
      vars{i}{2} = varinfo{i}.minmax(2);
    end
  end

  % if it is an array, check to see how bit it is because really big arrays are
  % not able to be displayed
  if strcmp(varinfo{i}.type,'array')
    % first we should not display variables that are larger than 2 dims
    arrayDims = ndims(vars{i}{2});
    if arrayDims > 2
      mrWarnDlg(sprintf('(mrParamsParse) Cannot handle arrays that are greater than 2 dims: %s',vars{i}{1}));
      keyboard
    end
    % now check that the product of dimensions does not exceed
    % the pref: maxArrayWidthForParamsDialog and maxArrayHeightForParamsDialog
    arrayToBig = false;
    arraySize = size(vars{i}{2});
    maxArrayWidthForParamsDialog = mrGetPref('maxArrayWidthForParamsDialog');
    maxArrayHeightForParamsDialog = mrGetPref('maxArrayHeightForParamsDialog');
    if arraySize(2) > maxArrayWidthForParamsDialog
      maxArrayWarning = sprintf('(mrParamsParse) The array for variable ''%s'' is of size [%i %i] which is too big to be displayed by mrParamsDialog (the maximum displayable width is set to %i). You will not be able to set or view this variable. If you would like to try to display anyway, you should set mrSetPref(''maxArrayWidthForParamsDialog'',%i) - or greater. Or you can set the variable by hand after calling mrParamsDialog',vars{i}{1},arraySize(1),arraySize(2),maxArrayWidthForParamsDialog,arraySize(2));
      disp(maxArrayWarning);
      arrayToBig = true;
    end
    if arraySize(1) > maxArrayHeightForParamsDialog
      maxArrayWarning = sprintf('(mrParamsParse) The array for variable ''%s'' is of size [%i %i] which is too big to be displayed by mrParamsDialog (the maximum displayable height is set to %i). You will not be able to set or view this variable. If you would like to try to display anyway, you should set mrSetPref(''maxArrayHeightForParamsDialog'',%i) - or greater. Or you can set the variable by hand after calling mrParamsDialog',vars{i}{1},arraySize(1),arraySize(2),maxArrayHeightForParamsDialog,arraySize(1));
      disp(maxArrayWarning);
      arrayToBig = true;
    end
    % if the array was to big, then we set the variable to display as an uneditable string
    if arrayToBig
      varinfo{i}.tooBigArrayValue = varinfo{i}.value;
      varinfo{i}.value = 'Array too large to display';
      varinfo{i}.editable = 0;
      varinfo{i}.type = 'string';
      varinfo{i}.description = sprintf('%s: %s',varinfo{i}.description,maxArrayWarning);
      % remove the following fields
      fieldsToRemove = {'incdec','incdecType','minmax','passCallbackOutput'};
      for iField = 1:length(fieldsToRemove)
	if isfield(varinfo{i},fieldsToRemove{iField})
	  varinfo{i} = rmfield(varinfo{i},fieldsToRemove{iField});
	end
      end
    end
  end
end

%----------------------------------- now check any contingencies
for i = 1:length(varinfo)
  % groups are handled just like contingent
  if isfield(varinfo{i},'group')
    varinfo{i}.contingent = varinfo{i}.group;
  end
  if any(isfield(varinfo{i},{'contingent','contingentNot'}))
    % get the name of the variable this one is contingent/or contingentNot on
    if isfield(varinfo{i},'contingent')
      contingent = varinfo{i}.contingent;
      contingentNot = false;
    else
      contingent = varinfo{i}.contingentNot;
      contingentNot = true;
    end
    if isempty(contingent),break,end
    % go look for the control variable and set it to
    % have a pointer to this variable
    foundControlVariable = 0;
    for j = 1:length(varinfo)
      if strcmp(contingent,varinfo{j}.name)
        foundControlVariable = 1;
        varinfo{i}.contingentOn = j;
        if ~isfield(varinfo{j},'controls')
          varinfo{j}.controls = i;
	  varinfo{j}.controlsNot = contingentNot;
        else
          varinfo{j}.controls(end+1) = i;
	  varinfo{j}.controlsNot(end+1) = contingentNot;
        end
      end
    end
    % if not found, then complain
    if ~foundControlVariable
      disp(sprintf('Control variable for %s (%s) not found, ignoring contingency',varinfo{i}.name,contingent));
      % otherwise set up the variable to be contingent.
      % that is keep an allValues field of all possible
      % values
    else
      % first deal with popupmenus that have a cell array of values usually
      if strcmp(varinfo{i}.type,'popupmenu')
	% now, if it has a cell array of cell arrays then it wants 
        % to switch between these values
	if iscell(varinfo{i}.value{1})
	  varinfo{i}.allValues = varinfo{i}.value;
	  varinfo{i}.value = varinfo{i}.allValues{1};
	else
	  % otherwise, just a single cell array
	  varinfo{i}.allValues{1} = varinfo{i}.value;
	end
      % other types
      elseif iscell(varinfo{i}.value)
        varinfo{i}.allValues = varinfo{i}.value;
        varinfo{i}.value = varinfo{i}.allValues{1};
      else
        varinfo{i}.allValues{1} = varinfo{i}.value;
      end
      % for numeric values, make sure that the value is set to a number
      if isnumeric(varinfo{i}.value) && ~strcmp(varinfo{i}.type,'checkbox') && ~strcmp(varinfo{i}.type,'array')
        varinfo{i}.value = num2str(varinfo{i}.value);
      end
      % and set a default for the control value. This will get
      % reset later when the dialog starts up to have the
      % value of the control variable
      varinfo{i}.oldControlVal = 0;
    end
  end
end


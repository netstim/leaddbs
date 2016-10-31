% mrParamsSet.m
%
%        $Id$
%      usage: params = mrParamsSet(params)
%         by: justin gardner
%       date: 03/18/08
%    purpose: pass in params and this will reset the fields in an
%             open mrParamsDialog, useful for when you have a
%             non-modal dialog. Note that you do not have to havee
%             all the fields in the params structure, only the ones
%             you want to set.
%             if executeCallback=1, any callback associated with the entry or the dialog box will be executed
%
function retval = mrParamsSet(params,executeCallback)

retval = [];
% check arguments
if ~any(nargin == [1 2])
  help mrParamsSet
  return
end

global gParams;

if ieNotDefined('executeCallback')
  executeCallback = 0;
end

% go through each one of the passed in parameters
paramFields = fieldnames(params);
for fnum = 1:length(paramFields)
  % look for the parameter in varinfo
  match = 0;
  for vnum = 1:length(gParams.varinfo)
    if strcmp(paramFields{fnum},gParams.varinfo{vnum}.name)
      match = vnum;
    end
  end
  % found the match, go set the field
  if match
    % numeric
    if strcmp(gParams.varinfo{match}.type,'numeric')
      % check if this is a group
      if isfield(gParams.varinfo{match},'group')
	groupVar = gParams.varinfo{match}.group;
	% keep values for each one of the group
	for i = 1:length(params.(paramFields{fnum}))
	  gParams.varinfo{match}.allValues{i} = params.(paramFields{fnum})(i);
	end
	% and set the current value in the gui
	currentVal = params.(paramFields{fnum})(params.(groupVar));
	set(gParams.ui.varentry{match},'String',thisNum2str(currentVal));
      
      else
	% simple numeric
	currentVal = params.(paramFields{fnum});
	set(gParams.ui.varentry{match},'String',thisNum2str(currentVal));
      end
      % check for minmax, so that we can turn arrows off and on appropriately
      if isfield(gParams.varinfo{match},'minmax')
	minmax = gParams.varinfo{match}.minmax;
	if isfield(gParams.ui,'incdec')
	  incdecUI = gParams.ui.incdec{match};
	  if ~isempty(incdecUI)
	    if currentVal == min(minmax)
	      set(incdecUI{1},'Enable','off');
	    else
	      set(incdecUI{1},'Enable','on');
	    end
	    if currentVal == max(minmax)
	      set(incdecUI{2},'Enable','off');
	    else
	      set(incdecUI{2},'Enable','on');
	    end
	  end
	end
      end	  
    % checkbox
    elseif strcmp(gParams.varinfo{match}.type,'checkbox')
      if ~isfield(gParams.varinfo{match},'group')
	set(gParams.ui.varentry{match},'Value',params.(paramFields{fnum}));
      else
	groupVar = gParams.varinfo{match}.group;
	% keep values for each one of the group
	for i = 1:length(params.(paramFields{fnum}))
	  gParams.varinfo{match}.allValues{i} = params.(paramFields{fnum})(i);
	end
	% and set the current value in the gui
	currentVal = params.(paramFields{fnum})(params.(groupVar));
	set(gParams.ui.varentry{match},'Value',currentVal);
      end
      % array
    elseif strcmp(gParams.varinfo{match}.type,'array')
      if isequal(size(gParams.varinfo{match}.value),size(params.(paramFields{fnum})))
	for matrixX = 1:size(params.(paramFields{fnum}),1)
	  for matrixY = 1:size(params.(paramFields{fnum}),2)
      set(gParams.ui.varentry{match}(matrixX,matrixY),'String',thisNum2str(params.(paramFields{fnum})(matrixX,matrixY)));
	  end
	end
      else
	disp(sprintf('(mrParamsSet) Array size of variable %s does not match',paramsFields{fnum}));
      end
    % popupmenu
    elseif strcmp(gParams.varinfo{match}.type,'popupmenu')
      % check if this is a group
      if isfield(gParams.varinfo{match},'group')
	% Don't do anything for a group parameter--fix later if you need this
      else
	value = find(strcmp(params.(paramFields{fnum}),gParams.varinfo{match}.value));
	if ~isempty(value)
	  set(gParams.ui.varentry{match},'Value',value);
	else
	  disp(sprintf('(mrParamsSet) %s is not a valid option for variable %s',params.(paramFields{fnum}),paramFields{fnum}));
	end
      end
    % push button. 
    elseif strcmp(gParams.varinfo{match}.type,'pushbutton')
      gParams.varinfo{match}.value = params.(paramFields{fnum});
    % string
    elseif strcmp(gParams.varinfo{match}.type,'string')
      % check if we are group'd
      if isfield(gParams.varinfo{match},'group')
	groupVar = gParams.varinfo{match}.group;
	% keep values for each one of the group
	for i = 1:length(params.(paramFields{fnum}))
	  gParams.varinfo{match}.allValues{i} = params.(paramFields{fnum}){i};
	end
	stringValue = params.(paramFields{fnum}){params.(groupVar)};
      % check if we are contingentOn
      elseif isfield(gParams.varinfo{match},'contingentOn')
	contingentValue = params.(gParams.varinfo{gParams.varinfo{match}.contingentOn}.name);
	stringValue = gParams.varinfo{match}.allValues{contingentValue};
      else
	stringValue = params.(paramFields{fnum});
      end
      set(gParams.ui.varentry{match},'String',stringValue);
    % unimplemented type
    else
      disp(sprintf('(mrParamsSet) Setting of type %s not implemented yet',gParams.varinfo{match}.type));
    end
    
    if executeCallback 
      % update params
      if isfield(gParams, 'callback')
        if ~isempty(gParams.callback)
          gParams.params = mrParamsGet(gParams.vars);
          if isfield(gParams,'callbackArg')
            feval(gParams.callback,gParams.params,gParams.callbackArg);
          else
            feval(gParams.callback,gParams.params);
          end
        end
      end
      % handle callbacks for non-push buttons
      if ~fieldIsNotDefined(gParams.varinfo{match},'callback')...
          && ~gParams.varinfo{match}.passCallbackOutput && ~strcmp(gParams.varinfo{match}.type,'pushbutton')
        callbackArgs={};
        if isfield(gParams.varinfo{match},'callbackArg')
          callbackArgs{end+1} = gParams.varinfo{match}.callbackArg;
        end
        callbackArgs{end+1} = mrParamsGet(gParams.vars);
        callbackEval(gParams.varinfo{match}.callback,callbackArgs);
      end
    end

  else
    if ~strcmp(paramFields{fnum},'paramInfo')
      disp(sprintf('(mrParamsSet) Could not find var %s',paramFields{fnum}));
    end
  end
end


if nargout == 1
  retval = mrParamsGet;
end


%modified num2str to increase the number of decimals for reals
function value = thisNum2str(value)

  if rem(value,1)~=0
    value = num2str(value,'%.6f');
  else
    value = num2str(value);
  end

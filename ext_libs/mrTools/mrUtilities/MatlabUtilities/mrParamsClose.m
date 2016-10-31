% mrParamsClose.m
%
%      usage: mrParamsClose()
%         by: justin gardner
%       date: 09/04/11
%    purpose: Closes an open mrParamsDialog without any arguments. With argument set to true, closes
%             as if you hit ok.
%
%             mrParamsClose(true);
%
function retval = mrParamsClose(tf)

% check arguments
if ~any(nargin == [0 1])
  help mrParamsClose
  return
end

global gParams;

% act as if ok was hit
if (nargin >= 1)
  if tf
    % try to hit the ok button
    if isfield(gParams,'okButton')
      if ~isempty(get(gParams.okButton,'Callback'))
	feval(get(gParams.okButton,'Callback'));
	return
      end
    end
    % if we got here, we could not call okButtoin
    disp(sprintf('(mrParamsClose) Called to simulate ok button close, but could not find okButton handler'));
    return
  else
    % try to hit the cancel button
    if isfield(gParams,'cancelButton')
      if ~isempty(get(gParams.cancelButton,'Callback'))
	feval(get(gParams.cancelButton,'Callback'));
	return
      end
    end
    % if we got here, we could not call okButtoin
    disp(sprintf('(mrParamsClose) Called to simulate cancel button close, but could not find cancelButton handler'));
    return
  end
end

% just close the figure without doing anything (no ok/cancel calling)
if isfield(gParams,'fignum')
  for i = 1:length(gParams.fignum)
    if ishandle(gParams.fignum(i))
      close(gParams.fignum(i));
    end
    if ishandle(gParams.fignum(i))
      delete(gParams.fignum(i));
    end
  end
end


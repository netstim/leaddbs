% mlrLoadLastView.m
%
%        $Id:$ 
%      usage: mlrLoadLastView()
%         by: justin gardner
%       date: 04/02/15
%    purpose: This function became necessary since in version 8.2, mathworks 
%             inexplicably busted figure handles. Now they are handles and
%             if you load an old variable with a handle, it causes massive
%             convulsions in matlab. It also tries to pop up a figure - which
%             is not what we want - eeks. Anyway, the solution is to load
%             part of the mrLastView variable (viewSettings) and check
%             its version - if it is an old style version then warn the user
%             that we can't load
%
function [v, viewSettings]  = mlrLoadLastView(filename)

v = [];
viewSettings = [];

% check arguments
if ~any(nargin == [0 1])
  help mlrLoadLastView
  return
end

if nargin == 0;
  filename = 'mrLastView.mat';
end
filename = setext(filename,'mat');

if ~isfile(filename);
  mrWarnDlg('(mlrLoadLastView) Could not find %s',filename);
  return
end

% new version of matlab needs to check if we are trying to load an old file
if ~verLessThan('matlab','8.1')
  if verLessThan('matlab','8.4') % for versions before 8.5, if mrLastView.mat has previously been saved with version >=8.5,
    % loading it issues a lot of warnings (even loading just viewSettings when it contains panel handles)
    % but everything seems to be loading fine, so disabling warnings
    warning('off','MATLAB:class:mustReturnObject'); 
  end
  % load the viewSettings part and look to see if it is a new file
  check = load(filename,'viewSettings');
  if isfield(check,'viewSettings') && isfield(check.viewSettings,'version') && (check.viewSettings.version>=2.0)
    % then we are ok, load the view part
    l = load(filename,'view');
    % return them both
    if nargout == 1
      % return as single argument
      v = l;
    else
      % or as two
      v = l.view;
      viewSettings = check.viewSettings;
    end
  else
    if verLessThan('matlab','8.4')
      % it can load, but will give lots of warnings, so tell user what is going on 
      l = load(filename);
      %JB: in fact I think the warnings appear only if mrLastView has been saved 
      % with matlab >=8.5 and is now being loaded with version <8.5. I already turned
      % them off above because warnings are issued also when loading only viewSettings
%       mrWarnDlg(sprintf('(mlrLoadLastView) The mrLastView found: %s is from an older version of matlab, you will likely see a bunch of weird warnings here, but ignore them - they have to do with the latest matlab not having the ability to load old mat files that had figure handles in them. Send complaints to Mathworks!',filename));
      % this causes lots of weird warnings, but doesn't seem to crash
      if isfield(l,'view') && isfield(l.view,'figure')
	if ishandle(l.view.figure)
	  close(l.view.figure);
	end
	l.view.figure = [];
      end
      % return as either one or two arguments
      if nargout == 1
	% return as single argument
	v = l;
      else
	% or as two
	v = l.view;
	viewSettings = l.viewSettings;
      end
    else
      % serious problems occur in 8.5
      mrWarnDlg(sprintf('(mlrLoadLastView) The mrLastView found: %s is from an older version of matlab which allowed saving figure handles. The geniuses at Mathworks have busted that, so loading this file will no longer work. Moving mrLastView.mat to mrLastView.mat.old You will lose any rois that were loaded but not saved and mrLoadRet will start up without bases and analyses loaded. If you really need what was in the viewer we suggest running on an earlier version of matlab - you just then need to copy mrLastView.mat.old back to mrLastView.mat, open mrLoadRet and then quit - this will save the file back w/out the offending figure handles. Send complaints to Mathworks!',filename),'Yes');
      movefile(filename,sprintf('%s.old',filename));
      return
    end
  end
  if verLessThan('matlab','8.4')
    warning('on','MATLAB:class:mustReturnObject');  % turning warning back on
  end
else
  % otherwise just load
  l = load('mrLastView');
  if nargout == 1
    v = l;
  else
    v = l.view;
    viewSettings = l.viewSettings;
  end
end



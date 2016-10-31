% mlrSmartfig.m
%
%        $Id:$ 
%      usage: mlrSmartfig(figname,<reuse>)
%         by: justin gardner
%       date: 09/06/11
%    purpose: puts up a figure at the same position and size of where it was
%             last opened. That is, if you call
% 
%             mlrSmartfig('testSmartfig');
%
%             move and resize the figure and close, the next time you call the same
%             it will put the figure up in the same place.
%
%             If reuse is set to 1 then if there is already an open figure
%             with the same name then it will change current figure focus
%             to that open figure instead of opening a new figure.
%
function f = mlrSmartfig(figname,reuse)

% check arguments
if ~any(nargin == [1 2])
  help mlrSmartfig
  return
end


% global that stores all figure handles
global gMLRSmartfig;

% if this is a number then it means to close
if isnumeric(figname)
  entryNum = figname;
  if ~isempty(gMLRSmartfig) && ~isequal(gMLRSmartfig(entryNum).fignum,-1)
    % remember position
    mrSetFigLoc(gMLRSmartfig(entryNum).figname,get(gMLRSmartfig(entryNum).fignum,'Position'));
    % close the figure
    delete(gMLRSmartfig(entryNum).fignum);
    % clear out its entry
    gMLRSmartfig(entryNum).fignum = -1;
    gMLRSmartfig(entryNum).figname = '';
  end
  return
end

% remove any problem characters
figname = fixBadChars(figname);

% get the figloc
figloc = mrGetFigLoc(figname);

% if we can reuse, then look for the a match
if ~ieNotDefined('reuse') && ~isequal(reuse,0) && ~isempty(gMLRSmartfig) && any(strcmp({gMLRSmartfig.figname},figname))
  f = gMLRSmartfig(first(find(strcmp({gMLRSmartfig.figname},figname)))).fignum;
  figure(f);
else
  % open figure
  f = figure;
  % find an entry location
  if isempty(gMLRSmartfig) || ~any([gMLRSmartfig.fignum] == -1)
    entryNum = length(gMLRSmartfig)+1;
  else
    entryNum = first(find([gMLRSmartfig.fignum] == -1));
  end
  % remember figure number
  gMLRSmartfig(entryNum).fignum = mlrGetFignum(f);
  % remember figure name
  gMLRSmartfig(entryNum).figname = figname;
  % set the position
  if ~isempty(figloc),set(f,'Position',figloc);end
  % set the close callback
  set(f,'CloseRequestFcn',sprintf('mlrSmartfig(%i)',entryNum));
end






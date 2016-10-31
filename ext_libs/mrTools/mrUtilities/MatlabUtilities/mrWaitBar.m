function h = mrWaitBar(x,t,newMessage)
%
% h = mrWaitBar([x],[title or handle])
%
% Calls Matlab's waitbar or disp depending on verbose preference.
% x: percentage of waiting done
% tile or handle: when first opening a waitbar this should be a string
%    which is used as the title of the window. Thereafter, it should be the
%    handle to the waitbar.
%
% To set the 'verbose' preference:
%    mrSetPref('verbose','Yes');
%    mrSetPref('verbose','No');
%
% djh, 5/2005
%
% djh, 5/2007, modified to use mrGetPref instead of Matlab's getpref

if ~exist('x','var')
    x=0;
end
if ~exist('t','var')
    t = 'Please wait';
end

% update of wait bar
if ishandle(t)
  if ieNotDefined('newMessage')
    waitbar(x,t);
  else
    waitbar(x,t,newMessage);
  end
  drawnow;
elseif isfield(t,'disppercent')
  if ieNotDefined('newMessage')
    disppercent(x);
  else
    disppercent(x,newMessage);
  end
    
% initial call
elseif ischar(t)
  % check the verbose preference
  verbose = mrGetPref('verbose');
  if strcmp(verbose,'Yes')
    % if verbose, make a window wait bar
    if length(t)<=35
      if ieNotDefined('newMessage')
         h = waitbar(x,'','name',t);
      else
         h = waitbar(x,newMessage,'name',t);
      end
    else
      if ieNotDefined('newMessage')
        h = waitbar(x,t,'name','mrWaitBar');
      else
        h = waitbar(x,{t,newMessage},'name','mrWaitBar');
      end
    end
    drawnow;
  else
    % otherwise use disppercent
    if ieNotDefined('newMessage')
      disppercent(-inf,t);
    else
      disppercent(-inf,[t '; ' newMessage]);
    end 
    h.disppercent = 1;
  end
end
return

% Test/debug
mrSetPref('verbose','Yes');
mrSetPref('verbose','No');

startTime = mglGetSecs;
h = mrWaitBar(0,'Test. Please wait...');
for i=1:100
    mrWaitBar(i/100,h);
end
mrCloseDlg(h);
mglGetSecs(startTime)

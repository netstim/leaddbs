% mrLoadRet, version 4.0 
%
% AUTHOR:   DJH
% PURPOSE:  Everything.
% DATE:     June, 2004
%
% usually, just run as
% mrLoadRet
%
% to run mrLoadRet without loading mrLastView (i.e. this
% will run mrLoadRet fresh without remembering your settings
% and keeping your loaded analyses, base volumes and ROIs from
% your last session).
% mrLoadRet([]);
%
% to load a different mrLastView
% mrLoadRet('mrLastView2');
%
function [v]= mrLoadRet(mrLastView)

% default to 
if nargin < 1
  mrLastView = 'mrLastView.mat';
end
if ~isfile('mrSession.mat')
  disp('(mrLoadRet) No mrSession.mat found in current directory');
  return
end

mlrPath mrTools;

% Define and initialize global variable MLR.
mrGlobals

% check to make sure we are not being run from another data directory
% since mrLoadRet is designed only to be run on one data directory at a time
if isempty(strfind(pwd,viewGet([],'homeDir')))
  % check to see if there is a session mismatch
  [thisSession thisGroups] = loadSession(pwd);
  if ~isempty(thisSession) && (~isequal(thisSession,MLR.session) || ~isequal(thisGroups,MLR.groups))
    disp(sprintf('(mrLoadRet) Current path: %s does not match',pwd));
    disp(sprintf('(mrLoadRet) homeDir: %s in MLR global',viewGet([],'homeDir')));
    disp(sprintf('(mrLoadRet) If you are trying to run two mrLoadRet sessions on different'));
    disp(sprintf('(mrLoadRet) datasets, you should run two separate matlab processes instead'));
    disp(sprintf('(mrLoadRet) Otherwise type mrQuit to quit your current session and try running mrLoadRet again'));
    return
  end
end

% Open inplane window
v = mrOpenWindow('Volume',mrLastView);

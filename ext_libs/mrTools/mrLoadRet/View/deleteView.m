function deleteView(view)
%
% deleteView(view)
%
% Removes the view from MLR.views
% view can be a view or a viewNum
%
% djh, 6/2004

if nargin ~= 1
  help deleteView;
  return
end

% get mrGlobals
mrGlobals

% flag to clear input variable
clearViewInCallerSpace = 0;

% passed in number, then it is a viewNum, not a view
if isnumeric(view)
  if ~isempty(viewGet([],'view',view))
    viewNum = view;
  else
    % not a valid viewNum, print error message and return
    if isempty(view)
      mrWarnDlg(sprintf('(%s) Empty view.',mfilename));
    else
      mrWarnDlg(sprintf('(%s) View number %i does not exist.',mfilename,view));
    end
    return
  end
% make sure it is a view
elseif isview(view)
  viewNum = view.viewNum;
  % if the input name is a variable name, flag it to be cleared int
  % the calling space
  if ~isempty(inputname(1))
    clearViewInCallerSpace = 1;
  end
else
  return;
end

% Delete it
MLR.views{viewNum} = [];
MLR.caches{viewNum} = [];

% clear it in the calling space
if clearViewInCallerSpace
  assignin('caller',inputname(1),[]);
end

% check to see if there are any valid views, if not, clear MLR
for i = 1:length(MLR.views)
  if ~isempty(MLR.views{i}),return,end
end
clear global MLR;

% Remove it and update viewNum for all of the remaining views
% (busted because doens't update local variables) 
%
% numviews = length(MLR.views);
% MLR.views = {MLR.views{find(viewNum ~= [1:numviews])}};
% for v = 1:length(MLR.views)
%   MLR.views{v}.viewNum = v;
% 	mlrGuiSet(v,'viewNum',v);
% end


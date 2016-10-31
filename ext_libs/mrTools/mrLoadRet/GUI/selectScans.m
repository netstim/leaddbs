function scanList = selectScans(view,title,scanNums)
% scanList = selectScans(view,[title],[scanNums]);
%
%   this function is deprecated, use scanList = selectInList(thisView,'scans',title,preselected)
%
%   Gather a list of scans available in Inplane/TSeries
%   and query the user for a sub-selection.
%
%   An alternate functionchooseScans uses numScans(view)
%   to determine the number of scans to choose from
%   Use selectScans if you will be analyzing the tSeries. 
%   Use chooseScans, if your code does not depend on the 
%   presence/absence of the tSeries files.
%
%   If scanNums is present, then will only allow the user
%   to select a scan from the scanNums list
%
% Output:
%  scanList: list of selected scans.
%
% 4/16/99  dbr Initial code
% 3/30/2001, djh, added optional title string
% 11/9/06 jlg mrLoadRet 4 conversion
%
% $Id$	

if ieNotDefined('title')
  title = 'Choose scans';
end

% get number of scans
nScans = viewGet(view,'nScans');

% get scanNums
if ieNotDefined('scanNums')
  scanNums = 1:nScans;
end

%Check for zero:
if nScans == 0
  mrErrorDlg('No scans found!');
  return
end

for i = 1:length(scanNums)
  scanNames{i} = sprintf('%i:%s (%s)',scanNums(i),viewGet(view,'description',scanNums(i)),viewGet(view,'tSeriesFile',scanNums(i)));
end

% Which scans to analyze?
iSel = buttondlg(title, scanNames);
scanList = scanNums(find(iSel));

return;

% deleteScans.m
%
%        $Id$
%      usage: deleteScans(<v>,<scanList>,<groupName/Num>)
%         by: justin gardner
%       date: 09/19/07
%    purpose: delete scans
%
function v = deleteScans(v,scanList,group)

% check arguments
if ~any(nargin == [0 1 2 3])
  help deleteScans
  return
end

% get a view if not passed in one.
if ieNotDefined('v')
  % get a new view
  v = newView;
  % remember to delete it when we exit this funciton
  % unless the user wants it
  if nargout ==0
    deleteViewWhenDone = 1;
  end
  % choose a group
  groupNames = viewGet(v,'groupNames');
  params = mrParamsDialog({{'groupName',groupNames,'type=popupmenu','Choose a group from which to delete scans'}},'Choose a group from which to delete scans');
  if isempty(params)
    if deleteViewWhenDone
      deleteView(v);
    end   
    return
  end
  v = viewSet(v,'curGroup',params.groupName);
else
  deleteViewWhenDone = 0;
end

% set group if passed in
if ~ieNotDefined('group')
  v = viewSet(v,'curGroup',group);
end

% get a scan list if it has not been passed in
if ieNotDefined('scanList')
  scanList = selectInList(v,'scans');
end

% check the scan list
badScanNums = setdiff(scanList,1:viewGet(v,'nScans'));
if ~isempty(badScanNums)
  disp(sprintf('(mrDeleteScans) No scans %s',num2str(badScanNums)));
  return
end

for iScan = 1:length(scanList)
    % first get time series name for each one of these scans
    % since as we delete them, then numbers stop making sense
    tSeriesFile{iScan} = viewGet(v,'tSeriesFile',scanList(iScan));
end
% now go through and delete
for iScan = 1:length(scanList)
    % get the scan number
    scanNum = viewGet(v,'scanNum',tSeriesFile{iScan});
    if ~isempty(scanNum)
        scanNum = scanNum(1);
        v = viewSet(v,'deleteScan',scanNum);
        disp(sprintf('Scan for file %s deleted.',tSeriesFile{iScan}));
    else
        disp(sprintf('(mrDeleteScans) Could not delete scan for file %s',tSeriesFile{iScan}));
    end
    if ~isempty(v.figure)
      v = viewSet(v,'overlayCache','init'); %need to empty the cache because scan numbers have changed
      refreshMLRDisplay(v.viewNum);
    end
end
if ~isempty(scanList)
    disp(sprintf('To remove the nifti files for these deleted scans run mrCleanDir'));
end

if deleteViewWhenDone
  deleteView(v);
end
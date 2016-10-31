% defaultReconcileParams.m
%
%        $Id$
%      usage: defaultReconcileParams(groupName,params,data)
%         by: justin gardner
%       date: 03/13/07
%    purpose: A default params reconcile. This does not necessarily
%             have to be called for params created by mrParamsDialog,
%             but if it is, then there will be a paramInfo field (for format
%             see wiki), and variables will be validated against that structure
%
%             The main purpose of this function is to match
%             timeseries filenames with scan numbers
%
function [params data] = defaultReconcileParams(groupName,params,data)

% check arguments
if ~any(nargin == [2 3])
  help defaultReconcileParams
  return
end

if ieNotDefined('data'),data=[];,end

% get group name
if ieNotDefined('groupName')
  if isfield(params,'groupName')
    groupName = params.groupName;
  else
    groupName = '';
  end
end
% set group number
if ~ieNotDefined('groupName')
  groupNum = viewGet([],'groupNum',groupName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if params is a cell array, then call defaultReconcileParams on each
% element of the cell array one at a time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iscell(params)
  newData = cell(1,viewGet([],'numScans',groupNum));
  newParams = cell(1,viewGet([],'numScans',groupNum));
  for i = 1:length(params)
    if ~isempty(params{i})
      if length(data) >= i
	[thisParams thisData] = defaultReconcileParams(groupName,params{i},data);
	if isfield(thisParams,'scanNum')
	  if ~isempty(thisParams.scanNum)
	    newParams{thisParams.scanNum} = thisParams;
	    newData{thisParams.scanNum} = thisData{thisParams.scanNum};
	  end
	elseif ~isempty(thisParams)
	  newParams{i} = thisParams;
	  newData{i} = thisData{i};
	end
      else
	params{i} = defaultReconcileParams(groupName,params{i});
      end
    end
  end
  data = newData;
  params = newParams;
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if this has the field paramInfo then it comes from'
% mrDefaultParamsGUI and we can check the parameters 
% for consistency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(mrGetPref('checkParamsConsistency'),'Yes') %JB: this is annoyingly slow, too complicated to fix and probably unnecessary, so I added an option to slip
  params = checkParams(params);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% look for a description, and see if it has something like [x...x], that should be replaced by the scan numbers selected in scanList and groupName 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params = fixDescription(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Next we do the meat of this program, we match scanNumbers
% with tseriesFilenames - so that the data will stay with
% the actual timeseries that it was generated with
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get scan numbers
if isfield(params,'scanList')
  scanListName = 'scanList';
elseif isfield(params,'scanNum')
  scanListName = 'scanNum';
elseif isfield(params,'tseriesFile')
  scanListName = 'scanList';
  params.scanList = 1:length(params.tseriesFile);
else
  scanListName = '';
end

% get a list of fields that are associated with scan numbers
% these will get matched to the new scan number that we get
% below when we reconcile against the filename
scanFields = {};
if isfield(params,'scanParams')
  scanFields{end+1} = 'scanParams';
end

if ~isempty(scanListName) && ~isempty(params.(scanListName))
  % if we don't have a tseriesFile field then generate it
  % this will usually happen the first time this function
  % is called on some params
  if ~isfield(params,'tseriesFile')
    % get the tseries name
    for iscan = 1:length(params.(scanListName))
      params.tseriesFile{iscan} = viewGet([],'tseriesFile',params.(scanListName)(iscan),groupNum);
    end
  % the second time through, we will have tseriesFile
  % so now we have to get the correct scan numbers for those
  % tseriesFiles.
  else
    % starting values
    newScanNums = [];
    newTSeriesFile = {};
    newdata = cell(1,viewGet([],'numScans',groupNum));
    for iFields = 1:length(scanFields)
      newScanFields.(scanFields{iFields}) = cell(1,viewGet([],'numScans',groupNum));
    end
    % see if the tseriesFile name matches
    for tnum = 1:length(params.tseriesFile)
      % get the scan number associated with the tseriesFile
      thisScanNum = viewGet([],'scanNum',params.tseriesFile{tnum},groupNum);
      % check for empty
      if isempty(thisScanNum)
	disp(sprintf('(defaultReconcileParams) Ignoring scan with tseriesFilename: %s because it no longer exists',params.tseriesFile{tnum}));
      % Not sure why the below happens--essentially there is no scanNum associated
      % with the tseries filename. Since there isn't any associated info
      % for that timeseries filename, we need to drop it as well (otherwise
      % the code will crash below because it can't find the correct
      % scanparams for the tseriesfile.
      elseif length(params.(scanListName)) < tnum
	disp(sprintf('(defaultReconcileParms) No associated scanParams/scanNum for analysis on %s',params.tseriesFile{tnum}));
      else
	% found it, set the scan number
	newScanNums(end+1) = thisScanNum;
	% keep the fields that go with each scan
	for iFields = 1:length(scanFields)
	  newScanFields.(scanFields{iFields}){thisScanNum} = ...
	      params.(scanFields{iFields}){params.(scanListName)(tnum)};
	end
	% keep the tseries filename
	newTSeriesFile{end+1} = params.tseriesFile{tnum};
	% and if there is data, fix the data.
	if ~isempty(data)
	  newdata{thisScanNum} = data{params.(scanListName)(tnum)};
	  % fix filename and path fields
	  if isfield(newdata{thisScanNum},'filename')
	    newdata{thisScanNum}.filename = params.tseriesFile{tnum};
	    newdata{thisScanNum}.filepath = fullfile(fileparts(newdata{thisScanNum}.filepath),params.tseriesFile{tnum});
	  end
	end
      end
    end
    % reset the scan field
    params.(scanListName) = newScanNums;
    params.tseriesFile = newTSeriesFile;
    % if there are no valid scan numbers then warn and return empty data
    if isempty(newScanNums)
      disp(sprintf('(defaultReconcileParams) All the previous scans that had been analyzed in this analysis no longer exist.'));
      data = {};
    % and reset the data
    else
      if ~isempty(data)
	data = newdata;
      end
    end
    % and reset all the scan fields
    for iFields = 1:length(scanFields)
      if isempty(newScanNums)
	params.(scanFields{iFields}) = {};
      else
	params.(scanFields{iFields}) = newScanFields.(scanFields{iFields});
      end
    end
  end  
end

% if there is a field called scanParams then reconcile those as well
if isfield(params,'scanParams') && iscell(params.scanParams) && ~isempty(params.scanParams)
  defaultReconcileParams(groupName,params.scanParams);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to check params generated by mrParamsDialog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = checkParams(params)

if isfield(params,'paramInfo')
  [vars varinfo] = mrParamsParse(params.paramInfo);
  for i = 1:length(varinfo)
    useDefault = 0;
    % make sure the field exists
    if ~isfield(params,varinfo{i}.name)
      useDefault = 1;
    % check if it is good
    else
      % is it numeric
      if strcmp(lower(varinfo{i}.type),'numeric')
	if ~isnumeric(params.(varinfo{i}.name))
	  useDefault = 1;
	end
	% also check if there is a minmax
	if isfield(varinfo{i},'minmax')
	  if params.(varinfo{i}.name) < varinfo{i}.minmax(1)
	    disp(sprintf('(defaultReconcileParams) %s out of range. Setting to min=%f',varinfo{i}.name,varinfo{i}.minmax(1)));
	    params.(varinfo{i}.name) = varinfo{i}.minmax(1);
	  elseif params.(varinfo{i}.name) > varinfo{i}.minmax(2)
	    disp(sprintf('(defaultReconcileParams) %s out of range. Setting to max=%f',varinfo{i}.name,varinfo{i}.minmax(2)));
	    params.(varinfo{i}.name) = varinfo{i}.minmax(2);
	  end
	end
      % or a checkbox
      elseif strcmp(lower(varinfo{i}.type),'checkbox')
	if isempty(params.(varinfo{i}.name)) || ~any(params.(varinfo{i}.name) == [0 1])
	  useDefault = 1;
	end
      end
    end      
    % see if we have to switch it to default
    if useDefault
      % make sure it is not a contingent value that has been shut
      % off, first get value it is contingent on
      if isfield(varinfo{i},'contingent')
	if isfield(params,varinfo{varinfo{i}.contingentOn}.name)
	  contingentValue = params.(varinfo{varinfo{i}.contingentOn}.name);
	else
	  contingentValue = varinfo{varinfo{i}.contingentOn}.value;
	  if iscell(contingentValue)
	    contingentValue = contingentValue{1};
	  end
	end
	if isstr(contingentValue),contingentValue = str2num(contingentValue);,end
	% if it has been shut down, give the parameter an empty
        % value and continue on
	if contingentValue==0
	  params.(varinfo{i}.name) = [];
	  continue
	end
      end
      % otherwise set the parameter to the default value
      if ~iscell(varinfo{i}.value)
	params.(varinfo{i}.name) = varinfo{i}.value;
      else
	params.(varinfo{i}.name) = varinfo{i}.value{1};
      end
      if isnumeric(params.(varinfo{i}.name))
	disp(sprintf('(defaultReconcileParams) Using default value for %s (%s)',varinfo{i}.name,num2str(params.(varinfo{i}.name))));
      else
	disp(sprintf('(defaultReconcileParams) Using default value for %s (%s)',varinfo{i}.name,params.(varinfo{i}.name)));
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to fix description field to contain
% info about group and scans
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = fixDescription(params)

if isfield(params,'description') && ~isempty(strfind(params.description,'[x...x]'))
  swaploc = strfind(params.description,'[x...x]');
  swaploc = swaploc(1);
  % get the group name
  if isfield(params,'groupName')
    groupName = params.groupName;
  else
    groupName = '';
  end
  % scan numbers
  if isfield(params,'scanList')
    scanNames = num2str(params.scanList);
  elseif isfield(params,'scanNum')
    scanNames = num2str(params.scanNum);
  else
    scanNames = '';
  end
  % now make the description string
  if ~isempty(groupName)
    params.description = sprintf('%s%s:%s%s',params.description(1:swaploc-1),groupName,scanNames,params.description(swaploc+7:end));
  else
    params.description = sprintf('%s%s%s',params.description(1:swaploc-1),scanNames,params.description(swaploc+7:end));
  end
end

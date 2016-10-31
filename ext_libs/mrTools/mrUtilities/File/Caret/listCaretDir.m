% listCaretFiles.m
%
%        $Id$ 
%      usage: listCaretFiles(<dirname=dirname>,<numNodes=numNodes>)
%         by: justin gardner
%       date: 12/03/09
%    purpose: lists the caret files in a directory, with node information
%      usage: listCaretFiles
%             To only list files with a certain number of nodes, e.g. 73730:
%             listCaretFiles('numNodes=73730');
%             Only display files with matching string
%             listCaretFiles('matchStr=Human');
%
function [outTopoFiles outCoordFiles]  = listCaretDir(varargin)

% check arguments
if ~any(nargin == [0 1 2 3 4])
  help listCaretFiles
  outCoordFiles = [];
  outTopoFiles = [];
  return
end

dirname=[];
numNodes=[];
matchStr=[];
noDisplay=[];
getArgs(varargin,{'dirname=.','numNodes=[]','matchStr=[]','noDisplay=0'});

% look for all topo files
topoFiles = dir(fullfile(dirname,'*.topo'));

if ~isempty(topoFiles)
  if ~noDisplay,disp(sprintf('================ TOPO FILES ================'));end
  % gather info about topo files
  for i = 1:length(topoFiles)
    % default sort index
    topoSortIndex(i) = inf;
    topoDateIndex(i) = inf;
    % load the topo file
    if isempty(matchStr) || ~isempty(strfind(topoFiles(i).name,matchStr))
      topoFile = openCaretFile(fullfile(dirname,topoFiles(i).name));
      if ~isempty(topoFile)
	% get tiles and nodes
	topoFiles(i).numTiles = topoFile.num_tiles;
	topoFiles(i).numNodes = max(topoFile.data(:));
	% get other info
	topoFiles(i).comment = topoFile.comment;
	topoFiles(i).date = topoFile.date;
	% get arrays of numNodes and datenum for sorting
	topoSortIndex(i) = topoFiles(i).numNodes;
	topoDateIndex(i) = topoFiles(i).datenum;
      end
    end
  end

  % sort by time and then by numNodes
  [topoDateIndex topoDateOrder] = sort(topoDateIndex);
  [topoSortIndex topoSortOrder] = sort(topoSortIndex(topoDateOrder));
  % now print out all topo files
  outIndex = 1;
  for i = 1:length(topoSortOrder)
    thisTopo = topoFiles(topoDateOrder(topoSortOrder(i)));
    if isempty(numNodes) || isequal(numNodes,thisTopo.numNodes)
      if isempty(matchStr) || ~isempty(strfind(thisTopo.name,matchStr))
	if ~noDisplay,disp(sprintf('%i: (%i nodes, %i tiles) %s [%s]',i,thisTopo.numNodes,thisTopo.numTiles,thisTopo.name,datestr(thisTopo.datenum)));end
	% save for output structure
	outTopoFiles(outIndex) = thisTopo;
	outIndex = outIndex + 1;
      end
    end
  end
end

% look for all coord files
coordFiles = dir(fullfile(dirname,'*.coord'));

if ~isempty(coordFiles)
  if ~noDisplay,disp(sprintf('================ COORD FILES ================'));end
  % gather info about coord files
  for i = 1:length(coordFiles)
    % default sort order
    coordSortIndex(i) = inf;
    coordDateIndex(i) = inf;
    % load the coord file
    if isempty(matchStr) || ~isempty(strfind(coordFiles(i).name,matchStr))
      coordFile = openCaretFile(fullfile(dirname,coordFiles(i).name));
      if ~isempty(coordFile)
	% get numNodes
	coordFiles(i).numNodes = coordFile.num_nodes;
	% get other info
	coordFiles(i).comment = coordFile.comment;
	coordFiles(i).date = coordFile.date;
	% get arrays of numNodes and datenum for sorting
	coordSortIndex(i) = coordFiles(i).numNodes;
	coordDateIndex(i) = coordFiles(i).datenum;
      end
    end
  end

  % sort by time and then by numNodes
  [coordDateIndex coordDateOrder] = sort(coordDateIndex);
  [coordSortIndex coordSortOrder] = sort(coordSortIndex(coordDateOrder));

  % print out all the coord files
  outIndex = 1;
  for i = 1:length(coordFiles)
    thisCoord = coordFiles(coordDateOrder(coordSortOrder(i)));
    if isempty(numNodes) || isequal(numNodes,thisCoord.numNodes)
      if isempty(matchStr) || ~isempty(strfind(thisCoord.name,matchStr))
	if ~noDisplay,disp(sprintf('%i: (%i nodes) %s [%s]',i,thisCoord.numNodes,thisCoord.name,datestr(thisCoord.datenum)));end
	% save for output structure
	outCoordFiles(outIndex) = thisCoord;
	outIndex = outIndex + 1;
      end
    end
  end
end

if ieNotDefined('outCoordFiles'),outCoordFiles = [];end
if ieNotDefined('outTopoFiles'),outTopoFiles = [];end

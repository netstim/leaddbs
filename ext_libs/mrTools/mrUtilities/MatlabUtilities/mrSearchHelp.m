% mrSearchHelp - grep mrLoadRet code for occurrences of string
%
%      usage: [  ] = mrSearchHelp( expression )
%         by: denis schluppeck
%       date: 2008-07-31
%        $Id$:
%     inputs: expression ; string or cell array of strings
%             options ; command line options to grep
%    outputs: 
%
%    purpose: it's sometimes helpful to search through all the code in
%    mrLoadRet  and  mrAlign to  look  for  occurrences of  functions,
%    strings, etc. This function  makes use of a matlab implementation
%    of the unix utility GREP.
%
%    help grep % for further info on regexp and options
%
%        e.g: mrSearchHelp('scanparams')
%             mrSearchHelp({'defaultInterrogator', 'glm'},'-r -u -n -l')
%            
function [  ]=mrSearchHelp( expression, options )

if nargin < 2
  % '-r -l ;
  options = '-If \.m$ -r -u -n';
end

if nargin < 1
  help mrSearchHelp
  return
end

if ~exist('grep')==2
  disp('(mrSearchHelp) the matlab grep function is not installed')
  disp('see http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=9647')
  return
end

% find out where mrLoadRet is installed on the current machine
loadretlocation = which('mrLoadRet');

% directory above mrLoadRet, to include mrAlign, mrUtilities, etc.
loadretpath = [fileparts(loadretlocation) '/../' ];


% create a line separator
sep=@(x) disp(sprintf(repmat('-',1,75)));

sep();
if iscell(expression)
  fl = grep([ options ' -e'],expression, loadretpath);
elseif isstr(expression)
  fl = grep([ options ' -e ' expression], loadretpath);
else
  disp('(uhoh) expression needs to be a string or cell array of strings!')
end
sep();

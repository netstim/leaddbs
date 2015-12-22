%MRSTRUCT_READ read mrStruct from file.
%		mrStruct = mrstruct_read(fileStr)
%
%		If nargin=0 or the file does not exists, a gui will prompt for user input.
%		'mrStruct': MR-Structure, which will be [] in case of error.
%
%		Examples:
%		mrStruct = mrstruct_read('c:\test\myMrStruct');
%
%		Harald Fischer
%		1/00
%
%		PC




function [mrStruct, fName] = mrstruct_read(fileStr)



%%%%% init and error check
mrStruct    = [];
fName= [];

if nargin==1
   if isempty(findstr(fileStr,'.mat')) & isempty(findstr(fileStr,'.MAT'))
      file1Str = sprintf('%s.mat',fileStr);
      file2Str = sprintf('%s.MAT',fileStr);
      if exist(file1Str)==2
         fileStr = file1Str;
      elseif exist(file2Str)==2
         fileStr = file1Str;
      else
         fileStr = '';
      end
   end
end


if nargin==0 | exist(fileStr)~=2
   [fileStr,dirStr]=uigetfile('*.mat','Load a mrStruct');
   if ischar(fileStr) & ischar(dirStr)
      fileStr = local_full_filename(dirStr,fileStr);
   else
      return; % choosing a file was canceled
   end
end


if isempty(fileStr) | exist(fileStr)~=2
   warning('cannot load mat file');
   return;
end
%%%%% End of: init and error check



%%%%% load
try
   tmpStruct   = load(fileStr);
   fieldsMxStr = fieldnames(tmpStruct);
   fieldStr    = char(fieldsMxStr(1,:));
   mrStruct    = getfield(tmpStruct,fieldStr);
catch
   warning('cannot load mrStruct due to inconcistencies');
end
%%%%% End of: load



%%%%% check for valid struct
isValid = mrstruct_istype(mrStruct);
if ~isValid
   warning('not a valid mrStruct');
   mrStruct = [];
   return;
end
fName= fileStr;
%%%%% check for valid struct




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  local functions                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%% 
%%%%% create full filename
function filePathStr = local_full_filename(pathStr,fileStr)

filePathStr = '';
pathStr     = deblank(pathStr);
fileStr     = deblank(fileStr);

if filesep=='\'
   if (pathStr(length(pathStr))=='\')
      filePathStr = sprintf('%s%s',pathStr,fileStr);
   else
      filePathStr = sprintf('%s\\%s',pathStr,fileStr);
   end
else
   if (pathStr(length(pathStr))=='/')
      filePathStr = sprintf('%s%s',pathStr,fileStr);
   else
      filePathStr = sprintf('%s/%s',pathStr,fileStr);
   end
end
%%%%% End of: create full filename
%%%%% 



% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics

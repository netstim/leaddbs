%MRSTRUCT_WRITE save mrStrut to file.
%		outFileStr = mrstruct_write(mrStruct,fileStr)
%
%		Save mrStruct of file. If nargin==1 or fileStr is not valid
%		a dialog box will prompt the user for input.
%
%		Examples:
%		fileNameStr = mrstruct_write(myExperimentStruct);
%		mrstruct_write(myExperimentStruct,'c:\my_test\myExperiment');
%
%		Harald Fischer
%		1/00
%
%		PC




function outFileStr = mrstruct_write(mrStruct,fileStr)



%%%%% init and error check
outFileStr = '';
if nargin==1,
   [fileStr,dirStr]=uiputfile('*.mat','Save a mrStruct');
   if fileStr~=0
       fileStr = local_full_filename(dirStr,fileStr);   
   end;
elseif nargin==0 || nargin>2,
   warning('wrong number of input arguments');
   return;
end;

isValid = mrstruct_istype(mrStruct);
if ~isValid,
   warning('not a valid mrStruct');
   return;
end;
%%%%% End of: init and error check



%%%%% save
if fileStr~=0,
	save(fileStr,'mrStruct', '-v7.3');
	extensionIndex = findstr(fileStr,'.mat');
	if isempty(extensionIndex),
       fileStr = sprintf('%s.mat',fileStr);
    end;
	if exist(fileStr)==2,
       outFileStr = fileStr;
	else
       warning('saving file failed');
    end;
end;
%%%%% End of: save



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  local functions                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%% 
%%%%% create full filename
function filePathStr = local_full_filename(pathStr,fileStr)

filePathStr = '';
pathStr     = deblank(pathStr);
fileStr     = deblank(fileStr);

if filesep=='\',
   if (pathStr(length(pathStr))=='\'),
      filePathStr = sprintf('%s%s',pathStr,fileStr);
   else
      filePathStr = sprintf('%s\\%s',pathStr,fileStr);
   end;
else
   if (pathStr(length(pathStr))=='/'),
      filePathStr = sprintf('%s%s',pathStr,fileStr);
   else
      filePathStr = sprintf('%s/%s',pathStr,fileStr);
   end;
end;
%%%%% End of: create full filename
%%%%% 



% Copyright (c) May 15th, 2001 by University of Freiburg, Dept. of Radiology, Section of Medical Physics

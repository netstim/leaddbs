function [nFiles, fileList] = countNiftiFiles(dirName,extension)
% [nFiles, fileList] = countNiftiFiles(dirName,extension)
%
% Uses the generic Matlab DIR function to count the Nifti data files
% in the input directory [dirName], and returns a sorted list
%
% Outputs: nFiles      number of Nifti data files in directory
%          fileList    cell array of file names

dS = dir(dirName);
fileList = {};
nList = length(dS);
if nList == 0, return; end

nFiles = 0;
for iList=1:nList
	fName = dS(iList).name;
	[pathstr,bname,ext] = fileparts(fName);
	if (strcmp(lower(ext), extension))
		fileList = [fileList {fName}];
		nFiles = nFiles + 1;
	end
end

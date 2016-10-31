% mlrImageGetAllFilenames.m
%
%        $Id:$ 
%      usage: imageFilenames = mlrImageGetAllFilenames(<dirname>,<'mustNotHaveDotInFilename=0'>)
%         by: justin gardner
%       date: 07/23/09
%    purpose: Function to return all the valid image filenames in a directory. This looks for all nifti
%             files and will only return the .hdr filename for dual files. It will also validate the 
%             headers and not return filenames with invalid nifti headers. It will return either hdr or nii files.
%
%             dirname is the name of the directory to operate on. defaults to current directory
%
function imageFilenames = mlrImageGetAllFilenames(dirname,varargin)

imageFilenames = {};

% check arguments
if ~any(nargin == [0 1 2 3 4])
  help mlrImageGetAllFilenames
  return
end

mustNotHaveDotInFilename = [];
getArgs(varargin,{'mustNotHaveDotInFilename=0'});

if ieNotDefined('dirname'),dirname = pwd;end

% get the direcory
d = dir(dirname);

% list of valid extensions
validExtensions = {'hdr','nii'};

% list of valid extensions that have been gzip'd
validZippedExtensions = {'nii'};

for i = 1:length(d)
  % check for valid extension
  if any(strcmp(getext(d(i).name),validExtensions))
    % check for valid filename
    if mlrImageIsImage(fullfile(dirname,d(i).name))
      imageFilenames{end+1} = d(i).name;
    end
  % check for valid gz extension
  elseif strcmp(getext(d(i).name),'gz')
    uncompressedFilename = stripext(d(i).name);
    % check if the extension (stripped of gz matches validZIppedExtensions
    if any(strcmp(getext(uncompressedFilename),validZippedExtensions))
      % check for valid filenames, also in this check here
      % do not use files in which there exists an uncompressed filename
      thisFilename = fullfile(dirname,d(i).name);
      if ~isfile(stripext(thisFilename)) &&  mlrImageIsImage(thisFilename)
	imageFilenames{end+1} = d(i).name;
      end
    end
  end
end

% remove any filenames that have dots in the middle of them, since this causes weird problems later
% since you can't depend on the dot marking extensions
if mustNotHaveDotInFilename
  imageFilenamesWithoutDot = {};
  for i = 1:length(imageFilenames)
    % check for a ziped extension
    thisFilename = imageFilenames{i};
    if strcmp(getext(thisFilename),'gz')
      thisFilename = stripext(thisFilename);
    end
    % now check for imageFilenames with multiple dots
    if isempty(strfind(stripext(thisFilename),'.'))
      imageFilenamesWithoutDot{end+1} = imageFilenames{i};
    else
      mrWarnDlg(sprintf('(mlrImageGetAllFilenames) Ignoring file %s because it has a . in the filename that does not mark the file extension. If you want to use this file, consider renaming to %s',imageFilenames{i},setext(fixBadChars(stripext(imageFilenames{i}),{'.','_'}),'hdr')));
    end
  end
  imageFilenames = imageFilenamesWithoutDot;
end




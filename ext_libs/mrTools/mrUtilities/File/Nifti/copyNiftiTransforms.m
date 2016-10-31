% copyNiftiTransforms.m
%
%        $Id$
%      usage: copyNiftiTransforms(<copyQform>,<copySform>)
%         by: justin gardner
%       date: 05/04/07
%    purpose: Function to copy qform and sform form one nifti to
%             another. Set copyQform to 1 or 0 to copy qform
%             Likewise for copySform
%
% 
%
function retval = copyNiftiTransforms(copyQform,copySform)

% check arguments
if ~any(nargin == [0 1 2])
  help copyQform
  return
end

if ~exist('copyQform','var'),copyQform = 1;,end
if ~exist('copySform','var'),copySform = 1;,end

% get the filenames
[srcFilename srcPath] = uigetfile('*.hdr','Get source nifti file');
if srcFilename==0,return,end
[destFilename destPath] = uigetfile('*.hdr','Get destination nifti file');
if destFilename==0,return,end

% get the headers
srcHeader = cbiReadNiftiHeader(fullfile(srcPath,srcFilename));
destHeader = cbiReadNiftiHeader(fullfile(destPath,destFilename));

% get fields to copy
copyFields = {};
if copyQform
  copyFields = {'qform44','qform_code','quatern_b','quatern_c','quatern_d','qoffset_x','qoffset_y','qoffset_z'};
end
if copySform
  copyFields = union(copyFields,{'sform44','sform_code','srow_x','srow_y','srow_z'});
end

% show the subjects the qforms and ask if they want to copy
paramsInfo = {};
% set source fields to display
for i = 1:length(copyFields)
  paramsInfo{end+1} = {sprintf('src_%s',copyFields{i}),srcHeader.(copyFields{i}),'editable=0',sprintf('%s from source',copyFields{i})};
end
% set dest fields to display
for i = 1:length(copyFields)
  paramsInfo{end+1} = {sprintf('dest_%s',copyFields{i}),destHeader.(copyFields{i}),'editable=0',sprintf('%s from destination',copyFields{i})};
end
% put up dialog
params = mrParamsDialog(paramsInfo,sprintf('Copy qform from %s to %s?',srcFilename,destFilename));
if isempty(params),return,end

% copy the current dest to a backup
destBackupFilename = sprintf('%s.%s',destFilename,datestr(now,'yymmdd-HHMMSS'));
disp(sprintf('Making backup of %s to %s',destFilename,destBackupFilename));
system(sprintf('cp %s %s',fullfile(destPath,destFilename),fullfile(destPath,destBackupFilename)));

% now copy qform and write out dest header
for i = 1:length(copyFields)
  destHeader.(copyFields{i}) = srcHeader.(copyFields{i});
end
cbiWriteNiftiHeader(destHeader,fullfile(destPath,destFilename));


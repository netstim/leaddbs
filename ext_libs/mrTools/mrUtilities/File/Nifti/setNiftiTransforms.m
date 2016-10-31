% setNiftiTransforms.m
%
%        $Id$
%      usage: setNiftiTransforms(<setQform>,<setSform>)
%         by: justin gardner
%       date: 05/04/07
%    purpose: Function to set qform and sform
%
% 
%
function retval = setNiftiTransforms(setQform,setSform)

% check arguments
if ~any(nargin == [0 1 2])
  help copyQform
  return
end

if ~exist('setQform','var'),setQform = 1;,end
if ~exist('setSform','var'),setSform = 1;,end

% get the filenames
[filename filepath] = uigetfile('*.hdr','Select a nifti file');
if filename==0,return,end

% get the headers
header = cbiReadNiftiHeader(fullfile(filepath,filename));

% show the subjects the qforms and ask if they want to change it
if setQform
  paramsInfo = {};
  paramsInfo{end+1} = {'qform44',header.qform44,'The qform to change'};
  paramsInfo{end+1} = {'qform_code',header.qform_code,'The qform code'};
  paramsInfo{end+1} = {'transform',eye(4),'A transformation matrix to multiply the above matrix by'};
  params = mrParamsDialog(paramsInfo,sprintf('Change qform?'));
  if isempty(params),return,end

  % set the fields
  header = cbiSetNiftiQform(header,params.qform44*params.transform);
  header.qform_code = params.qform_code;
end

% show the subjects the qforms and ask if they want to change it
if setSform
  paramsInfo = {};
  paramsInfo{end+1} = {'sform44',header.sform44,'The sform to change'};
  paramsInfo{end+1} = {'sform_code',header.sform_code,'The sform code'};
  paramsInfo{end+1} = {'transform',eye(4),'A transformation matrix to multiply the above matrix by'};
  params = mrParamsDialog(paramsInfo,sprintf('Change sform?'));
  if isempty(params),return,end

  % set the fields
  header = cbiSetNiftiSform(header,params.sform44*params.transform);
  header.sform_code = params.sform_code;
end

% copy the current dest to a backup
backupFilename = sprintf('%s.%s',filename,datestr(now,'yymmdd-HHMMSS'));
disp(sprintf('Making backup of %s to %s',filename,backupFilename));
system(sprintf('cp %s %s',fullfile(filepath,filename),fullfile(filepath,backupFilename)));

% and write it out
cbiWriteNiftiHeader(header,fullfile(filepath,filename));

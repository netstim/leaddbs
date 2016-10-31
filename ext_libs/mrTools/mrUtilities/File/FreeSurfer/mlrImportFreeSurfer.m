% mlrImportFreeSurfer.m
%
%        $Id$
%      usage: mlrImportFreeSurfer(varargin)
%         by: eli merriam
%       date: 07/11/07
%    purpose: 
%
function retval = mlrImportFreeSurfer(varargin)

% check arguments
% if ~any(nargin == [0])
%   help mlrImportFreeSurfer
%   return
% end

if ~exist('mrParamsDialog')
  disp(sprintf('(mlrImportFreeSurfer) You must have mrTools in your path to run this'));
  return
end

mriConvert = 'mri_convert';
[retval retstr] = system('which mri_convert');
if retval == 1
  mrWarnDlg('(mlrImportFreeSurfer) Could not find FreeSurfer command mri_convert which is needed to convert the FreeSurfer anatomy file to a nifti file. This is usually in the bin directory under your freesurfer installation. You may need to install freesurfer and add that to your path. See instructions on wiki http://gru.stanford.edu/doku.php/mrTools/howTo#installation','Yes');
  mriConvert = [];
end

% evaluate the arguments
eval(evalargs(varargin));

% check directories
if ~isdir('surf'),disp(sprintf('(mlrImportFreeSurfer) Could not find surf directory'));return,end
if ~isdir('mri'),disp(sprintf('(mlrImportFreeSurfer) Could not find mri directory'));return,end

if ieNotDefined('baseName'), baseName = getLastDir(pwd);end

wmFile    = 'smoothwm';
gmFile    = 'pial';
infFile   = 'inflated';
curvFile  = 'curv';
anatFile  = 'T1.mgz';
hemi      = {'lh', 'rh'};
hemiNames = {'left', 'right'};
outDir    = fullfile(pwd,'surfRelax');
volumeCropSize = [176 256 256];

paramsInfo = {...
    {'freeSurferDir',pwd,'directory where the freeSurfer files live'}, ...
    {'outDir', outDir,'directory that OFF surfaces will be written to'}, ...
    {'wmFile', wmFile, 'name of the surface defining the white/gray boundry'}, ...
    {'gmFile', gmFile, 'name of the surface defining the gray/pial boundry'}, ...
    {'curvFile', curvFile, 'name of the file specifing the curvature'}, ...
    {'infFile', infFile, 'name of the inflated surface'}, ...
    {'anatFile', anatFile, 'name of the base anatomy file'}, ...
    {'baseName', baseName, 'subject initials'}, ...
    {'volumeCropSize',volumeCropSize, 'Size to crop the volume anatomy to'},...
             };

% get the parameters from the user
params = mrParamsDialog(paramsInfo,'mlrImportFreeSurfer',1.5);

% return if user hit cancel
if isempty(params)
  return
end

if ~isdir(outDir)
  mkdir(outDir);
end

% import the white and gray matter surface, as well as the inflated surface and curvature
disp(sprintf('(mlrImportFreeSurfer) Converting FreeSurfer surfaces to OFF format'))
for i = 1:length(hemi)
  % convert inner surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.wmFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_WM.off'));
  if isfile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find inner (white matter) surface %s',getLastDir(surfFile,2)));
  end

  % convert outer surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.gmFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_GM.off'));
  if isfile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find outer (pial/gray matter) surface %s',getLastDir(surfFile,2)));
  end

  % convert inflated surface
  surfFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.infFile));
  outFile = fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Inf.off'));
  if isfile(surfFile)
    freeSurfer2off(surfFile, outFile, params.volumeCropSize);
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find inflated surface %s',getLastDir(surfFile,2)));
  end

  curvFile = fullfile(params.freeSurferDir, 'surf', strcat(hemi{i}, '.', params.curvFile));
  if isfile(curvFile)
    [curv, fnum] = freesurfer_read_curv(curvFile);
    saveVFF(fullfile(params.outDir, strcat(params.baseName, '_', hemiNames{i}, '_Curv.vff')), -curv)
  else
    disp(sprintf('(mlrImportFreeSurfer) Could not find curvature file %s',getLastDir(curvFile)));
  end
end

% convert the volume to mlr volume
anatFile = fullfile(params.freeSurferDir, 'mri', params.anatFile);
niftiExt = mrGetPref('niftiFileExtension');
switch niftiExt
  case '.nii'
    out_type='nii';
  case '.img'
    out_type='nifti1';
end
outFile = fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt));
if isfile(anatFile)
  if ~isempty(mriConvert)
    disp(sprintf('(mlrImportFreeSurfer) Converting the volume anatomy to nifti format'))
    setenv('LD_LIBRARY_PATH', '/usr/pubsw/packages/tiffjpegglut/current/lib:/opt/local/lib:/usr/local/lib:/opt/local/lib')
    system(sprintf('%s --out_type %s --out_orientation RAS --cropsize %i %i %i %s %s',mriConvert,out_type,params.volumeCropSize(1),params.volumeCropSize(2),params.volumeCropSize(3),anatFile,outFile));
  else
    disp(sprintf('(mlrImportFreeSurfer) !!!! No mlr_convert command available, so not making canonical anatomy !!!!'));
  end
  if ~isfile(outFile)
    mrWarnDlg(sprintf('(mlrImportFreeSurfer) !!!! Canonical anatomy not created !!!!'));
  end
    
%   h = mlrImageReadNiftiHeader(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt)));
%   h.qform44(1,4) = -(h.dim(2)-1)/2;
%   h = cbiSetNiftiQform(h, h.qform44);
%   h = cbiSetNiftiSform(h, h.qform44);
%   cbiWriteNiftiHeader(h, fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp', niftiExt)));
else
  disp(sprintf('(mlrImportFreeSurfer) Could not find volume %s',getLastDir(anatFile,2)));
end


return;

disp(sprintf('(mlrImportFreeSurfer) Constructing the gray graph'))
setenv('DYLD_LIBRARY_PATH', '/Users/eli/src/TFI/sw/lib/')
for i = 1:length(hemi)
    system(sprintf('surf2graph -verbose -gimage %s -gmin 40 -gmax 120 -expand -1 -gsurface %s -gray %s -savelabels %s %s', ...
                   strcat(params.baseName, '_', 'mprage_pp.hdr'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_GM.off'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '.Gray'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_graylabels.img'), ...
                   strcat(params.baseName, '_', hemiNames{i}, '_WM.off')...
                   ));
end


return;
% disp(sprintf('(mlrImportFreeSurfer) Fixing nifti header so that the volume is centered'))
% h = mlrImageReadNiftiHeader(sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));
% h.qform44(1,4) = (h.dim(2)-1)/2;
% h.qoffset_x = -1*(h.dim(2)-1)/2;
% h.sform_code = 0;
% cbiWriteNiftiHeader(h, sprintf(fullfile(params.outDir, strcat(params.baseName, '_', 'mprage_pp.hdr'))));


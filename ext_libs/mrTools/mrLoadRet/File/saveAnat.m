function pathStr = saveAnat(view,anatomyName,confirm,saveAs,savePath)
%
% saveAnat(view,[anatomyName],[confirm],[saveAs],[savePath])
%
% Saves a base anatomy as a nifti file: anatomyName.img
%
% anatomyName: Can be either the name or the number (in which case it is
%          converted to the name). Default: current base volume.
% confirm: If filename already exists, prompt user to over-write. 
%          Default: uses 'overwritePolicy' preference.
%  saveAs: If 1, then asks user for where to put anatomy (default=0)
% savePath: If set then saves the anatomy to the specified path
% 
% 
%
% djh, 7/2006
% $Id$	

pathStr = [];
if ieNotDefined('anatomyName')
  anatomyNum = viewGet(view,'currentBase');
  anatomyName = viewGet(view,'baseName',anatomyNum);
end
if isstr(anatomyName)
  anatomyNum = viewGet(view,'baseNum',anatomyName);
elseif isnumeric(anatomyName)
  anatomyNum = anatomyName;
  anatomyName = viewGet(view,'baseName',anatomyNum);
else
  myErrorDlg(['Invalid base volume (anatomy) name: ',anatomyName]);
end

if ieNotDefined('confirm')
  pref = mrGetPref('overwritePolicy');
  confirm = ~strcmp(pref,'Overwrite');
end
if ieNotDefined('saveAs')
  saveAs = 0;
end
if ieNotDefined('savePath')
  savePath = [];
end

% Extract data and hdr
baseVolume = viewGet(view,'baseVolume',anatomyNum);
data = baseVolume.data;
hdr = baseVolume.hdr;

% also get the base structure, and remove the data and hdr
% this will get saved along with the nifti file
base = viewGet(view,'base',anatomyNum);
base.data = [];
base.hdr = [];

% set the file extension
niftiFileExtension = mrGetPref('niftiFileExtension');
if isempty(niftiFileExtension)
  niftiFileExtension = '.img';
end

% Path
if saveAs
  % get default path
  global saveAnatPath;
  if isempty(saveAnatPath)
    saveAnatPath = mrGetPref('volumeDirectory');
  end
  % ask user for where to put it
  [anatomyName,saveAnatPath] = uiputfile({sprintf('*%s',niftiFileExtension),'Nifti files'},sprintf('Save base anatomy %s as...',base.name),fullfile(saveAnatPath,[anatomyName niftiFileExtension]));
  % if cancel
  if anatomyName == 0
    saveAnatPath = mrGetPref('volumeDirectory');
    return
  end
  anatomyName = stripext(anatomyName);
  % make into a path
  pathStr = fullfile(saveAnatPath,[anatomyName,niftiFileExtension]);
elseif ~isempty(savePath)
  pathStr = fullfile(savePath,[anatomyName,niftiFileExtension]);
else
  pathStr = fullfile(viewGet(view,'anatomydir'),[anatomyName,niftiFileExtension]);
end  
% Write, though check for over-writing
saveFlag = 'Yes';
if exist([pathStr,'.mat'],'file')
  if confirm
    saveFlag = questdlg([pathStr,' already exists. Overwrite?'],...
			'Save Overlay?','Yes','No','No');
  end
end
if strcmp(saveFlag,'Yes')
  fprintf('(saveAnat) Saving %s...',pathStr);
  [byteswritten,hdr] = cbiWriteNifti(pathStr,data,hdr,'float32');
  % also write out the base structure as a .mat file
  eval(sprintf('save %s.mat base',stripext(pathStr)));
  fprintf('done\n');
else
  fprintf('Anatomy not saved...');
end
return;

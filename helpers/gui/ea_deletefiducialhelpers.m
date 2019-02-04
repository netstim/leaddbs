function ea_deletefiducialhelpers(~,~,handles)
answ=questdlg('Are you sure you want to delete fiducial marker(s) for selected patient(s)?','Delete Fiducial Markers','Yes','Abort','Abort');
if strcmp(answ,'Abort')
    return
end
answ=questdlg('Do you also wish to delete corresponding fiducials in template space?','Delete Fiducial Markers','Yes','No','No');
switch lower(answ)
    case 'yes'
   deltemptoo=1;     
    case 'no'
           deltemptoo=0;
end

uipatdir=getappdata(handles.leadfigure,'uipatdir');
cnt=1;
for pt=1:length(uipatdir)
    fids=dir([uipatdir{pt},filesep,'fiducials',filesep,'*.nii.gz']);
    
    for fi=1:length(fids)
        % pt folder
        delete([uipatdir{pt},filesep,'fiducials',filesep,fids(fi).name]);
        
        % templates
        if exist([ea_space,'fiducials',filesep,fids(fi).name],'file') && deltemptoo
            delete([ea_space,'fiducials',filesep,fids(fi).name]);
        end
    end
    
end


function smoothgzip(pathn,filen)
matlabbatch{1}.spm.spatial.smooth.data = {fullfile(pathn,filen)};
matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
matlabbatch{1}.spm.spatial.smooth.dtype = 8;
matlabbatch{1}.spm.spatial.smooth.im = 1;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',{matlabbatch});
movefile(fullfile(pathn,['s',filen]),fullfile(pathn,filen));
gzip(fullfile(pathn,filen));
delete(fullfile(pathn,filen));



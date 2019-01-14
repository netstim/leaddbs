function ea_addfiducialhelper(~,~,handles)

uipatdir=getappdata(handles.leadfigure,'uipatdir');

% choose in template:

spm('defaults', 'fmri');
Fgraph = spm_figure('GetWin', 'Graphics');
Finter = spm('FnUIsetup','Select Fiducials', 0);

figure(Fgraph); clf;

spm_orthviews('Reset');
spm_orthviews('Image', [ea_space,'t1.nii']);
colormap('gray');
cameratoolbar('resetcamera')
cameratoolbar('close')
rotate3d off;

if spm_input(['Select fiducial point on template and click'] , 1,'OK|Retry', [1,0], 1)
    fid.template = spm_orthviews('Pos')';
end
uuid=ea_generate_uuid;


% choose in patient(s):

for pt=1:length(uipatdir)
    figure(Fgraph); clf;
    options=ea_getptopts([uipatdir{pt},filesep]);
    [~,presentfiles]=ea_assignpretra(options);
    spm_orthviews('Reset');
    ptspace{pt}=[uipatdir{pt},filesep,presentfiles{1}];
    spm_orthviews('Image', ptspace{pt});
    colormap('gray');
    cameratoolbar('resetcamera')
    cameratoolbar('close')
    rotate3d off;
    
    if spm_input(['Select corresponding fiducial point on patient anchor modality and click'] , 1,'OK|Retry', [1,0], 1)
        fid.patient(pt,:) = spm_orthviews('Pos')';
    end
    
end
close(Fgraph);
close(Finter);
% define in template:
ea_mkdir([ea_space,'ficudials']);
ea_spherical_roi([ea_space,'ficudials',filesep,uuid,'.nii'],fid.template,10,0,[ea_space,'t1.nii']);
smoothgzip([ea_space,'ficudials'],[uuid,'.nii']);

% define in pt:
for pt=1:length(uipatdir)
    ea_mkdir([uipatdir{pt},filesep,'ficudials']);
    ea_spherical_roi([uipatdir{pt},filesep,'ficudials',filesep,uuid,'.nii'],fid.patient(pt,:),10,0,ptspace{pt});
    smoothgzip([uipatdir{pt},filesep,'ficudials'],[uuid,'.nii']);
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



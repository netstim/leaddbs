function ea_gencoregcheckfigs(directory,scrf,options)
disp(['Exporting coregistration check images to ',directory,'checkreg/scrf...']);
[~,filespresent]=ea_assignpretra(options);

primanat=[directory,'scrf',filesep,filespresent{1}];
setenv('FSLOUTPUTTYPE','NIFTI');


fis2anat={[directory,'scrf',filesep,scrf,'movim.nii']};

basedir=[ea_getearoot,'ext_libs',filesep,'fsl',filesep];
if ispc
    SLICER = [basedir, 'slicer.exe'];
else
    SLICER = [basedir, 'slicer.', computer('arch')];    
end
cnt=1;
try
    nm=ea_cleanmethodname(['_',ea_whichnormmethod(directory)]);
catch
    nm='';
end
try
    cm=load([directory,'ea_coregctmethod_applied.mat']);
    cm=ea_cleanmethodname(['_',cm.coregct_method_applied{end}]);
catch
    cm='';
end
try
    mm=load([directory,'ea_coregmrmethod_applied.mat']);
    mm=ea_cleanmethodname(['_',mm.coregmr_method_applied{end}]);
catch
    mm='';
end
for fi=1:length(fis2anat)
    [~,fname]=fileparts(fis2anat{fi});
    [~,rfname]=fileparts(primanat);
    switch [fname,'.nii'] % cannot use options.modality here since also in CT imaging, e.g. anat_t1 or anat_pd could be used using coregmrmethod applied.
        case options.prefs.rawctnii_unnormalized
            suffx=cm; % CT suffix
        otherwise
            suffx=mm; % MR suffix
    end
    if isempty(scrf)
    ofname{cnt}='standard.png';
    else
       ofname{cnt}='scrf.png';
    end
    if exist(fis2anat{fi},'file')
        cmd{cnt}=[SLICER,' ',ea_path_helper(fis2anat{fi}),' ',ea_path_helper(primanat),' -a ',ea_path_helper([directory,'scrf',filesep,ofname{cnt}])];
        cnt=cnt+1;
    end
end

if ~exist([directory,'scrf',],'file')
    mkdir([directory,'scrf',]);
end

setenv('FSLOUTPUTTYPE','NIFTI')
if exist('cmd','var')
for c=1:length(cmd)
    if ~ispc
        system(['bash -c "', cmd{c}, '"']);
    else
        system(cmd{c});
    end
end
disp('Done.');
end



function name=ea_cleanmethodname(name)
name(strfind(name,':'))=[];
name(strfind(name,' '))='_';


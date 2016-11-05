function ea_gencoregcheckfigs(options)
directory=[options.root,options.patientname,filesep];
[~,filespresent]=ea_assignpretra(options);

primanat=[directory,filespresent{1}];
mnihires=[ea_getearoot,'templates',filesep,'mni_hires.nii'];
try
    oanat=filespresent(2:end);
catch
	oanat={}; 
end
for oa=1:length(oanat)
	oanat{oa}=[directory,oanat{oa}]; 
end
switch options.modality
    case 1
        fis2anat={[directory,options.prefs.tranii_unnormalized],...
            [directory,options.prefs.cornii_unnormalized],...
            [directory,options.prefs.sagnii_unnormalized]};
        fis2mni={[directory,options.prefs.gprenii],...
            [directory,options.prefs.gtranii],...
            [directory,options.prefs.gcornii],...
            [directory,options.prefs.gsagnii]};
    case 2
        fis2anat={[directory,options.prefs.ctnii_coregistered]};
        fis2mni={[directory,options.prefs.gprenii],...
            [directory,options.prefs.gctnii]};
end
fis2anat=[fis2anat,oanat];

basedir=[ea_getearoot,'ext_libs',filesep,'fsl',filesep];
if ispc
    SLICER = [basedir, 'slicer.exe'];
else
    SLICER = [basedir, 'slicer.', computer('arch')];    
end
cnt=1;
try
    nm=load([directory,'ea_normmethod_applied.mat']);
    nm=nm.norm_method_applied{end};
catch
    nm='';
end
try
    cm=load([directory,'ea_coregctmethod_applied.mat']);
    cm=cm.coregct_method_applied{end};
catch
    cm='';
end
try
    mm=load([directory,'ea_coregmrmethod_applied.mat']);
    mm=mm.coregmr_method_applied{end};
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
    ofname{cnt}=[fname,'2',rfname,'_',cm,'.png'];
    
    cmd{cnt}=[SLICER,' ',ea_path_helper(fis2anat{fi}),' ',ea_path_helper(primanat),' -a ',ea_path_helper([directory,'checkreg',filesep,ofname{cnt}])];    

    cnt=cnt+1;
end
for fi=1:length(fis2mni)
    [~,fname]=fileparts(fis2mni{fi});
    [~,rfname]=fileparts(mnihires);
    ofname{cnt}=[fname,'2',rfname,'_',nm,'.png'];
    cmd{cnt}=[SLICER,' ',ea_path_helper(fis2mni{fi}),' ',ea_path_helper(mnihires),' -a ',ea_path_helper([directory,'checkreg',filesep,ofname{cnt}])];

    cnt=cnt+1;
end
if ~exist([directory,'checkreg'],'file')
    mkdir([directory,'checkreg']);
end

setenv('FSLOUTPUTTYPE','NIFTI')
for c=1:length(cmd)
    if ~ispc
        system(['bash -c "', cmd{c}, '"']);
    else
        system(cmd{c});
    end
end

function ea_gencoregcheckfigs(options)
fprintf('\nExporting coregistration check images to %scheckreg...\n', [options.root,options.patientname,filesep]);
directory=[options.root,options.patientname,filesep];
[~,filespresent]=ea_assignpretra(options);

primanat=[directory,filespresent{1}];
mnihires=[ea_space(options),options.primarytemplate,'.nii'];
setenv('FSLOUTPUTTYPE','NIFTI');
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
        if exist([directory,options.prefs.ctnii_coregistered],'file')
            try ea_tonemapct_file(options,'native'); end
        end

        if exist([directory,options.prefs.gctnii],'file')
            try ea_tonemapct_file(options,'mni'); end
        end

        fis2anat={[directory,'tp_',options.prefs.ctnii_coregistered]};
        fis2mni={[directory,options.prefs.gprenii],...
            [directory,'tp_',options.prefs.gctnii]};
end
fis2anat=[fis2anat,oanat'];


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
if ~exist([directory,'checkreg'],'file')
    mkdir([directory,'checkreg']);
end
for fi=1:length(fis2anat)
    [~,fname]=fileparts(fis2anat{fi});
    [~,rfname]=fileparts(primanat);
    switch [fname,'.nii'] % cannot use options.modality here since also in CT imaging, e.g. anat_t1 or anat_pd could be used using coregmrmethod applied.
        case options.prefs.tp_ctnii_coregistered
            suffx=cm; % CT suffix
        otherwise
            suffx=mm; % MR suffix
    end
    ofname{cnt}=[fname,'2',rfname,suffx,'.png'];

    if exist(fis2anat{fi},'file')
        
        ea_gencheckregpair(stripext(fis2anat{fi}),stripext(primanat),[directory,'checkreg',filesep,ofname{cnt}]);
    end
end
for fi=1:length(fis2mni)
    [~,fname]=fileparts(fis2mni{fi});
    [~,rfname]=fileparts(mnihires);
    ofname{cnt}=[fname,'2',rfname,nm,'.png'];
    if exist(fis2mni{fi},'file')
                ea_gencheckregpair(stripext(fis2mni{fi}),stripext(mnihires),[directory,'checkreg',filesep,ofname{cnt}]);
    end
end

disp('Done.');

function fn=stripext(fn)
[pth,fn,ext]=fileparts(fn);
fn=fullfile(pth,fn);

function name=ea_cleanmethodname(name)
name(strfind(name,':'))=[];
name(strfind(name,' '))='_';


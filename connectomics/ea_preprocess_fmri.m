function ea_preprocess_fmri(options)


directory=[options.root,options.patientname,filesep];

V=spm_vol([directory,options.prefs.rest]);
signallength=length(V);

%% run sequence of proxyfunctions (below):
ea_realign_fmri(signallength,options); % realign fMRI

ea_newseg(directory,options.prefs.prenii_unnormalized,0,options,1); % Segment anat

ea_coreg_pre2fmri(options); % register pre 2 fmri (for timecourse-extraction).

ea_smooth_fmri(signallength,options); % slightly smooth fMRI data


function ea_realign_fmri(signallength,options)
%% realign fmri.
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'r',options.prefs.rest],'file')

    tmpdir=ea_getleadtempdir;
    uuid=ea_generate_uuid;
    copyfile([directory,options.prefs.rest],[tmpdir,uuid,'.nii']);


    disp('Realignment of rs-fMRI data...');
    filetimepts=cell(signallength,1);
    for i = 1:signallength
        filetimepts{i}=[tmpdir,uuid,'.nii',',',num2str(i)];
    end

    matlabbatch{1}.spm.spatial.realign.estwrite.data = {filetimepts};
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight = {''};
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
    spm_jobman('run',{matlabbatch});
    clear matlabbatch;
    disp('Done.');

    movefile([tmpdir,'r',uuid,'.nii'],[directory,'r',options.prefs.rest]);
    movefile([tmpdir,'mean',uuid,'.nii'],[directory,'mean',options.prefs.rest]);
    movefile([tmpdir,'rp_',uuid,'.txt'],[directory,'rp_',ea_stripex(options.prefs.rest),'.txt']);

    ea_delete([tmpdir,uuid,'.nii']);
end


function ea_coreg_pre2fmri(options)
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'r',ea_stripex(options.prefs.rest),'_',options.prefs.prenii_unnormalized],'file')
    ea_coreg2images_generic(options,...
        [directory,options.prefs.prenii_unnormalized],...
        [directory,'r',options.prefs.rest],...
        [directory,'r',ea_stripex(options.prefs.rest),'_',options.prefs.prenii_unnormalized],...
        {[directory,'c1',options.prefs.prenii_unnormalized];...
        [directory,'c2',options.prefs.prenii_unnormalized];...
        [directory,'c3',options.prefs.prenii_unnormalized]},1,[],1);
    for t=1:3
        movefile([directory,'rc',num2str(t),options.prefs.prenii_unnormalized],[directory,'r',ea_stripex(options.prefs.rest),'_c',num2str(t),options.prefs.prenii_unnormalized]);
    end
end


function ea_smooth_fmri(signallength,options)
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'sr',options.prefs.rest],'file')
    filetimepts=cell(signallength,1);
    for i = 1:signallength
        filetimepts{i}=[directory,'r',options.prefs.rest,',',num2str(i)];
    end

    matlabbatch{1}.spm.spatial.smooth.data = filetimepts;
    matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear jobs matlabbatch
end

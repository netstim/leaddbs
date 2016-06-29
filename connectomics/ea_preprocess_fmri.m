function ea_preprocess_fmri(options)


directory=[options.root,options.patientname,filesep];

V=spm_vol([directory,options.prefs.rest]);
signallength=length(V);

%% run sequence of proxyfunctions (below):
ea_realign_fmri(signallength,options); % realign fMRI

% new segment (fMRI, preoperative anatomical image):
fis={['r',options.prefs.rest,',1'],options.prefs.prenii_unnormalized};
for fi=1:length(fis)
    ea_newsegment_proxy(fis{fi},options); % new Segment fMRI
end

ea_coreg_pre2fmri(options); % register pre 2 fmri (for timecourse-extraction).
ea_smooth_fmri(signallength,options); % slightly smooth fMRI data



function ea_realign_fmri(signallength,options)
%% realign fmri.
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'r',options.prefs.rest],'file');

    disp('Realignment of rs-fMRI data...');
    filetimepts=cell(signallength,1);
    for i = 1:signallength
        filetimepts{i}=[directory,options.prefs.rest,',',num2str(i)];
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
    delete([directory,'mean',options.prefs.rest]);
    disp('Done.');
end


function ea_newsegment_proxy(file,options)
directory=[options.root,options.patientname,filesep];

%% new segment.
if ~exist([directory,'c1',file],'file');
    ea_newseg(directory,file,0,options);
    delete([directory,'c4',file]);
    delete([directory,'c5',file]);
    [~,rf]=fileparts(file);
    delete([directory,rf,'_seg8.mat']);

end

function ea_coreg_pre2fmri(options)
directory=[options.root,options.patientname,filesep];
[~,rf]=fileparts(options.prefs.rest);
if ~exist([directory,'rr',rf,options.prefs.prenii_unnormalized],'file')
    %% coreg mprage to fMRI (for GM-map)

    copyfile([directory,options.prefs.prenii_unnormalized],[directory,'k',options.prefs.prenii_unnormalized])
    copyfile([directory,'c1',options.prefs.prenii_unnormalized],[directory,'kc1',options.prefs.prenii_unnormalized]) %% use copies for coregistration to leave original files untouched.
    copyfile([directory,'c2',options.prefs.prenii_unnormalized],[directory,'kc2',options.prefs.prenii_unnormalized]) %% use copies for coregistration to leave original files untouched.
    copyfile([directory,'c3',options.prefs.prenii_unnormalized],[directory,'kc3',options.prefs.prenii_unnormalized]) %% use copies for coregistration to leave original files untouched.

    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[directory,'r',options.prefs.rest,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[directory,'k',options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {[directory,'kc1',options.prefs.prenii_unnormalized,',1'];
        [directory,'kc2',options.prefs.prenii_unnormalized,',1'];
        [directory,'kc3',options.prefs.prenii_unnormalized,',1']
        };
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = ['rr',rf];


    try
        jobs{1}=matlabbatch;
        spm_jobman('run',jobs);
    end

    clear matlabbatch jobs;
    % cleanup:

    movefile([directory,'rr',rf,'k',options.prefs.prenii_unnormalized],[directory,'rr',rf,options.prefs.prenii_unnormalized],'f')
    movefile([directory,'rr',rf,'kc1',options.prefs.prenii_unnormalized],[directory,'rr',rf,'c1',options.prefs.prenii_unnormalized],'f') %% restore original files..
    movefile([directory,'rr',rf,'kc2',options.prefs.prenii_unnormalized],[directory,'rr',rf,'c2',options.prefs.prenii_unnormalized],'f') %% restore original files..
    movefile([directory,'rr',rf,'kc3',options.prefs.prenii_unnormalized],[directory,'rr',rf,'c3',options.prefs.prenii_unnormalized],'f') %% restore original files..

    movefile([directory,'k',options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized])
    movefile([directory,'kc1',options.prefs.prenii_unnormalized],[directory,'c1',options.prefs.prenii_unnormalized]) %% restore original files..
    movefile([directory,'kc2',options.prefs.prenii_unnormalized],[directory,'c2',options.prefs.prenii_unnormalized]) %% restore original files..
    movefile([directory,'kc3',options.prefs.prenii_unnormalized],[directory,'c3',options.prefs.prenii_unnormalized]) %% restore original files..
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

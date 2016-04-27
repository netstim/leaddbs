function varargout=ea_normalize_ants_multimodal(options)
% This is a function that normalizes both a copy of transversal and coronar
% images into MNI-space. The goal was to make the procedure both robust and
% automatic, but still, it must be said that normalization results should
% be taken with much care because all reconstruction results heavily depend
% on these results. Normalization of DBS-MR-images is especially
% problematic since usually, the field of view doesn't cover the whole
% brain (to reduce SAR-levels during acquisition) and since electrode
% artifacts can impair the normalization process. Therefore, normalization
% might be best archieved with other tools that have specialized on
% normalization of such image data.
%
% The procedure used here uses the ANTs Syn approach to map a patient's
% brain to MNI space directly.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='Advanced Normalization Tools (ANTs) SyN - multimodal (T1,T2, PD & FA)';
    varargout{2}={'SPM8','SPM12'};
    return
end

uset1=1; % set to zero if you do not wish to use T1 data for normalization even if present.
usepd=0; % set to zero if you do not wish to use PD data for normalization even if present.
usefa=1; % set to zero if you do not wish to use FA data for normalization even if present.

directory=[options.root,options.patientname,filesep];

cnt=1;

% T1
if uset1
    if exist([directory,options.prefs.prenii_unnormalized_t1],'file')
                disp('Including T1 data for (grey-matter) normalization');
        ea_coreg2images(options,[directory,options.prefs.prenii_unnormalized_t1],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized_t1]);
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_t1.nii'];
        from{cnt}=[directory,options.prefs.prenii_unnormalized_t1];
        weights(cnt)=1.25;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

% PD
if usepd
    if exist([directory,options.prefs.prenii_unnormalized_pd],'file')
        disp('Including PD data for (grey-matter) normalization');
        ea_coreg2images(options,[directory,options.prefs.prenii_unnormalized_pd],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized_pd]);
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_pd.nii'];
        from{cnt}=[directory,options.prefs.prenii_unnormalized_pd];
        weights(cnt)=1.25;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

if usefa
    % check for presence of FA map
    if ~exist([directory,options.prefs.fa2anat],'file')
        if ~exist([directory,options.prefs.fa],'file')
            if ~exist([directory,options.prefs.dti],'file')
                disp('No dMRI data has been found. Proceeding without FA');
            else
                ea_isolate_fa(directory,options);
            end
        end
        
        if exist([directory,options.prefs.fa],'file') % recheck if has been built.
            ea_coreg2images(options,[directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.fa2anat]);
        end
    end
    if exist([directory,options.prefs.fa2anat],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_fa.nii'];
        from{cnt}=[directory,options.prefs.fa2anat];
        weights(cnt)=1;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

masks=segmentall(from,options); % masks not yet implemented.

% The convergence criterion for the multivariate scenario is a slave to the last metric you pass on the ANTs command line.
to{cnt}=[options.earoot,'templates',filesep,'mni_hires.nii'];
from{cnt}=[directory,options.prefs.prenii_unnormalized];
weights(cnt)=15;
metrics{cnt}='MI';
cnt=cnt+1;

% Do the coreg part for postoperative images:

ea_coregmr(options,options.prefs.normalize.coreg);

weights=weights/sum(weights); % normalize to sum of 1

ea_ants_nonlinear(to,...
    from,...
    [directory,options.prefs.prenii],weights,metrics,options);

ea_apply_normalization(options);

function masks=segmentall(from,options)
directory=[fileparts(from{1}),filesep];
for fr=1:length(from)
    [~,fn,ext]=fileparts(from{fr});
    switch [fn,ext]
        case options.prefs.fa2anat
            if ~exist([directory,'tc2',options.prefs.prenii_unnormalized],'file')
                if ~exist([directory,'c2',options.prefs.prenii_unnormalized],'file')
                    ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
                end
                % assume that tc2 doesn't exist
                nii=ea_load_nii([directory,'c2',options.prefs.prenii_unnormalized]);
                nii.img=nii.img>0.7;
                nii.fname=[directory,'tc2',options.prefs.prenii_unnormalized];
                spm_write_vol(nii,nii.img);
            end
            masks{fr,1}=[options.earoot,'templates',filesep,'mni_hires_c2mask.nii'];
            masks{fr,2}=[directory,'tc2',options.prefs.prenii_unnormalized];
        otherwise
            if ~exist([directory,'tc1',options.prefs.prenii_unnormalized],'file')
                if ~exist([directory,'c1',options.prefs.prenii_unnormalized],'file')
                    ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
                end
                % assume that tc1 doesn't exist
                nii=ea_load_nii([directory,'c1',options.prefs.prenii_unnormalized]);
                nii.img=nii.img>0.3;
                nii.fname=[directory,'tc1',options.prefs.prenii_unnormalized];
                spm_write_vol(nii,nii.img);
            end
            masks{fr,1}=[options.earoot,'templates',filesep,'mni_hires_c1mask.nii'];
            masks{fr,2}=[directory,'tc1',options.prefs.prenii_unnormalized];
    end
end

if ~exist([directory,'tc1c2',options.prefs.prenii_unnormalized],'file')
    Vc1=ea_load_nii([directory,'tc1',options.prefs.prenii_unnormalized]);
    Vc2=ea_load_nii([directory,'tc2',options.prefs.prenii_unnormalized]);
    Vc1.img=Vc1.img+Vc2.img;
    Vc1.fname=[directory,'tc1c2',options.prefs.prenii_unnormalized];
    spm_write_vol(Vc1,Vc1.img);
end

try delete([directory,'c1',options.prefs.prenii_unnormalized]); end
try delete([directory,'c2',options.prefs.prenii_unnormalized]); end
try delete([directory,'c3',options.prefs.prenii_unnormalized]); end
try delete([directory,'c4',options.prefs.prenii_unnormalized]); end
try delete([directory,'c5',options.prefs.prenii_unnormalized]); end
try delete([directory,'c6',options.prefs.prenii_unnormalized]); end
try delete([directory,'y_',options.prefs.prenii_unnormalized]); end
try delete([directory,'iy_',options.prefs.prenii_unnormalized]); end
[~,fna]=fileparts(options.prefs.prenii_unnormalized);
try delete([directory,fna,'_seg8.mat']); end

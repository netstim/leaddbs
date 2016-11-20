function varargout=ea_coreg_all_mri(options,usebrainmask)

% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn


if ischar(options) % return name of method.
    varargout{1}='Advanced Normalization Tools (ANTs) SyN - multimodal (T1,T2, PD & FA)';
    varargout{2}={'SPM8','SPM12'};
    return
end

if ~exist('coregonly','var')
    coregonly=0;
end
if ~exist('includeatlas','var')
    includeatlas=0;
end

uset1=1; % set to zero if you do not wish to use T1 data for normalization even if present.
usepd=1; % set to zero if you do not wish to use PD data for normalization even if present.
usefa=1; % set to zero if you do not wish to use FA data for normalization even if present.


directory=[options.root,options.patientname,filesep];
cnt=1;

if usebrainmask
    bprfx='b';
else
    bprfx='';
end


% T1
if uset1 && ~strcmp(options.primarytemplate,'_t1')
    if exist([directory,options.prefs.prenii_unnormalized_t1],'file')
        disp('Including T1 data for (grey-matter) normalization');
        if ~includeatlas

        ea_coreg2images(options,[directory,options.prefs.prenii_unnormalized_t1],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized_t1]);
        end
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_t1.nii'];
        if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
            ea_maskimg(options,[directory,options.prefs.prenii_unnormalized_t1],bprfx);
        end
        from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized_t1];
        weights(cnt)=1.25;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

% PD
if usepd && ~strcmp(options.primarytemplate,'_pd')
    if exist([directory,options.prefs.prenii_unnormalized_pd],'file')
        disp('Including PD data for (grey-matter) normalization');
        if ~includeatlas
        ea_coreg2images(options,[directory,options.prefs.prenii_unnormalized_pd],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.prenii_unnormalized_pd]);
        end
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_pd.nii'];
        if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
            ea_maskimg(options,[directory,options.prefs.prenii_unnormalized_pd],bprfx);
        end
        from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized_pd];
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
                ea_isolate_fa(options);
            end
        end
        if exist([directory,options.prefs.fa],'file') % check again since could have been built above
            if ~includeatlas % if includeatlas is set we can assume that images have been coregistered and skulstripped already
 %               ea_rocrop([directory,options.prefs.fa]);
                if exist([directory,options.prefs.fa],'file') % recheck if has been built.
                    ea_coreg2images(options,[directory,options.prefs.fa],[directory,options.prefs.prenii_unnormalized],[directory,options.prefs.fa2anat]);
                end
            end
        end
    end
    if exist([directory,options.prefs.fa2anat],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
            ea_maskimg(options,[directory,options.prefs.fa2anat],bprfx);
        end
        to{cnt}=[options.earoot,'templates',filesep,'mni_hires_fa.nii'];
        from{cnt}=[directory,bprfx,options.prefs.fa2anat];
        weights(cnt)=0.5;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

%masks=segmentall(from,options); % masks not yet implemented.

% The convergence criterion for the multivariate scenario is a slave to the last metric you pass on the ANTs command line.
to{cnt}=[options.earoot,'templates',filesep,'mni_hires',options.primarytemplate,'.nii'];
if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
    ea_maskimg(options,[directory,options.prefs.prenii_unnormalized],bprfx);
end
from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized];
weights(cnt)=1.5;
metrics{cnt}='MI';
cnt=cnt+1;

if includeatlas % append as last to make criterion converge on this one.
   to{cnt}=[options.earoot,'templates',filesep,'mni_hires_distal.nii'];
   from{cnt}=[directory,'anat_atlas.nii.gz'];
   weights(cnt)=1.5;
   metrics{cnt}='MI'; % could think about changing this to CC
   cnt=cnt+1;
end


% Do the coreg part for postoperative images:
try
    ea_coregmr(options,options.prefs.normalize.coreg);
end

function masks=segmentall(from,options)
directory=[fileparts(from{1}),filesep];
for fr=1:length(from)
    [~,fn,ext]=fileparts(from{fr});
    switch [fn,ext]
        case options.prefs.fa2anat
            if ~exist([directory,'tc2',options.prefs.prenii_unnormalized],'file')
            	ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
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
                ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
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


function ea_genbrainmask(options)
directory=[options.root,options.patientname,filesep];
ea_newseg(directory,options.prefs.prenii_unnormalized,0,options);
c1=ea_load_nii([directory,'c1',options.prefs.prenii_unnormalized]);
c2=ea_load_nii([directory,'c2',options.prefs.prenii_unnormalized]);
c3=ea_load_nii([directory,'c3',options.prefs.prenii_unnormalized]);
bm=c1;
bm.img=c1.img+c2.img+c3.img;
bm.fname=[directory,'brainmask.nii'];
bm.img=bm.img>0.5;
spm_write_vol(bm,bm.img);


function ea_maskimg(options,file,prefix)
directory=[options.root,options.patientname,filesep];
if ~exist([directory,'brainmask.nii'],'file')
	ea_genbrainmask(options);
end
[pth,fn,ext]=fileparts(file);
if ~exist([pth,filesep,prefix,fn,ext],'file')
    nii=ea_load_nii(file);
    bm=ea_load_nii([directory,'brainmask.nii']);
    nii.img=nii.img.*double(bm.img);
    nii.fname=[pth,filesep,prefix,fn,ext];
    spm_write_vol(nii,nii.img);
end

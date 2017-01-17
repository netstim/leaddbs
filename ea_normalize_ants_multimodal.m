function varargout=ea_normalize_ants_multimodal(options,includeatlas)
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
    varargout{1}='Advanced Normalization Tools (ANTs) SyN - multimodal (Avants 2008)';
    varargout{2}=1;
    return
end


usebrainmask=0;
uset1=1;
uset2=1;
usepd=1;
usefa=1;
if ~exist('includeatlas','var')
    includeatlas=0;
end

if ~includeatlas % second run from maget-brain segment
ea_checkcoregallmri(options,usebrainmask)
end



directory=[options.root,options.patientname,filesep];
options.coregmr.method = 'Coreg MRIs: ANTs'; % hard-code to ANTs for now here.
cnt=1;

if usebrainmask
    bprfx='b';
else
    bprfx='';
end


% T1
if uset1 && ~strcmp(options.primarytemplate,'t1')
    if exist([directory,options.prefs.prenii_unnormalized_t1],'file')
        disp('Including T1 data for (grey-matter) normalization');
        to{cnt}=[ea_space(options),'t1.nii'];
        if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
            ea_maskimg(options,[directory,options.prefs.prenii_unnormalized_t1],bprfx);
        end
        from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized_t1];
        weights(cnt)=1.25;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

% T2
if uset2 && ~strcmp(options.primarytemplate,'t2')
    if exist([directory,options.prefs.prenii_unnormalized],'file')
        disp('Including T1 data for (grey-matter) normalization');
        to{cnt}=[ea_space(options),'t2.nii'];
        if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
            ea_maskimg(options,[directory,options.prefs.prenii_unnormalized],bprfx);
        end
        from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized];
        weights(cnt)=1.25;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end


% PD
if usepd && ~strcmp(options.primarytemplate,'pd')
    if exist([directory,options.prefs.prenii_unnormalized_pd],'file')
        disp('Including PD data for (grey-matter) normalization');
        to{cnt}=[ea_space(options),'pd.nii'];
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
    if exist([directory,options.prefs.fa2anat],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        to{cnt}=[ea_space(options),'fa.nii'];
        from{cnt}=[directory,bprfx,options.prefs.fa2anat];
        weights(cnt)=0.5;
        metrics{cnt}='MI';
        cnt=cnt+1;
    end
end

%masks=segmentall(from,options); % masks not yet implemented.

% The convergence criterion for the multivariate scenario is a slave to the last metric you pass on the ANTs command line.

if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
    ea_maskimg(options,[directory,options.prefs.prenii_unnormalized],bprfx);
end

to{cnt}=[ea_space(options),options.primarytemplate,'.nii'];
from{cnt}=[directory,bprfx,options.prefs.prenii_unnormalized];
weights(cnt)=1.5;
metrics{cnt}='MI';
cnt=cnt+1;

if includeatlas % append as last to make criterion converge on this one.
   to{cnt}=[ea_space(options),'distal.nii'];
   from{cnt}=[directory,'anat_atlas.nii.gz'];
   weights(cnt)=1.5;
   metrics{cnt}='MI'; % could think about changing this to CC
   cnt=cnt+1;
end



    ea_ants_nonlinear(to,from,[directory,options.prefs.gprenii],weights,metrics,options);
    ea_apply_normalization(options);

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
            masks{fr,1}=[ea_space(options),'c2mask.nii'];
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
            masks{fr,1}=[ea_space(options),'c1mask.nii'];
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
bm.img=bm.img>0.6;
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

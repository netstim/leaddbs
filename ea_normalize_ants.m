function varargout=ea_normalize_ants(options,includeatlas)
% This is a function that normalizes both a copy of transversal and coronal
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
    varargout{1}='ANTs (Avants 2008)';
    varargout{2}=1; % dummy output
    varargout{3}=1; % hassettings.
    varargout{4}=1; % is multispectral
    return
end

usefa=options.prefs.machine.normsettings.ants_usefa;
usebrainmask=0;

if ~exist('includeatlas','var')
    includeatlas=0;
end

if ~includeatlas % second run from maget-brain segment
 %   ea_checkcoregallmri(options,usebrainmask)
end

directory=[options.root,options.patientname,filesep];
cnt=1;
if usebrainmask
    bprfx='b';
else
    bprfx='';
end

spacedef=ea_getspacedef; % get definition of current space we are working in
[~,anatpresent]=ea_assignpretra(options);

if usefa && spacedef.hasfa % first put in FA since least important (if both an FA template and an fa2anat file is available)

    if exist([directory,options.prefs.fa2anat],'file') % recheck if now is present.
        disp('Including FA information for white-matter normalization.');
        to{cnt}=[ea_space(options),'fa.nii'];
        from{cnt}=[directory,bprfx,options.prefs.fa2anat];
        weights(cnt)=0.5;
        cnt=cnt+1;
    elseif exist([directory,options.prefs.fa],'file') % recheck if now is present.
            disp('Including FA information for white-matter normalization (weight = 0.5).');
            ea_coregimages(options,[directory,bprfx,options.prefs.fa],[directory,anatpresent{1}],[directory,bprfx,options.prefs.fa2anat],{},0,[],1);
            to{cnt}=[ea_space(options),'fa.nii'];
            from{cnt}=[directory,bprfx,options.prefs.fa2anat];
            weights(cnt)=0.5;
            cnt=cnt+1;
    end
end

anatpresent=flip(anatpresent); % reverse order since most important transform should be put in last.
% The convergence criterion for the multivariate scenario is a slave to the last metric you pass on the ANTs command line.
for anatf=1:length(anatpresent)
    disp(['Including ',anatpresent{anatf},' data for (grey-matter) normalization (weight = 1.25)']);

    to{cnt}=ea_niigz([ea_space(options),ea_det_to(anatpresent{anatf},spacedef)]);
    if usebrainmask && (~includeatlas) % if includeatlas is set we can assume that images have been coregistered and skulstripped already
        ea_maskimg(options,[directory,anatpresent{anatf}],bprfx);
    end
    from{cnt}=[directory,bprfx,anatpresent{anatf}];

        weights(cnt)=1.25;

    cnt=cnt+1;
end

if exist([directory,'segmentations'],'dir')
    segs=dir([directory,'segmentations']);
    for seg=1:length(segs)
        if ~strcmp(segs(seg).name(1),'.')
            if strfind(segs(seg).name,'.nii')
                if exist(ea_niigz([ea_space,'segmentations',filesep,segs(seg).name]),'file') % check if matching template exists
                    disp(['Including segmentations/',segs(seg).name,' for segment based assistance (weight = 20).']);
                    from=[{[directory,'segmentations',filesep,segs(seg).name]},from]; % append to front (since last one is convergence critical)
                    to=[{ea_niigz([ea_space,'segmentations',filesep,segs(seg).name])},to];
                    weights=[20,weights]; % set weight to 20 - DO NOT CHANGE THIS VALUE BELOW 3. IF VALUE IS CHANGED, SEGMENTATIONS WILL BE CONSIDERED SLABS IN ea_ants_nonlinear ~line 63 - would need to be changed there, as well.
                end
            end
        end
    end
end

if exist([directory,'fiducials'],'dir')
    segs=dir([directory,'fiducials',filesep,'native']);
    for seg=1:length(segs)
        if ~strcmp(segs(seg).name(1),'.')
            if strfind(segs(seg).name,'.nii')
                if exist(ea_niigz([directory,'fiducials',filesep,ea_getspace,filesep,segs(seg).name]),'file') % check if matching template exists
                    disp(['Including fiducials/',segs(seg).name,' for fiducial based assistance (weight = 5).']);
                    from=[{[directory,'fiducials',filesep,'native',filesep,segs(seg).name]},from]; % append to front (since last one is convergence critical)
                    to=[{ea_niigz([directory,'fiducials',filesep,ea_getspace,filesep,segs(seg).name])},to];
                    weights=[5,weights]; % set weight to 5 - DO NOT CHANGE THIS VALUE BELOW 3. IF VALUE IS CHANGED, FIDUCIALS WILL BE CONSIDERED SLABS IN ea_ants_nonlinear ~line 63 - would need to be changed there, as well.
                end
            end
        end
    end
end

if includeatlas % append as last to make criterion converge on this one.
    to{cnt}=ea_niigz([ea_space(options),'atlas']);
    from{cnt}=ea_niigz([directory,'anat_atlas']);
    weights(cnt)=1.5;
    cnt=cnt+1;
end

ea_ants_nonlinear(to,from,[directory,options.prefs.gprenii],weights,options);
ea_apply_normalization(options);

%% add methods dump:
if options.prefs.machine.normsettings.ants_scrf
    [scit,lcit]=ea_getspacedefcit;
    cits={
        'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'
        'Schoenecker, T., Kupsch, A., Kuehn, A. A., Schneider, G.-H., & Hoffmann, K. T. (2009). Automated Optimization of Subcortical Cerebral MR Imaging-Atlas Coregistration for Improved Postoperative Electrode Localization in Deep Brain Stimulation. AJNR Am J Neuroradiol, 30(10), 1914?1921. http://doi.org/10.3174/ajnr.A1741'
        };
    if ~isempty(lcit)
        cits=[cits;{lcit}];
    end

    ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',ea_cell2strlist(anatpresent),') using the'...
        ' SyN registration approach as implemented in Advanced Normalization Tools (Avants 2008; http://stnava.github.io/ANTs/).',...
        ' Nonlinear deformation into template space was achieved in five stages: After two linear (rigid followed by affine) steps, ',...
        ' A nonlinear (whole brain) SyN-registration stage was followed by two nonlinear SyN-registrations that consecutively focused on the area of interest ',...
        ' as defined by subcortical masks in Schoenecker 2008.'],...
        cits);
else
    [scit,lcit]=ea_getspacedefcit;
    cits={
        'Avants, B. B., Epstein, C. L., Grossman, M., & Gee, J. C. (2008). Symmetric diffeomorphic image registration with cross-correlation: evaluating automated labeling of elderly and neurodegenerative brain. Medical Image Analysis, 12(1), 26?41. http://doi.org/10.1016/j.media.2007.06.004'
        };
    if ~isempty(lcit)
        cits=[cits;{lcit}];
    end
    ea_methods(options,['Pre- (and post-) operative acquisitions were spatially normalized into ',ea_getspace,' space ',scit,' based on preoperative acquisition(s) (',ea_cell2strlist(anatpresent),') using the'...
        ' SyN registration approach as implemented in Advanced Normalization Tools (Avants 2008; http://stnava.github.io/ANTs/).',...
        ' Nonlinear deformation into template space was achieved in three stages: After two linear (rigid followed by affine) steps, ',...
        ' a nonlinear (whole brain) SyN registration stage was added.'],...
        cits);
end


function masks=segmentall(from,options)
directory=[fileparts(from{1}),filesep];
for fr=1:length(from)
    [~,fn,ext]=fileparts(from{fr});
    switch [fn,ext]
        case options.prefs.fa2anat
            if ~exist([directory,'tc2',options.prefs.prenii_unnormalized],'file')
            	ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0,options);
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
                ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0,options);
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
ea_newseg(fullfile(directory,options.prefs.prenii_unnormalized),0,options);
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


function templatetouse=ea_det_to(anatfile,spacedef)

anatfile=strrep(anatfile,'anat_','');
anatfile=ea_stripext(anatfile);

for avtpl=1:length(spacedef.templates)
    thistemp=spacedef.templates{avtpl};
    clear normmaptouse
    for normmap=1:size(spacedef.norm_mapping,1)
        thistocell=spacedef.norm_mapping{normmap,2};
        thisfromcell=spacedef.norm_mapping{normmap,1};
        if ~iscell(thistocell)
            thistocell={thistocell};
        end
        if ~iscell(thisfromcell)
            thisfromcell={thisfromcell};
        end
        if ismember(thistemp,thistocell)
            for fromfiles=1:length(thisfromcell)
                if contains(anatfile,thisfromcell{fromfiles})
                    templatetouse=thistemp;
                    return
                end
            end
        end
    end
end

% template still hasn't been assigned, use misfit template if not empty:
if ~isempty(spacedef.misfit_template)
    templatetouse=spacedef.misfit_template;
else
    templatetouse=spacedef.templates{1}; % just use first one.
end

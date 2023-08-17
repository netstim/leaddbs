function ea_ptspecific_atl(options)

if startsWith(options.atlasset, 'Local atlas: ') % manually installed atlas coded with this prefix by lead-dbs.
    return
end

troot=[options.earoot,'templates',filesep];
aroot=[ea_space(options,'atlases'),options.atlasset,filesep];
proot=[options.root,options.patientname,filesep];

if ~exist(ea_niigz([proot,'atlases',filesep,options.atlasset,filesep,'gm_mask']),'file') % check rebuild needed
    switch options.prefs.normalize.inverse.warp
        case 'tpm'
            generate_local_tpm(troot,aroot,proot,0,options)
        case 'inverse'
            ea_warp_atlas_to_native(troot,aroot,proot,0,options)
            load([proot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat']); % generates atlases variable
            options.atl.can=0;
            options.atl.ptnative=1;
            ea_genatlastable(atlases,[proot,'atlases',filesep],options);
    end
end

function ea_warp_atlas_to_native(troot,aroot,proot,force,options)

if ~exist([aroot,'atlas_index.mat'],'file')
    ea_error('Please visualize this atlas in MNI space once before visualizing the atlas in native space.');
else
    load([aroot,'atlas_index.mat']);
end

if ~exist([proot,'atlases'], 'dir')
    mkdir([proot,'atlases']);
end

copyfile(aroot, [proot,'atlases',filesep,options.atlasset]);
p=load([proot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat']);
p.atlases.rebuild=1;
save([proot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'], '-struct', 'p', '-v7.3');

ea_delete([proot,'atlases',filesep,options.atlasset,filesep,'gm_mask.nii*']);

if ismember(options.prefs.dev.profile,{'se'})
    interp=0;
else
    interp=1;
end

for atlas=1:length(atlases.names)
    subfolders = {};
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            subfolders{end+1} = 'lh';
        case 2 % right hemispheric atlas.
            subfolders{end+1} = 'rh';
        case 3 % both-sides atlas composed of 2 files.
            subfolders{end+1} = 'lh';
            subfolders{end+1} = 'rh';
        case 4 % mixed atlas (one file with both sides information.
            subfolders{end+1} = 'mixed';
        case 5 % midline atlas (one file with both sides information.
            subfolders{end+1} = 'midline';
    end
    
    for i = 1:length(subfolders)
        file = fullfile(proot, 'atlases', options.atlasset, subfolders{i}, atlases.names{atlas});
        
        if isnumeric(atlases.pixdim{atlas,1})
            ea_apply_normalization_tofile(options, {ea_niigz(file)}, {ea_niigz(file)}, 1, interp);
            ea_crop_nii(file);
            
        elseif ischar(atlases.pixdim{atlas,1})
            fib_load = load(file);
            src = fullfile(ea_space, 't1.nii');
            dest = options.subj.coreg.anat.preop.(options.subj.AnchorModality);
            transform = fullfile(options.subj.subjDir,'forwardTransform');
            
            switch atlases.pixdim{atlas,1}
                case 'fibers'
                    XYZ_mm = fib_load.fibers(:,1:3);
                    [~, XYZ_vox] = ea_map_coords([XYZ_mm'; ones(1,size(XYZ_mm,1))], src);
                    fib_load.fibers(:,1:3) = ea_map_coords(XYZ_vox, src, transform, dest)';
                case 'discfibers'
                    for side = 1:size(fib_load.fibcell,2)
                        XYZ_mm = vertcat(fib_load.fibcell{1,side}{:});
                        [~, XYZ_vox] = ea_map_coords([XYZ_mm'; ones(1,size(XYZ_mm,1))], src);
                        mapped = ea_map_coords(XYZ_vox, src, transform, dest)';
                        start_idx = 1;
                        for k = 1:length(fib_load.fibcell{1,side})
                            stop_idx = start_idx + size(fib_load.fibcell{1,side}{k},1) - 1;
                            fib_load.fibcell{1,side}{k} = mapped(start_idx:stop_idx,:);
                            start_idx = stop_idx + 1;
                        end
                    end
                otherwise
                    warning(['Unrecognized pixdim for ' atlases.names{atlas}]);
                    delete(file);
            end
            
            save(file, '-struct', 'fib_load', '-v7.3');
            
        else
            warning(['Unrecognized pixdim for ' atlases.names{atlas}]);
            delete(file);
        end
    end
end


function generate_local_tpm(troot,aroot,proot,force,options)

% make directories in patient folder
mkdir([proot,'atlases']);
mkdir([proot,'atlases',filesep,'native']);
if exist([proot,'atlases',filesep,options.atlasset],'file')
    return
end

mkdir([proot,'atlases',filesep,options.atlasset]);
mkdir([proot,'atlases',filesep,options.atlasset,filesep,'lh']);
mkdir([proot,'atlases',filesep,options.atlasset,filesep,'rh']);
mkdir([proot,'atlases',filesep,options.atlasset,filesep,'mixed']);
mkdir([proot,'atlases',filesep,options.atlasset,filesep,'midline']);

if ~exist([aroot,'atlas_index.mat'],'file')
    atlases=ea_genatlastable([],ea_space(options,'atlases'),options);
else
    load([aroot,'atlas_index.mat']);
end

cnt=1;

for atlas=1:length(atlases.names)
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            atlf=[aroot,'lh',filesep];
            patlf=[proot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            tpmf=[aroot,'tpm',filesep,'lh',filesep];
        case 2 % right hemispheric atlas.
            atlf=[aroot,'rh',filesep];
            patlf=[proot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];
            tpmf=[aroot,'tpm',filesep,'rh',filesep];
        case 3 % both-sides atlas composed of 2 files.
            ratlf=[aroot,'rh',filesep];
            pratlf=[proot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];

            latlf=[aroot,'lh',filesep];
            platlf=[proot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            rtpmf=[aroot,'tpm',filesep,'rh',filesep];
            ltpmf=[aroot,'tpm',filesep,'lh',filesep];
        case 4 % mixed atlas (one file with both sides information.
            atlf=[aroot,'mixed',filesep];
            patlf=[proot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep];
            tpmf=[aroot,'tpm',filesep,'mixed',filesep];
        case 5 % midline atlas (one file with both sides information.
            atlf=[aroot,'midline',filesep];
            patlf=[proot,'atlases',filesep,options.atlasset,filesep,'midline',filesep];
            tpmf=[aroot,'tpm',filesep,'midline',filesep];
    end

    for side=detsides(atlases.types(atlas))
        if atlases.types(atlas)==3
            switch side
                case 1
                    atlf=ratlf;
                    patlf=pratlf;
                    tpmf=rtpmf;
                case 2
                    atlf=latlf;
                    patlf=platlf;
                    tpmf=ltpmf;
            end
        end

        % gzip support
        if strcmp(atlases.names{atlas}(end-2:end),'.gz')
            gunzip([atlf,atlases.names{atlas}]);
            atln=atlases.names{atlas}(1:end-3);
            wasgz(cnt)=1;
        else
            atln=atlases.names{atlas};
            wasgz(cnt)=0;
        end

        if options.prefs.normalize.inverse.customtpm
            nii=ea_load_nii([atlf,atln]);
            nii.img=double(nii.img);
            nii.img=nii.img/max(nii.img(:)); % max 1
            nii.fname=[atlf,'t',atln];
            spm_write_vol(nii,nii.img);
            clear nii
        end

        atlasfile{cnt}=[atlf,'t',atln,',1'];
        oatlasfile{cnt}=[atlf,atln];

        tatlasfile{cnt}=[atlf,'t',atln];
        rawatlasfile{cnt}=[atlf,atln];

        atlfname{cnt}=atln;
        atlaspath{cnt}=patlf;
        cnt=cnt+1;
    end
end % collecting files loop


%% generate TPM
if options.prefs.normalize.inverse.customtpm

    if ~exist([aroot,'TPM_Lorio_Draganski.nii'],'file') || ~exist([aroot,'TPM_Lorio_Draganski.nii.gz'],'file') || force % check for pre-built TPM

        matlabbatch{1}.spm.util.imcalc.input = [
                                                {[troot,'TPM_Lorio_Draganski.nii,1'];}
                                                atlasfile'];
        matlabbatch{1}.spm.util.imcalc.output = 'TPM_Lorio_Draganski.nii';
        matlabbatch{1}.spm.util.imcalc.outdir = {aroot};
        matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = -1;
        matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;

        jobs{1}=matlabbatch;

        spm_jobman('run',jobs);
        clear jobs matlabbatch

        tnii=ea_load_untouch_nii([troot,'TPM_Lorio_Draganski.nii']);
        anii=ea_load_untouch_nii([aroot,'TPM_Lorio_Draganski.nii']);

        tnii.img(:,:,:,1)=0.5*tnii.img(:,:,:,1)+0.5*anii.img(:,:,:,1);
        tnii.img(:,:,:,2)=0.5*tnii.img(:,:,:,2)-0.5*anii.img(:,:,:,1);
        c1=tnii.img(:,:,:,1); c2=tnii.img(:,:,:,2);
        c1(c1>1)=1; c2(c2<0)=0;
        tnii.img(:,:,:,1)=c1; clear('c1'); tnii.img(:,:,:,2)=c2; clear('c2');
        ea_save_untouch_nii(tnii,[aroot,'TPM_Lorio_Draganski.nii']);

    end

    %% apply deformation fields:
    if exist([aroot,'TPM_Lorio_Draganski.nii.gz'],'file')
        gunzip([aroot,'TPM_Lorio_Draganski.nii.gz']);
    end

    tpmroot=aroot;

    matlabbatch{1}.spm.tools.preproc8.channel.vols = {[proot,options.prefs.prenii_unnormalized,',1']};
    matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
    matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
    matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,1']};
    matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,2']};
    matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,3']};
    matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,4']};
    matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,5']};
    matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[tpmroot,'TPM_Lorio_Draganski.nii,6']};
    matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
    matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
    matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
    matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
    matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear matlabbatch jobs

    gzip([aroot,'TPM_Lorio_Draganski.nii']);
    ea_delete([aroot,'TPM_Lorio_Draganski.nii']);
    movefile([proot,'iy_',options.prefs.prenii_unnormalized],[proot,'atlases',filesep,options.atlasset,filesep,'iy_warp.nii'])
    movefile([proot,'y_',options.prefs.prenii_unnormalized],[proot,'atlases',filesep,options.atlasset,filesep,'y_warp.nii'])
    warpfile=[proot,'atlases',filesep,options.atlasset,filesep,'iy_warp.nii'];
else
    warpfile=[proot,'y_ea_inv_normparams.nii'];
end

% warp atlas to patient space
for fi=1:length(oatlasfile)
    matlabbatch{1}.spm.util.defs.comp{1}.def = {warpfile};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = oatlasfile(fi);
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {atlaspath{fi}};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    jobs{1}=matlabbatch;
    spm_jobman('run',jobs);
    clear matlabbatch jobs
    movefile([atlaspath{fi},'w',atlfname{fi}],[atlaspath{fi},atlfname{fi}]);
    ea_crop_nii([atlaspath{fi},atlfname{fi}]);
end

% cleanup loop
for fi=1:length(atlasfile)
    if wasgz(fi)
        gzip(rawatlasfile{fi});
        ea_delete(rawatlasfile{fi});
    end
    if options.prefs.normalize.inverse.customtpm
        ea_delete(tatlasfile{fi});
    end
end


function sides=detsides(opt)

switch opt
    case 1 % left hemispheric atlas
        sides=1;
    case 2 % right hemispheric atlas
        sides=2;
    case 3
        sides=1:2;
    case 4
        sides=1:2;
    case 5
        sides=1; % midline
end

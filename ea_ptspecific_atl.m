%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

function ea_ptspecific_atl(options)

troot=[options.earoot,'templates',filesep];
aroot=[options.earoot,'atlases',filesep,options.atlasset,filesep];
proot=[options.root,options.patientname,filesep];
if ~exist([options.earoot,'templates',filesep,'TPM.nii'],'file')
   ea_generate_tpm;

end

generate_tpm(troot,aroot,proot,0,options)



function generate_tpm(troot,aroot,proot,force,options)

% make directories in patient folder
mkdir([proot,'atlases']);
mkdir([proot,'atlases',filesep,'native']);
if exist([proot,'atlases',filesep,'native',filesep,options.atlasset],'file')
    return
end
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset]);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'mixed']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'midline']);



if ~exist([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat'],'file')
    atlases=ea_genatlastable([],root,options);
else
    load([options.earoot,'atlases',filesep,options.atlasset,filesep,'atlas_index.mat']);
    atlases=ea_genatlastable(atlases,options.earoot,options);
end

cnt=1;
    
for atlas=1:length(atlases.names)
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            patlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'lh',filesep];
        case 2 % right hemispheric atlas.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];
            patlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'rh',filesep];
        case 3 % both-sides atlas composed of 2 files.
            ratlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];
            pratlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh',filesep];

            latlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            platlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh',filesep];
            rtpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'rh',filesep];
            ltpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'lh',filesep];
        case 4 % mixed atlas (one file with both sides information.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep];
            patlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'mixed',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'mixed',filesep];
        case 5 % midline atlas (one file with both sides information.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'midline',filesep];
            patlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'midline',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'midline',filesep];
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
    
    if ~exist([aroot,'TPM.nii'],'file') || ~exist([aroot,'TPM.nii.gz'],'file') || force % check for pre-built TPM
        
        
        
        
        matlabbatch{1}.spm.util.imcalc.input = [{
            [troot,'TPM.nii,1'];
            }
            atlasfile'];
        matlabbatch{1}.spm.util.imcalc.output = ['TPM.nii'];
        matlabbatch{1}.spm.util.imcalc.outdir = {aroot};
        matlabbatch{1}.spm.util.imcalc.expression = 'sum(X)';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = -1;
        matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        
        
        jobs{1}=matlabbatch;
        
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        
        
        
        tnii=ea_load_untouch_nii([troot,'TPM.nii']);
        anii=ea_load_untouch_nii([aroot,'TPM.nii']);
        
        tnii.img(:,:,:,1)=0.5*tnii.img(:,:,:,1)+0.5*anii.img(:,:,:,1);
        tnii.img(:,:,:,2)=0.5*tnii.img(:,:,:,2)-0.5*anii.img(:,:,:,1);
        c1=tnii.img(:,:,:,1); c2=tnii.img(:,:,:,2);
        c1(c1>1)=1; c2(c2<0)=0;
        tnii.img(:,:,:,1)=c1; clear('c1'); tnii.img(:,:,:,2)=c2; clear('c2');
        ea_save_untouch_nii(tnii,[aroot,'TPM.nii']);
        
    end
    
    
    
    %% apply deformation fields:
    
    if exist([aroot,'TPM.nii.gz'],'file')
        gunzip([aroot,'TPM.nii.gz']);
    end
    
    tpmroot=aroot;
else
    tpmroot=troot;
end

matlabbatch{1}.spm.tools.preproc8.channel.vols = {[proot,options.prefs.prenii_unnormalized,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[tpmroot,'TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[tpmroot,'TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[tpmroot,'TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[tpmroot,'TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[tpmroot,'TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[tpmroot,'TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
jobs{1}=matlabbatch;
cfg_util('run',jobs);
clear matlabbatch jobs
if options.prefs.normalize.inverse.customtpm
    gzip([aroot,'TPM.nii']);
    delete([aroot,'TPM.nii']);
    movefile([proot,'iy_',options.prefs.prenii_unnormalized],[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'iy_warp.nii'])
    movefile([proot,'y_',options.prefs.prenii_unnormalized],[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'y_warp.nii'])
    warpfile=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'iy_warp.nii'];
else
    warpfile=[proot,'y_ea_inv_normparams.nii'];
end
%apply deformation fields to respective atlas.

% warp atlas to patient space
for fi=1:length(oatlasfile)
    matlabbatch{1}.spm.util.defs.comp{1}.def = {warpfile};
    matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = oatlasfile(fi);
    matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {atlaspath{fi}};
    matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
    matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
    matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs
    movefile([atlaspath{fi},'w',atlfname{fi}],[atlaspath{fi},atlfname{fi}]);
    ea_crop_nii([atlaspath{fi},atlfname{fi}]);
end


% cleanup loop

for fi=1:length(atlasfile)
    if wasgz(fi)
        gzip(rawatlasfile{fi});
        delete(rawatlasfile{fi});
    end
    if options.prefs.normalize.inverse.customtpm
        delete(tatlasfile{fi});
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

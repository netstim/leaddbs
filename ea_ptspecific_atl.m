%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

function ea_ptspecific_atl(options)

troot=[options.earoot,'templates',filesep,'segment',filesep];
aroot=[options.earoot,'atlases',filesep,options.atlasset,filesep];
proot=[options.root,options.patientname,filesep];
ages={'E','M','Y'};
sides={'R','L'};


% first generate TPMs for each atlas.
generate_tpm(troot,aroot,options)


for age=1:length(ages)
    copyfile([proot,'pre_tra.nii'],[proot,ages{age},'_pre_tra.nii'])
    
    matlabbatch{1}.spm.tools.preproc8.channel.vols = {[proot,ages{age},'_pre_tra.nii,1']};
    matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
    matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
    matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[troot,'STN_',ages{age},'_TPM.nii,1']};
    matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [1 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[troot,'STN_',ages{age},'_TPM.nii,2']};
    matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[troot,'STN_',ages{age},'_TPM.nii,3']};
    matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[troot,'STN_',ages{age},'_TPM.nii,4']};
    matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[troot,'STN_',ages{age},'_TPM.nii,5']};
    matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[troot,'STN_',ages{age},'_TPM.nii,6']};
    matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
    matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
    matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
    matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
    matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
    matlabbatch{1}.spm.tools.preproc8.warp.write = [1 1];
    
    
    jobs{1}=matlabbatch;
    %cfg_util('run',jobs);
    clear matlabbatch jobs
    
    mkdir([proot,'atlases']);
    mkdir([proot,'atlases',filesep,'mixed']);
    matlabbatch{1}.spm.util.defs.comp{1}.def = {[proot,'y_',ages{age},'_pre_tra.nii']};
    matlabbatch{1}.spm.util.defs.ofname = '';
    matlabbatch{1}.spm.util.defs.fnames = {
        [troot,'tSTN_L_',ages{age},'.nii,1']
        [troot,'tSTN_R_',ages{age},'.nii,1']
        };
    matlabbatch{1}.spm.util.defs.savedir.saveusr = {proot};
    matlabbatch{1}.spm.util.defs.interp = 4;
    
    jobs{1}=matlabbatch;
    cfg_util('run',jobs);
    clear matlabbatch jobs
    
    
    
end





function generate_tpm(troot,aroot,options)


mkdir([aroot,'tpm']);
mkdir([aroot,'tpm','lh']);
mkdir([aroot,'tpm','rh']);
mkdir([aroot,'tpm','mixed']);
mkdir([aroot,'tpm','midline']);


atlases=ea_genatlastable(options);

for atlas=1:length(atlases.names)
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'lh',filesep];
        case 2 % right hemispheric atlas.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'rh',filesep];
        case 3 % both-sides atlas composed of 2 files.
            ratlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'lh',filesep];
            latlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'rh',filesep];
            rtpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'rh',filesep];
            ltpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'lh',filesep];
        case 4 % mixed atlas (one file with both sides information.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'mixed',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'mixed',filesep];
        case 5 % midline atlas (one file with both sides information.
            atlf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'midline',filesep];
            tpmf=[options.earoot,'atlases',filesep,options.atlasset,filesep,'tpm',filesep,'midline',filesep];
    end
    
    for side=detsides(atlases.types(atlas))
        if atlases.types(atlas)==3
        switch side
            case 1
                atlf=ratlf;
                tpmf=rtpmf;
            case 2
                atlf=latlf;
                tpmf=ltpmf;
        end
        end
        
        
        % gzip support
        
        if strcmp(atlases.names{atlas}(end-2:end),'.gz')
            gunzip([atlf,atlases.names{atlas}]);
            atlases.names{atlas}=atlases.names{atlas}(1:end-3);
            wasgz=1;
        else
            wasgz=0;
        end
        
        nii=load_untouch_nii([atlf,atlases.names{atlas}]);
        nii.img=double(nii.img);
        nii.img=nii.img/max(nii.img(:)); % max 1
        save_untouch_nii(nii,[atlf,'t',atlases.names{atlas}]);
        clear nii
        
        matlabbatch{1}.spm.util.imcalc.input = {
            [troot,'TPM.nii,1']
            [atlf,'t',atlases.names{atlas},',1']
            };
        matlabbatch{1}.spm.util.imcalc.output = ['t',atlases.names{atlas}];
        matlabbatch{1}.spm.util.imcalc.outdir = {tpmf};
        matlabbatch{1}.spm.util.imcalc.expression = 'i2';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = -4;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 16;
        
        
        jobs{1}=matlabbatch;
        cfg_util('run',jobs);
        clear jobs matlabbatch
        
        
        
    end
    tnii=load_untouch_nii([troot,'TPM.nii']);
    anii=load_untouch_nii([tpmf,'t',atlases.names{atlas}]);
    gzip([tpmf,'t',atlases.names{atlas}]);
    tnii.img(:,:,:,1)=tnii.img(:,:,:,1)+lnii.img+rnii.img;
    save_untouch_nii(tnii,[tpmf,'t',atlases.names{atlas}]);
    
    
    keyboard
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

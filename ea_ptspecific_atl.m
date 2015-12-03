%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 4252 $)
%-----------------------------------------------------------------------

function ea_ptspecific_atl(options)

troot=[options.earoot,'templates',filesep,'segment',filesep];
aroot=[options.earoot,'atlases',filesep,options.atlasset,filesep];
proot=[options.root,options.patientname,filesep];
if ~exist([options.earoot,'templates',filesep,'TPM.nii'],'file')
   ea_generate_tpm;

end

generate_tpm(troot,aroot,proot,1,options)








function generate_tpm(troot,aroot,proot,force,options)

% make directories in tpm folder
mkdir([aroot,'tpm']);
mkdir([aroot,'tpm',filesep,'lh']);
mkdir([aroot,'tpm',filesep,'rh']);
mkdir([aroot,'tpm',filesep,'mixed']);
mkdir([aroot,'tpm',filesep,'midline']);

% make directories in patient folder
mkdir([proot,'atlases']);
mkdir([proot,'atlases',filesep,'native']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset]);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'mixed']);
mkdir([proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'midline']);

atlases=ea_genatlastable(options.earoot,options);

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
            wasgz=1;
        else
            atln=atlases.names{atlas};
            wasgz=0;
        end

        if ~exist([tpmf,atlases.names{atlas}],'file') || force % check for pre-built TPM


            nii=load_untouch_nii([atlf,atln]);
            nii.img=double(nii.img);
            nii.img=nii.img/max(nii.img(:)); % max 1
            save_untouch_nii(nii,[atlf,'t',atln]);
            clear nii

            matlabbatch{1}.spm.util.imcalc.input = {
                [troot,'TPM.nii,1'];
                [atlf,'t',atln,',1']
                };
            matlabbatch{1}.spm.util.imcalc.output = [atln];
            matlabbatch{1}.spm.util.imcalc.outdir = {tpmf};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1+i2';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = -4;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 16;


            jobs{1}=matlabbatch;

            cfg_util('run',jobs);
            clear jobs matlabbatch
            delete([atlf,'t',atln]);




            tnii=load_untouch_nii([troot,'TPM.nii']);
            anii=load_untouch_nii([tpmf,atln]);

            tnii.img(:,:,:,1)=tnii.img(:,:,:,1)+anii.img(:,:,:,1);
            save_untouch_nii(tnii,[tpmf,atln]);

        end

        if ~exist([patlf,atlases.names{atlas}],'file') || force

            % tpm file for atlas is now in folder. Begin segmentation..

            % gzip support..
            if ~exist([tpmf,atln],'file') && exist([tpmf,atln,'.gz'],'file')
                gunzip([tpmf,atln,'.gz']);
                wasgz=1;
            end


            matlabbatch{1}.spm.tools.preproc8.channel.vols = {[proot,options.prefs.prenii_unnormalized,',1']};
            matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
            matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
            matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[tpmf,atln,',1']};
            matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[tpmf,atln,',2']};
            matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[tpmf,atln,',3']};
            matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[tpmf,atln,',4']};
            matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[tpmf,atln,',5']};
            matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 2;
            matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
            matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[tpmf,atln,',6']};
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



            %apply deformation fields to respective atlas.

            % warp atlas to patient space

            matlabbatch{1}.spm.util.defs.comp{1}.def = {[proot,'iy_',options.prefs.prenii_unnormalized]};
            matlabbatch{1}.spm.util.defs.ofname = '';
            matlabbatch{1}.spm.util.defs.fnames = {[atlf,atln,',1']};
            matlabbatch{1}.spm.util.defs.savedir.saveusr = {patlf};
            matlabbatch{1}.spm.util.defs.interp = 1;
            jobs{1}=matlabbatch;
            cfg_util('run',jobs);
            clear matlabbatch jobs
            movefile([patlf,'w',atln],[patlf,atln]);


%             % apply original transformation back to warped atlas.
%             matlabbatch{1}.spm.util.defs.comp{1}.def = {[proot,'y_ea_normparams.nii']};
%             matlabbatch{1}.spm.util.defs.ofname = '';
%             matlabbatch{1}.spm.util.defs.fnames = {[patlf,atln,',1']};
%             matlabbatch{1}.spm.util.defs.savedir.saveusr = {patlf};
%             matlabbatch{1}.spm.util.defs.interp = 1;
%             jobs{1}=matlabbatch;
%             cfg_util('run',jobs);
%             clear matlabbatch jobs
%
%             movefile([patlf,'w',atln],[patlf,atln]);


            % re-gzip tpm, patl and atlas file.
            if wasgz
                try gzip([tpmf,atln]); end
                try delete([tpmf,atln]); end
                try gzip([atlf,atln]); end
                try delete([atlf,atln]); end
                try gzip([patlf,atln]); end
                try delete([patlf,atln]); end
            end
        end



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

function ea_norm_ptspecific_atl(options)

% troot=[options.earoot,'templates',filesep,'segment',filesep];
% aroot=[ea_space(options,'atlases'),options.atlasset,filesep];
proot=[options.root,options.patientname,filesep];

force=0; % always re-process..


% make directories in patient folder
mkdir([proot,'atlases']);
mkdir([proot,'atlases']);
if exist([proot,'atlases',filesep,'mni',filesep,options.atlasset],'file');
    return
end
mkdir([proot,'atlases',filesep,'mni',filesep,options.atlasset]);
mkdir([proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'lh']);
mkdir([proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'rh']);
mkdir([proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'mixed']);
mkdir([proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'midline']);


    if ~exist([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat'],'file')
        atlases=ea_genatlastable([],ea_space('atlases'),options);
    else
        load([ea_space(options,'atlases'),options.atlasset,filesep,'atlas_index.mat']);
        atlases=ea_genatlastable(atlases,ea_space('atlases'),options);
    end

for atlas=1:length(atlases.names)
    switch atlases.types(atlas)
        case 1 % left hemispheric atlas.
            natlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh',filesep];
            patlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'lh',filesep];
        case 2 % right hemispheric atlas.
            natlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh',filesep];
            patlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'rh',filesep];
        case 3 % both-sides atlas composed of 2 files.
            nratlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'rh',filesep];
            pratlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'rh',filesep];
            nlatlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'lh',filesep];
            platlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'lh',filesep];
        case 4 % mixed atlas (one file with both sides information.
            natlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'mixed',filesep];
            patlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'mixed',filesep];
        case 5 % midline atlas (one file with both sides information.
            natlf=[proot,'atlases',filesep,'native',filesep,options.atlasset,filesep,'midline',filesep];
            patlf=[proot,'atlases',filesep,'mni',filesep,options.atlasset,filesep,'midline',filesep];
    end

    for side=detsides(atlases.types(atlas))
        if atlases.types(atlas)==3
            switch side
                case 1
                    natlf=nratlf;
                    patlf=pratlf;
                case 2
                    natlf=nlatlf;
                    patlf=platlf;
            end
        end
        % gzip support

        if strcmp(atlases.names{atlas}(end-2:end),'.gz')
            try
                gunzip([natlf,atlases.names{atlas}]);
            end
            atln=atlases.names{atlas}(1:end-3);
            wasgz=1;
        else
            atln=atlases.names{atlas};
            wasgz=0;
        end



        if ~exist([patlf,atlases.names{atlas}],'file') || force

            % apply original transformation back to warped atlas.

            matlabbatch{1}.spm.util.defs.comp{1}.def = {[proot,'y_ea_normparams.nii']};
            matlabbatch{1}.spm.util.defs.out{1}.pull.fnames = {[natlf,atln]};
            matlabbatch{1}.spm.util.defs.out{1}.pull.savedir.saveusr = {patlf};
            matlabbatch{1}.spm.util.defs.out{1}.pull.interp = 4;
            matlabbatch{1}.spm.util.defs.out{1}.pull.mask = 1;
            matlabbatch{1}.spm.util.defs.out{1}.pull.fwhm = [0 0 0];
            jobs{1}=matlabbatch;
            spm_jobman('run',jobs);
            clear matlabbatch jobs

            movefile([patlf,'w',atln],[patlf,atln]);


            % re-gzip tpm, patl and atlas file.
            if wasgz
                try gzip([natlf,atln]); end
                ea_delete([natlf,atln]);
                try gzip([patlf,atln]); end
                ea_delete([patlf,atln]);
            end
        end



    end




end

%% methods dump:
ea_methods(directory,...
            ['Subcortical atlases were projected into native space using the inverse deformation field mapping from native to template space (estimated in the normalization step based on pre-operative acquisitions) ',...
            'as implemented in Lead-DBS software',...
            ' (Horn & Kuehn 2005; SCR_002915; http://www.lead-dbs.org).'],...
            {'Horn, A., & Kühn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});






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

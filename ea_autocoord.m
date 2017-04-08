function ea_autocoord(options)
% This function is the main function of LEAD-DBS. It will generate a
% vector of coordinates.
% Trajectory{1} will be the right trajectory, trajectory{2} the
% left one.
% For each hemisphere of the brain, this function will call the
% reconstruction routine ea_autocoord_side and lateron call functions for
% manual correction of the results, render and slice views.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn

% get accurate electrode specifications and save it in options.
options=ea_resolve_elspec(options);
options.prefs=ea_prefs(options.patientname);

directory = [options.root,options.patientname,filesep];

if options.dicomimp % do DICOM-Import.
    ea_dicom_import(options);
    return
end

% check connectome-mapper tags
if isfield(options,'lcm')
        ea_lcm(options);
end

if ~strcmp(options.patientname,'No Patient Selected') % only 3D-rendering viewer can be opened if no patient is selected.

    % move files for compatibility
    try  ea_compat_patfolder(options); end

    % assign/order anatomical images
    [options,presentfiles]=ea_assignpretra(options);

    % generate grid file
    if ~exist(ea_niigz([directory,'grid.nii']),'file')
        try
            ea_gengrid(options);
        end
    end

    % anat preprocess, only do once.
    % a small hidden file '.pp' inside patient folder will show this has been done before.
    if ~exist([directory,'.pp'],'file') && ~exist([directory,'ea_normmethod_applied.mat'],'file')
        % apply reorientation/cropping and biasfieldcorrection
        for fi=1:length(presentfiles)
            ea_anatpreprocess([directory,presentfiles{fi}]);
        end

        % Reslice(interpolate) preoperative anatomical image if needed
try        ea_resliceanat(options); end

        try
            fs = fopen([directory,'.pp'],'w');
            fprintf(fs,'%s','anat preprocess done');
            fclose(fs);
        end
    end

    if options.modality==2 % CT support
        options.prefs.tranii=options.prefs.ctnii;
        options.prefs.tranii_unnormalized=options.prefs.rawctnii_unnormalized;
    end

    if options.coregct.do
        eval([options.coregct.method,'(options)']); % triggers the coregct function and passes the options struct to it.
        ea_dumpnormmethod(options,options.coregct.method,'coregctmethod');
        ea_tonemapct_file(options,'native'); % (Re-) compute tonemapped (native space) CT
        ea_gencoregcheckfigs(options); % generate checkreg figures
    end

    if options.coregctcheck
        % export "control" niftis with wireframe of normal anatomy..
        ea_show_ctcoregistration(options);
    end


    if options.normalize.do

        % 1. coreg all available preop MRI
        ea_checkcoregallmri(options,0,1); % check and coregister all preoperative MRIs here.

        % 2. then coreg post to pre MRI:
        %try % fix me - can we get rid of this try/catch here?
            ea_coregmr(options);

        %end

        % 3. finally perform normalization based on dominant or all preop
        % MRIs:
        ea_dumpnormmethod(options,options.normalize.method,'normmethod'); % has to come first due to applynormalization.
        eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.

        if options.modality==2 % (Re-) compute tonemapped (normalized) CT
              ea_tonemapct_file(options,'mni');
        end
        % 4. generate coreg-check figs (all to all).
        ea_gencoregcheckfigs(options); % generate checkreg figures
    end

    if isfield(options,'gencheckreg') % this case is an exception when calling from the Tools menu.
        if options.gencheckreg
                    ea_gencoregcheckfigs(options); % generate checkreg figures
        end
    end

    if options.dolc % perform lead connectome subroutine..
        ea_perform_lc(options);
    end

    if options.atl.genpt % generate patient specific atlas set
        ea_ptspecific_atl(options);
    end

    if options.atl.normalize % normalize patient-specific atlas-set.
        ea_norm_ptspecific_atl(options)
    end

    if options.normalize.check
        % export "control" niftis with wireframe of normal anatomy..
        ea_show_normalization(options);
    end

    if options.doreconstruction
        ea_checkfiles(options);

        for side=1:length(options.sides)
            %try
            % call main routine reconstructing trajectory for one side.
            [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(options.patientname,options,options.sides(side));

            % refit electrodes starting from first electrode (this is redundant at this point).
            coords_mm{side} = ea_map_coords(coords', [directory,options.prefs.tranii])';

            [~,distmm]=ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),[directory,options.prefs.tranii]);

            comp = ea_map_coords([0,0,0;trajvector{side}]', [directory,options.prefs.tranii])'; % (XYZ_mm unaltered)

            trajvector{side}=diff(comp);

            normtrajvector{side}=trajvector{side}./norm(trajvector{side});

            for electrode=2:4
                coords_mm{side}(electrode,:)=coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
            end
            markers(side).head=coords_mm{side}(1,:);
            markers(side).tail=coords_mm{side}(4,:);

            orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);

            markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
            markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality

            coords_mm=ea_resolvecoords(markers,options);
        end

        %     if length(coords_mm)==4 % only one side was processed.
        %         if options.sides==1
        %             coords_mm=[nan(4,3);coords_mm];
        %         elseif options.sides==2
        %             coords_mm=[coords_mm;nan(4,3)];
        %         end
        %     end

        % transform trajectory to mm space:
        for side=1:length(options.sides)
            try
                if ~isempty(trajectory{side})
                    trajectory{side}=ea_map_coords(trajectory{side}', [directory,options.prefs.tranii])';
                end

            end
        end
        % save reconstruction results
        elmodel=options.elmodel;
        ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
        ea_methods(options,...
            ['DBS-Electrodes were automatically pre-localized in native & template space using Lead-DBS software',...
            ' (Horn & Kuehn 2005; SCR_002915; http://www.lead-dbs.org).'],...
            {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});
    end

    if options.manualheightcorrection
        % load reconstruction results
        % try
        %     [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
        % catch
        %     ea_error([patientname,': No reconstruction information found. Please run reconstruction first.']);
        % end
        % ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
        mcfig=figure('name',[options.patientname,': Manual Height Correction'],'numbertitle','off');
        warning('off');
        try
            ea_maximize(mcfig);
        end
        ea_manualreconstruction(mcfig,options.patientname,options);
    else
        ea_write(options)
    end

else
    ea_write(options)
end

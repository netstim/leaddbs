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

% set patientdir
options.prefs.patientdir = options.patientname;

% get accurate electrode specifications and save it in options.
options = ea_resolve_elspec(options);

directory = [options.root,options.patientname,filesep];

if options.dicomimp || options.assignnii % do DICOM-Import.
    if options.dicomimp
        if strcmp(options.patientname, 'No Patient Selected')
            msgbox('Please choose patient directory first!','Error','error');
        else
            ea_dicom_import(options);
        end
    end

    if options.assignnii
        if strcmp(options.patientname, 'No Patient Selected')
            msgbox('Please choose patient directory first!','Error','error');
        else
            outdir = [options.root, options.patientname, filesep];
            % assign image type here
            di = dir([outdir,'*.nii']);
            di = ea_sortbytes(di);
            for d=1:length(di)
                dcfilename=[outdir,di(d).name];
                ea_imageclassifier({dcfilename});
            end
            figs=allchild(0);
            ids={figs.Tag};
            [~,imclassfids]=ismember(ids,'imclassf');
            if ~any(imclassfids)
                msgbox('All NIfTI files have been assigned already.');
            else
                set(figs(logical(imclassfids)),'Visible','on');
            end
            if isempty(di)
                msgbox('Could not find any NIfTI files to rename/assign.');
            end
        end
    end

    return % For now we recommend to do import & processing in separate run calls.
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
        try ea_resliceanat(options); end

        try
            fs = fopen([directory,'.pp'],'w');
            fprintf(fs,'%s','anat preprocess done');
            fclose(fs);
        end
    end

    if options.modality == 2 % CT support
        options.prefs.tranii=options.prefs.ctnii;
        options.prefs.tranii_unnormalized=options.prefs.rawctnii_unnormalized;

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
    end

    if options.normalize.do
        % 1. coreg all available preop MRI
        ea_checkcoregallmri(options,0,1); % check and coregister all preoperative MRIs here.

        % 2. then coreg postop MRI to preop MRI
        ea_coregmr(options);

        % 3. finally perform normalization based on dominant or all preop MRIs
        ea_dumpnormmethod(options,options.normalize.method,'normmethod'); % has to come first due to applynormalization.
        eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.

        if options.modality == 2 % (Re-) compute tonemapped (normalized) CT
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
        if ~options.prefs.env.dev % hard set to TRAC/CORE if not in dev mode.
            options.reconmethod=1;
        end

        switch options.reconmethod
            case 1 % TRAC/CORE
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


                % transform trajectory to mm space:
                for side=1:length(options.sides)
                    try
                        if ~isempty(trajectory{side})
                            trajectory{side}=ea_map_coords(trajectory{side}', [directory,options.prefs.tranii])';
                        end

                    end
                end
                % save reconstruction results
                ea_methods(options,...
                    ['DBS-Electrodes were automatically pre-localized in native & template space using Lead-DBS software',...
                    ' (Horn & Kuehn 2015; SCR_002915; http://www.lead-dbs.org).'],...
                    {'Horn, A., & Kuehn, A. A. (2015). Lead-DBS: a toolbox for deep brain stimulation electrode localizations and visualizations. NeuroImage, 107, 127?135. http://doi.org/10.1016/j.neuroimage.2014.12.002'});

            case 2 % PaCER
                options.prefs.ctnii_coregistered

                elecmodels=PaCER([options.root,options.patientname,filesep,options.prefs.ctnii_coregistered],'finalDegree',1,'electrodeType',ea_mod2pacermod(options.elmodel));

                for side=options.sides
                    coords_mm{side}=elecmodels{side}.getContactPositions3D;
                    for dim=1:3
                        trajectory{side}(:,dim)=linspace(coords_mm{side}(1,dim),coords_mm{side}(1,dim)+10*(coords_mm{side}(1,dim)-coords_mm{side}(end,dim)),20);
                    end

                    markers(side).head=coords_mm{side}(1,:);
                    markers(side).tail=coords_mm{side}(4,:);
                    normtrajvector{side}=(coords_mm{side}(1,:)-coords_mm{side}(end,:))/...
                    norm((coords_mm{side}(1,:)-coords_mm{side}(end,:)));
                    orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);

                    markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
                    markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality

                end

                options.native=1;
                options.hybridsave=1;
                ea_methods(options,...
                    ['DBS-Electrodes were automatically pre-localized in native & template space using the PaCER algorithm',...
                    ' (Husch 20xx; http://adhusch.github.io/PaCER/).'],...
                    {'Husch (20xx). PaCER - A fully automated method for electrode trajectory and contact reconstruction in deep brain stimulation.'});
        end
        elmodel=options.elmodel;
        ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
        if isfield(options,'hybridsave')
            options=rmfield(options,'hybridsave');
        end
 
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
        %warning('off');
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

function di=ea_sortbytes(di)
if isempty(di)
    return
end
for d=1:length(di)
    bytesc(d)=di(d).bytes;
end
[~,order]=sort(bytesc,'ascend');
di=di(order);

function model=ea_mod2pacermod(model)
% current dictionary to translate between Lead-DBS and PaCER nomenclature.
% Hoping to standardize this in the future.
switch model
    case 'Medtronic 3389'
       % pass through (same nomenclature)
    case 'Medtronic 3387'
       % pass through (same nomenclature)
    case 'Boston Scientific Vercise Directed'
        model='Boston Vercise Directional';
    otherwise
        model='Unkown Electrode Type';
end

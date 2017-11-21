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

if isfield(options,'predict')
   ea_predict(options); 
end

if ~strcmp(options.patientname,'No Patient Selected') % only 3D-rendering viewer can be opened if no patient is selected.

    % move files for compatibility
    try ea_compat_patfolder(options); end

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

        try ea_precoreg([directory,presentfiles{1}],options.primarytemplate); end

        try
            fs = fopen([directory,'.pp'],'w');
            fprintf(fs,'%s','anat preprocess done');
            fclose(fs);
        end
    end

    % NEED FURTHER TUNE: auto detection of MRCT modality for the patient
    try
        modality = ea_checkctmrpresent(directory);
        modality = find(modality);
        if isempty(modality)    % no postop image present
            options.modality = 1;    % set to MR to work it around
        elseif length(modality) == 2    % both MR and CT image present
            options.modality = options.prefs.preferMRCT;  % set the modality according to 'prefs.preferMRCT'
        else    % only one modality present
            options.modality = modality;
        end
    end

    if options.modality == 2 % CT support
        options.prefs.tranii=options.prefs.ctnii;
        options.prefs.tranii_unnormalized=options.prefs.rawctnii_unnormalized;

        if options.coregct.do && ~ea_coreglocked(options,['tp_',options.prefs.ctnii_coregistered])
            diary([directory, 'coregCT_', datestr(now, 'yyyymmddTHHMMss'), '.log']);
            eval([options.coregct.method,'(options)']); % triggers the coregct function and passes the options struct to it.
            ea_dumpnormmethod(options,options.coregct.method,'coregctmethod');
            ea_tonemapct_file(options,'native'); % (Re-) compute tonemapped (native space) CT
            ea_gencoregcheckfigs(options); % generate checkreg figures
            diary off
        end

    end

    if options.coregmr.do
        diary([directory, 'coregMR_', datestr(now, 'yyyymmddTHHMMss'), '.log']);
        % 1. coreg all available preop MRI
        ea_checkcoregallmri(options,0,1); % check and coregister all preoperative MRIs here.

        % 2. then coreg postop MRI to preop MRI
        ea_coregmr(options);
        diary off
    end

    if options.coregmr.check
        options.normcoreg='coreg';
        ea_checkcoreg(options);
    end

    if options.normalize.do
        diary([directory, 'normalize_', datestr(now, 'yyyymmddTHHMMss'), '.log']);
        if ~(ea_coreglocked(options,'glanat')==2) || strcmp(options.normalize.method,'ea_normalize_apply_normalization') % =2 means permanent lock for normalizations and only happens if all preop anatomy files were approved at time of approving normalization.
            if ea_coreglocked(options,'glanat')==1 && ~strcmp(options.normalize.method,'ea_normalize_apply_normalization') % in this case, only perform normalization if using a multispectral approach now.
                [~,~,~,doit]=eval([options.normalize.method,'(''prompt'')']);
            else
                doit=1;
            end
            if doit || strcmp(options.normalize.method,'ea_normalize_apply_normalization')
                clear doit
                % 3. finally perform normalization based on dominant or all preop MRIs
                ea_dumpnormmethod(options,options.normalize.method,'normmethod'); % has to come first due to applynormalization.
                eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.

                if options.modality == 2 % (Re-) compute tonemapped (normalized) CT
                    ea_tonemapct_file(options,'mni');
                end
                % 4. generate coreg-check figs (all to all).
                ea_gencoregcheckfigs(options); % generate checkreg figures
            end
        end
        diary off
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
        options.normcoreg='normalize';
        ea_checkcoreg(options);

    end

    if options.doreconstruction
        switch options.reconmethod
            case 1 % TRAC/CORE
                [coords_mm,trajectory,markers]=ea_runtraccore(options);
                        options.native=0;

            case 2 % PaCER
                try
                    [coords_mm,trajectory,markers]=ea_runpacer(options);
                            options.native=1;

                catch % revert to TRAC/CORE
                    [coords_mm,trajectory,markers]=ea_runtraccore(options);
                            options.native=0;

                end
        end
        options.hybridsave=1;
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



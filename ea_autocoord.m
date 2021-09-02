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
options.tic=tic;

% set patientdir
options.prefs.patientdir = options.patientname;

% Get BIDS fetcher and subj struct
bids = getappdata(options.leadfigure, 'bids');
subjId = getappdata(options.leadfigure, 'subjId');
options.subj = bids.getSubj(subjId{1}, options.modality);

% get accurate electrode specifications and save it in options.
options = ea_resolve_elspec(options);
directory = [options.root,options.patientname,filesep];

if strcmp(options.leadprod, 'dbs') || strcmp(options.leadprod, 'connectome')
    if options.dicomimp.do || options.assignnii % do DICOM-Import.
        if options.dicomimp.do
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
end

% check connectome-mapper tags
if isfield(options,'lcm')
    ea_lcm(options);
end

if isfield(options,'predict')
   ea_predict(options);
end

% only 3D-rendering viewer can be opened if no patient is selected.
if ~strcmp(options.patientname,'No Patient Selected') && ~isempty(options.patientname)
    % Copy post-op images to preprocessing folder, no preproc is done for now
    fields = fieldnames(options.subj.postopAnat);
    for i=1:length(fields)
        if ~isfile(options.subj.postopAnat.(fields{i}).preproc)
            ea_mkdir(fileparts(options.subj.postopAnat.(fields{i}).preproc));
            copyfile(options.subj.postopAnat.(fields{i}).raw, options.subj.postopAnat.(fields{i}).preproc);
        end
    end

    % Preprocessing pre-op images
    preprocessing = 0;
    fields = fieldnames(options.subj.preopAnat);
    for i=1:length(fields)
        if ~isfile(options.subj.preopAnat.(fields{i}).preproc)
            % Copy files
            ea_mkdir(fileparts(options.subj.preopAnat.(fields{i}).preproc));
            copyfile(options.subj.preopAnat.(fields{i}).raw, options.subj.preopAnat.(fields{i}).preproc);

            % Run reorientation, cropping and bias field correction
            ea_anatpreprocess(options.subj.preopAnat.(fields{i}).preproc);

            % Preprocessing only for pre-op anchor image
            if i==1
                ea_resliceanat(options.subj.preopAnat.(fields{i}).preproc);
                % ea_acpcdetect(options.subj.preopAnat.(fields{i}).preproc);
            end

            preprocessing = 1;
        end
    end

    if preprocessing
        fprintf('\nPreprocessing finished.\n\n');
    end

    % Set primary template
    if ismember(options.subj.AnchorModality, fieldnames(bids.spacedef.norm_mapping))
        options.primarytemplate = bids.spacedef.norm_mapping.(options.subj.AnchorModality);
    else
        options.primarytemplate = bids.spacedef.misfit_template;
    end

    % Pre-coregister pre-op anchor image
    ea_precoreg(options.subj.preopAnat.(fields{1}).preproc, ... % Input anchor image
        options.primarytemplate, ... % Template to use
        options.subj.preopAnat.(fields{1}).coreg, ... % Output pre-coregistered image
        options.subj.coreg.transform.(fields{1})); % % Pre-coregistration transform

    % NEED FURTHER TUNE: auto detection of MRCT modality for the patient
    try
        options.modality = ea_getmodality(directory);
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

                % cleanup already normalized versions:
                ea_delete([options.root,options.patientname,filesep,options.prefs.gprenii]);
                for fi=2:length(presentfiles)
                    ea_delete([options.root,options.patientname,filesep,'gl',presentfiles{fi}]);
                end

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

    if options.atl.normalize % normalize patient's atlas-set.
        ea_norm_ptspecific_atl(options)
    end

    if options.scrf.do
        if ~ea_coreglocked(options,'brainshift') || options.overwriteapproved
            options.autobrainshift=1;
            ea_subcorticalrefine(options);
            options=rmfield(options,'autobrainshift');
        end
    end

    if options.normalize.check %check box "Check Results" in "Volume Registrations" panel
        % export "control" niftis with wireframe of normal anatomy..
        if ~exist([directory,'checkreg'],'file')
            ea_gencoregcheckfigs(options); % generate checkreg figures if they not yet exist
        end
        options.normcoreg='normalize';
        ea_checkcoreg(options);
        drawnow; % this prevents the figure from changing the name with multiple subjects
        e=evalin('base', 'checkregempty');
        evalin('base',' clear checkregempty');
        if e && ~ea_coreglocked(options,'brainshift') ...
             && exist([directory,'scrf',filesep,'scrf_instore.mat'], 'file')
            ea_subcorticalrefine(options);
        end
    end

    if options.normalize.refine
        if options.prefs.env.dev
            ea_runwarpdrive(options);
        else
            ea_checkstructures(options);
        end
    end

    if options.ecog.extractsurface.do
       switch options.ecog.extractsurface.method
           case 1 % CAT 12
               hastb=ea_hastoolbox('cat');
               if ~hastb
                   ea_error('CAT12 needs to be installed to the SPM toolbox directory');
               end
               ea_cat_seg(options);
           case 2 % FS
               hastb=ea_hastoolbox('freesurfer');

               if ~hastb
                   ea_error('Freesurfer needs to be installed and connected to Lead-DBS');
               end
               hastb=ea_hastoolbox('fsl');
               if ~hastb
                   ea_error('FSL needs to be installed and connected to Lead-DBS');
               end

               options.prefs=ea_prefs;
               [options,presentfiles]=ea_assignpretra(options);
               setenv('SUBJECTS_DIR',[options.root,options.patientname,filesep]);
               if exist([options.root,options.patientname,filesep,'fs'],'dir')
                   rmdir([options.root,options.patientname,filesep,'fs'],'s');
               end
               system([options.prefs.fspath,filesep,'bin',filesep,...
                   'recon-all',...
                   ' -subjid fs',...
                   ' -i ',[options.root,options.patientname,filesep,presentfiles{1}],...
                   ' -all']);
       end
    end

    if options.doreconstruction
        wasnative=options.native;
        poptions=ea_checkmanapproved(options);
        if ~isempty(poptions.sides)
            switch options.reconmethod
                case 'Refined TRAC/CORE' % refined TRAC/CORE
                    [coords_mm,trajectory,markers]=ea_runtraccore(poptions);
                    options.native = 0; % Output in template space
                    options.hybridsave=1; % save output of TRAC/CORE before progressing
                    options.elside=options.sides(1);
                    elmodel=options.elmodel;
                    ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
                    [coords_mm,trajectory,markers] = ea_refinecoords(poptions); % experimental fiducial marker refine method
                    options.native = 1;

                case 'TRAC/CORE (Horn 2015)' % TRAC/CORE
                    [coords_mm,trajectory,markers]=ea_runtraccore(poptions);
                    options.native=0;

                case 'PaCER (Husch 2017)' % PaCER
                    try
                        [coords_mm,trajectory,markers]=ea_runpacer(poptions);
                        options.native=1;
                    catch e % revert to TRAC/CORE
                        warning('PaCER or the LeadDBS Wrapper for PaCER failed with the following error:');
                        disp(['Identifier: ' e.identifier]);
                        disp(['Message: ' e.message]);
                        disp(['In: ' e.stack(1).file]);
                        disp(['Method: ' e.stack(1).name]);
                        disp(['Line: ' num2str(e.stack(1).line)]);
                        disp('Please check your input data carefully.');
                        disp('If the error persists, please consider a bug report at <a href="https://github.com/adhusch/PaCER/issues">https://github.com/adhusch/PaCER/issues</a>.');
                        ea_error('PaCER failed. Potentially try running the TRAC/CORE Algorithm or a fully manual pre-reconstruction.');
                    end

                case 'Manual' % Manual
                    [coords_mm,trajectory,markers]=ea_runmanual(poptions);
                    options.native=1;

                case 'Slicer (Manual)' % Manually mark lead head/tail in Slicer 3D
                    [coords_mm,trajectory,markers]=ea_runmanualslicer(poptions);
                    options.native=1;
            end
            options.hybridsave=1;
            options.elside=options.sides(1);
            elmodel=options.elmodel;
            ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
            if isfield(options,'hybridsave')
                options=rmfield(options,'hybridsave');
            end
        end
        options.native=wasnative; % restore original setting.
    end

    if options.manualheightcorrection
        poptions=ea_checkmanapproved(options);
        if ~isempty(poptions.sides)
            mcfig=figure('name',[options.patientname,': Electrode Reconstruction'],'numbertitle','off');
            %warning('off');
            try
                ea_maximize(mcfig);
            end
            options.elside=options.sides(1);
            ea_manualreconstruction(mcfig,options.patientname,options);
        end
    else
        ea_write(options)
    end
else
    ea_write(options)
end


function poptions=ea_checkmanapproved(options)
poptions=options;
if ~options.overwriteapproved && exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file') % only re-open not manually approved reconstructions.
    load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
    todel=[]; cnt=1;
    for side=options.sides
        try % index may exceed entries in legacy saves
            if reco.props(side).manually_corrected
                todel(cnt)=side; cnt=cnt+1;
            end
        end
    end
    [~,ix]=ismember(todel,poptions.sides);
    poptions.sides(ix)=[]; % do not re-reconstruct the ones already approved.

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

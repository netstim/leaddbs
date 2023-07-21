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
if isfield(options, 'leadfigure')
    bids = getappdata(options.leadfigure, 'bids');
    subjId = getappdata(options.leadfigure, 'subjId');
    if ~isempty(bids)
        options.bids = bids;
        if ~isempty(subjId{options.subjInd})
            options.subj = bids.getSubj(subjId{options.subjInd}, options.modality);
        end
    end
end

if strcmp(options.leadprod, 'dbs')
    if options.importdcm.do && bids.subjDataOverview.hasSourcedata(options.subj.subjId)
        niftis = ea_dcm_to_nii(options.subj.sourcedataDir, fullfile(options.subj.rawdataDir, 'unsorted'), options.importdcm.tool);
        if isempty(niftis)
            ea_cprintf('CmdWinWarnings', 'No output found from DICOM to NIfTI conversion for "%s"!\n', options.subj.subjId);
            options.leadfigure.focus;
            return;
        end
    end

    if options.importdcm.do || options.importnii.do
        unsortedFiles = ea_regexpdir(fullfile(options.subj.rawdataDir, 'unsorted'), '.*\.nii(\.gz)?');
        if ~isempty(unsortedFiles)
            ea_nifti_to_bids(unsortedFiles, bids.datasetDir, ['sub-', options.subj.subjId]);
            ea_delete(fullfile(options.subj.rawdataDir, 'unsorted'));
            ea_genrawimagesjson(bids.datasetDir, options.subj.subjId);
        else
            ea_cprintf('CmdWinWarnings', 'No unsorted raw images found for "%s"!\n', options.subj.subjId);
        end

        options.leadfigure.focus;
        return;
    end
end

% get accurate electrode specifications and save it in options.
options = ea_resolve_elspec(options);

% check connectome-mapper tags
if isfield(options,'lcm')
    ea_lcm(options);
end

if isfield(options,'predict')
   ea_predict(options);
end

% only 3D-rendering viewer can be opened if no patient is selected.
if ~strcmp(options.patientname,'No Patient Selected') && ~isempty(options.patientname)
    if isfile(fullfile(options.bids.datasetDir, 'miniset.json'))
        isMiniset = 1;
    else
        isMiniset = 0;
    end

    % Copy post-op images to preprocessing folder, no preproc is done for now
    if  isMiniset || ~isfield(options.subj, 'postopAnat')
        fields = {};
    else
        fields = fieldnames(options.subj.postopAnat);
    end

    for i=1:length(fields)
        if ~isfile(options.subj.postopAnat.(fields{i}).preproc)
            ea_mkdir(fileparts(options.subj.postopAnat.(fields{i}).preproc));
            % Copy file to preproc, take care of .nii.gz raw image
            if strcmp(options.prefs.niiFileExt, '.nii')
                copyfile(options.subj.postopAnat.(fields{i}).raw, [options.subj.postopAnat.(fields{i}).preproc, '.gz']);
                gunzip([options.subj.postopAnat.(fields{i}).preproc, '.gz']);
                delete([options.subj.postopAnat.(fields{i}).preproc, '.gz']);
            else
                copyfile(options.subj.postopAnat.(fields{i}).raw, options.subj.postopAnat.(fields{i}).preproc);
            end
        end
    end

    % Preprocessing pre-op images
    preprocessing = 0;
    if isMiniset || ~isfield(options.subj, 'preopAnat')
        fields = {};
    else
        fields = fieldnames(options.subj.preopAnat);
    end

    for i=1:length(fields)
        if ~isfile(options.subj.preopAnat.(fields{i}).preproc)
            % Copy files
            ea_mkdir(fileparts(options.subj.preopAnat.(fields{i}).preproc));

            % Copy file to preproc, take care of .nii.gz raw image
            if strcmp(options.prefs.niiFileExt, '.nii')
                copyfile(options.subj.preopAnat.(fields{i}).raw, [options.subj.preopAnat.(fields{i}).preproc, '.gz']);
                gunzip([options.subj.preopAnat.(fields{i}).preproc, '.gz']);
                delete([options.subj.preopAnat.(fields{i}).preproc, '.gz']);
            else
                copyfile(options.subj.preopAnat.(fields{i}).raw, options.subj.preopAnat.(fields{i}).preproc);
            end

            % Run reorientation, cropping and bias field correction
            ea_anatpreprocess(options.subj.preopAnat.(fields{i}).preproc);

            % Preprocessing steps only for pre-op anchor image
            if i==1
                % Skip reslicing anchor image to 0.7^3 resolution when only pre-op images exist
                if isfield(options.subj, 'postopAnat')
                    ea_resliceanat(options.subj.preopAnat.(fields{i}).preproc);
                end
                % ea_acpcdetect(options.subj.preopAnat.(fields{i}).preproc);
            end

            preprocessing = 1;
        end
    end

    if preprocessing
        fprintf('\nPreprocessing finished.\n\n');
    end

    % Set primary template
    subjAnchor = regexprep(options.subj.AnchorModality, '[^\W_]+_', '');
    if ismember(subjAnchor, fieldnames(bids.spacedef.norm_mapping))
        options.primarytemplate = bids.spacedef.norm_mapping.(subjAnchor);
    else
        options.primarytemplate = bids.spacedef.misfit_template;
    end

    % Pre-coregister pre-op anchor image
    if ~isMiniset && isfield(options.subj, 'preopAnat')
        if ~isfile(options.subj.preopAnat.(fields{1}).coreg)
            ea_precoreg(options.subj.preopAnat.(fields{1}).preproc, ... % Input anchor image
                options.primarytemplate, ... % Template to use
                options.subj.preopAnat.(fields{1}).coreg, ... % Output pre-coregistered image
                options.subj.coreg.transform.(fields{1})); % % Pre-coregistration transform

            % Check if anchor image has been properly pre-coregistered. If
            % not, fallback to preproc image.
            try
                load_nii(options.subj.preopAnat.(fields{1}).coreg);
            catch
                ea_cprintf('CmdWinWarnings', 'Anchor image was not properly pre-coregistered. Fallback to preproc image instead.');
                ea_delete(options.subj.coreg.transform.(fields{1}));
                copyfile(options.subj.preopAnat.(fields{1}).preproc, options.subj.preopAnat.(fields{1}).coreg);
            end
        end
    end

    coregDone = 0;

    if options.coregmr.do
        % Coregister pre-op MRIs to pre-op anchor image
        % TODO: coreg_fa disabled currently
        coregDone = ea_coregpreopmr(options);
    end

    if strcmp(options.subj.postopModality, 'MRI') && options.coregmr.do
        % Coregister post-op MRI to pre-op MRI
        coregDone = ea_coregpostopmr(options) || coregDone;
    end

    if strcmp(options.subj.postopModality, 'CT') && options.coregct.do
        % Coregister post-op CT to pre-op MRI
        coregDone = ea_coregpostopct(options) || coregDone;
    end

    if coregDone
        % Generate checkreg figures for coregistration
        ea_gencheckregfigs(options, 'coreg');
    end

    if options.normalize.do
        
        normlock = ea_reglocked(options, options.subj.preopAnat.(options.subj.AnchorModality).norm);
        
        if contains(options.normalize.method, 'apply')
            doit = true;
        elseif normlock == 1 % Both pre-op images coreg and norm were approved.
            doit = false;
        elseif normlock == 0.5 % Pre-op images coreg changed, rerun when using multispectral norm method.
            [~, ~, ~, doit] = eval([options.normalize.method,'(''prompt'')']);
        else
            doit = true;
        end
        
        if doit
            ea_normalize(options);
            ea_gencheckregfigs(options, 'norm');
        end
    end

    if options.scrf.do
        if ~ea_reglocked(options, options.subj.brainshift.anat.scrf) || options.overwriteapproved
            options.autobrainshift = 1;
            ea_dumpmethod(options, 'brainshift');
            ea_subcorticalrefine(options);
            options = rmfield(options,'autobrainshift');
        end
    end

    if isfield(options,'gencheckreg') % Exception when calling from the Tools menu.
        if options.gencheckreg
            ea_gencheckregfigs(options); % generate checkreg figures
        end
    end

    if options.dolc % perform lead connectome subroutine..
        ea_perform_lc(options);
    end

    if options.d2.write || options.d3.write
        if options.atl.genpt % generate patient specific atlas set
            ea_ptspecific_atl(options);
        end
    end

    if options.checkreg
        % Export checkreg figures
        if isempty(ea_regexpdir([options.subj.coregDir, filesep, 'checkreg'], '^(?!\.).*\.png$'))
            ea_gencheckregfigs(options, 'coreg');
        end

        if isempty(ea_regexpdir([options.subj.normDir, filesep, 'checkreg'], '^(?!\.).*\.png$'))
            ea_gencheckregfigs(options, 'norm');
        end

        ea_checkreg(options);
        drawnow; % Prevents the figure from changing the name with multiple subjects

        e=evalin('base', 'checkregempty');
        evalin('base',' clear checkregempty');

        if e && ~ea_reglocked(options, options.subj.brainshift.anat.scrf) ...
             && isfile(options.subj.brainshift.transform.instore)
            ea_subcorticalrefine(options);
        end
    end

    if options.normalize.refine
        ea_runwarpdrive(options, '0');
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
               if exist([options.subj.freesurferDir,filesep,'sub-',options.subj.subjId,filesep],'dir')
                   if options.overwriteapproved
                       % for now still ask user to confirm recalculation
                       % since fs takes so long.
                       answ=questdlg('Existing FreeSurfer output folder found. Are you sure you want to recalculate results & overwrite?', ...
                          'FreeSurfer output found','Recalculate & Overwrite','Skip','Skip');

                       switch lower(answ)
                           case 'recalculate & overwrite'
                               ea_runfreesurfer(options)
                       end
                   end
               else
                    ea_runfreesurfer(options);
               end
       end
    end

    if options.doreconstruction
        wasnative = options.native;
        poptions = ea_checkmanapproved(options);
        ea_mkdir(options.subj.reconDir);
        if ~isempty(poptions.sides)
            if numel(options.uipatdirs) > 1
                % Override recon method in multiple patients case
                options.reconmethod = options.prefs.reco.method.(options.subj.postopModality);
            end
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

    if options.refinelocalization
        poptions=ea_checkmanapproved(options);
        if ~isempty(poptions.sides)
            mcfig=figure('name',[options.subj.subjId,': Electrode Reconstruction'],'numbertitle','off');
            mcfig.WindowState = 'maximized';
            options.elside=options.sides(1);
            ea_manualreconstruction(mcfig,options.subj.subjId,options);
        end
    else
        ea_write(options)
    end
else
    ea_write(options)
end


function options = ea_checkmanapproved(options)
if ~options.overwriteapproved && isfile(options.subj.recon.recon) % only re-open not manually approved reconstructions.
    load(options.subj.recon.recon, 'reco');
    to_delete = find([reco.props.manually_corrected]);
    [~,idx] = ismember(to_delete, options.sides);
    options.sides(idx) = []; % do not re-reconstruct the ones already approved.
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

function ea_calc_biophysical_lg(M, options, selection, selectedConn, selectedParc, handles)

% handles not necessary for most cases (just for patient specific fMRI).

for pt=selection
    % set pt specific options
    [options.root, options.patientname] = fileparts(M.patient.list{pt});
    options.root = [options.root, filesep];

    options = ea_getptopts(fullfile(options.root, options.patientname), options);

    fprintf('\nProcessing %s...\n\n', options.patientname);
    try
        options.numcontacts=size(M.elstruct(pt).coords_mm{1},1);
    catch % no localization present or in wrong format.
        ea_error(['Please localize ',options.patientname,' first.']);
    end
    options.d3.verbose='off';
    options.d3.elrendering=1;	% hard code to viz electrodes in this setting.
    options.d3.exportBB=0;	% don't export brainbrowser struct by default
    options.d3.colorpointcloud=0;

    options.d3.hlactivecontacts=0;
    options.d3.showactivecontacts=0;
    options.d3.showpassivecontacts=0;
    try
        options.d3.isomatrix=M.isomatrix;
    catch
        options.d3.isomatrix={};
    end
    try
        options.d3.isomatrix_name=M.isomatrix_name;
    catch
        options.d3.isomatrix_name={};
    end
    options.normregressor=M.ui.normregpopup;
    for reg=1:length(options.d3.isomatrix)
        try
            options.d3.isomatrix{reg}=ea_reformat_isomatrix(options.d3.isomatrix{reg},M,options);
        end
    end

    options.d3.isovscloud=M.ui.isovscloudpopup;
    options.d3.showisovolume=M.ui.showisovolumecheck;
    options.d3.exportBB=0;
    options.expstatvat.do=0;
    try
        options.expstatvat.vars=M.clinical.vars(M.ui.clinicallist);
        options.expstatvat.labels=M.clinical.labels(M.ui.clinicallist);
        options.expstatvat.pt=pt;
    end
    options.expstatvat.dir=M.root;

    % Step 1: Re-calculate closeness to subcortical atlases.
    options.leadprod = 'group';
    options.patient_list=M.patient.list;
    options.d3.mirrorsides=0;

    resultfig=ea_elvis(options, M.elstruct, pt);

    if ~isfield(options.subj, 'norm')
        ea_cprintf('CmdWinWarnings', 'Running in Miniset mode: %s...\n', options.subj.subjId);
        volumespresent=0;
    elseif isempty(dir([options.subj.norm.transform.inverseBaseName, '*']))
        ea_cprintf('CmdWinWarnings', 'Tranformation not found for %s...\n', options.subj.subjId);
        volumespresent=0;
    else
        volumespresent=1;
    end

    % Step 2: Re-calculate VAT
    if isfield(M,'S')
        try
            setappdata(resultfig,'curS',M.S(pt));
        catch
            ea_error(['Stimulation parameters for ', options.subj.subjId, ' are not set.']);
        end

        vfs = ea_regexpdir(ea_getearoot, 'ea_genvat_.*\.m$', 0);
        vfs = regexp(vfs, '(ea_genvat_.*)(?=\.m)', 'match', 'once');
        vfnames = cellfun(@(x) eval([x, '(''prompt'');']), vfs, 'Uni', 0);

        [~,ix]=ismember(M.vatmodel,vfnames);
        try
            ea_genvat=eval(['@',vfs{ix}]);
        catch
            keyboard
        end

        options=getappdata(resultfig,'options'); % selected atlas could have refreshed.

        options.orignative=options.native; % backup
        options.native=~ea_getprefs('vatsettings.estimateInTemplate'); % see whether VTAs should be directly estimated in template space or not
        if options.native && ~volumespresent
            ea_cprintf('CmdWinWarnings', 'Calculating VTA in template space since patient folder %s is incomplete.\n', options.subj.subjId);
            options.native=0;
        end

        %setappdata(handles.leadfigure,'resultfig',resultfig);
        setappdata(resultfig,'elstruct',M.elstruct(pt));
        setappdata(resultfig,'elspec',options.elspec);

        if options.native % Reload native space coordinates
            coords = ea_load_reconstruction(options);
        else
            coords = M.elstruct(pt).coords_mm;
        end

        vatCalcPassed = [0 0];
        stimparams = struct();
        if strcmp(M.S(pt).model, 'OSS-DBS (Butenko 2020)')
            if options.prefs.machine.vatsettings.butenko_calcAxonActivation
                feval(ea_genvat,M.S(pt),options);
                ea_cprintf('CmdWinWarnings', 'OSS-DBS axon activation mode detect, skipping calc stats for %s!\n', options.patientname);
                continue;
            else
                [vatCalcPassed, stimparams] = feval(ea_genvat,M.S(pt),options);
            end
        else
            for side=1:2
                try
                    [vtafv,vtavolume] = feval(ea_genvat,coords,M.S(pt),side,options,['gs_',M.guid]);
                    vatCalcPassed(side) = 1;
                catch
                    vtafv=[];
                    vtavolume=0;
                    vatCalcPassed(side) = 0;
                end
                stimparams(1,side).VAT(1).VAT = vtafv;
                stimparams(1,side).volume = vtavolume;
            end
        end

        options.native=options.orignative; % restore
        setappdata(resultfig,'stimparams',stimparams(1,:));
    end

    % Calc VAT stats (atlas intersection and volume)
    if all(vatCalcPassed)
        ea_calc_vatstats(resultfig,options);
    else
        ea_cprintf('CmdWinErrors', 'Failed to calculate VTA for patient %s side %s!\n', options.patientname, num2str(find(~vatCalcPassed)));
    end

    % Step 3: Re-calculate connectivity from VAT to rest of the brain.
    if all(vatCalcPassed)
        % Convis part:
        directory = [options.root, options.patientname, filesep];
        if ischar(selectedConn) && ~ismember(selectedConn, {'Do not calculate connectivity stats', ''})
            switch selectedConn 
                case 'Patient''s fMRI time courses'
                    ea_error('Group statistics for fMRI are not yet supported. Sorry, check back later!');
                    pV=spm_vol([ea_space(options,'labeling'),selectedParc,'.nii']);
                    pX=spm_read_vols(pV);
                    ea_cvshowvatfmri(resultfig,pX,directory,filesare,handles,pV,selectedParc,selectedConn,options);
                case 'Patient''s fiber tracts'
                    ea_cvshowvatdmri(resultfig,directory,{selectedConn,'gs'},selectedParc,options);
            end
        elseif isstruct(selectedConn)
            ea_cvshowvatdmri(resultfig,directory,{selectedConn,'gs'},selectedParc,options);
        end
    end

    close(resultfig);
end

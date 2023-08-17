classdef BIDSFetcher
    % BIDS dataset fetcher

    properties
        settings
        spacedef
        datasetDir
        subjId
        subjDataOverview
    end
    
    properties (Access = private, Constant)
        anchorSpace = 'anchorNative';
    end

    methods
        %% Constructor
        function obj = BIDSFetcher(datasetDir, verbose)
            if ~exist('verbose', 'var') || isempty(verbose)
                verbose = 0;
            end

            if ~isfolder(datasetDir)
                error('Specified dataset folder doesn''t exist!');
            end

            % Set up properties
            obj.settings = obj.leadPrefs('m');
            obj.spacedef = ea_getspacedef;
            obj.datasetDir = GetFullPath(datasetDir);

            % Check dataset description file
            if ~isfile(fullfile(obj.datasetDir, 'dataset_description.json'))
                ea_cprintf('CmdWinWarnings', 'Could not find dataset description file, generating one now...\n');
                ea_generate_datasetDescription(obj.datasetDir, 'root_folder');
            end

            obj.subjId = obj.getSubjId;
            obj.subjDataOverview = obj.getSubjDataOverview;

            % Check rawimages.json
            if isempty(obj.subjId)
                ea_cprintf('CmdWinWarnings', 'No subj found in the dataset: %s\n', obj.datasetDir);
            else
                subjWithoutRawimagesJson = obj.subjDataOverview.Row(~obj.subjDataOverview.hasRawimagesJson & obj.subjDataOverview.hasRawdata);
                if ~isempty(subjWithoutRawimagesJson)
                    for i=1:length(subjWithoutRawimagesJson)
                        ea_genrawimagesjson(obj.datasetDir, subjWithoutRawimagesJson{i});
                    end
                    obj.subjDataOverview = obj.getSubjDataOverview;
                end
            end

            % TODO: BIDS validation

            % Verbose
            if verbose
                fprintf('\nLoaded BIDS dataset at: %s.\nFound the following subjects:\n', obj.datasetDir);
                fprintf('%s\n', obj.subjId{:});
                fprintf('\n');
            end
        end

        %% Data fetching functions
        function subjId = getSubjId(obj)
            % Find subject folders: sub-*
            subjSourceDirs = ea_regexpdir([obj.datasetDir, filesep, 'sourcedata'], '^sub-[^\W_]+$', 0, 'dir');
            if ~isempty(subjSourceDirs)
                ea_mkdir(strrep(subjSourceDirs, 'sourcedata', 'rawdata'));
                ea_mkdir(strrep(subjSourceDirs, 'sourcedata', ['derivatives', filesep, 'leaddbs']));
            end

            subjRawDirs = ea_regexpdir([obj.datasetDir, filesep, 'rawdata'], '^sub-[^\W_]+$', 0, 'dir');
            if ~isempty(subjRawDirs)
                ea_mkdir(strrep(subjRawDirs, 'rawdata', ['derivatives', filesep, 'leaddbs']));
            end

            subjDirs = ea_regexpdir([obj.datasetDir, filesep, 'derivatives', filesep, 'leaddbs'], '^sub-[^\W_]+$', 0, 'dir');

            subjId = regexp(subjDirs, ['(?<=\', filesep, 'sub-)[^\W_]+'], 'match', 'once');
        end

        function subjDataOverview = getSubjDataOverview(obj, subjId)
            if ~exist('subjId', 'var')
                subjId = obj.subjId;
            end

            if ischar(subjId)
                subjId = {subjId};
            end

            subjDataOverview = false(numel(subjId), 4);

            for i=1:numel(subjId)
                subjDerivativesDir = fullfile(obj.datasetDir, 'derivatives', 'leaddbs', ['sub-', subjId{i}]);
                subjRawdataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId{i}]);
                subjSourcedataDir = fullfile(obj.datasetDir, 'sourcedata', ['sub-', subjId{i}]);

                if isfile(fullfile(subjDerivativesDir, 'prefs', ['sub-', subjId{i}, '_desc-rawimages.json']))
                    subjDataOverview(i, 1) = true;
                end

                bidsRawFiles = ea_regexpdir(subjRawdataDir, 'ses-(pre|post)op.*\.nii\.gz$');
                if ~isempty(bidsRawFiles)
                    subjDataOverview(i, 2) = true;
                end

                unsortedRawFiles = ea_regexpdir(fullfile(subjRawdataDir, 'unsorted'), '.*\.nii(\.gz)?$');
                if ~isempty(unsortedRawFiles)
                    subjDataOverview(i, 3) = true;
                end

                souceFiles = ea_regexpdir(subjSourcedataDir, '.*', 0, 'a');
                if ~isempty(souceFiles)
                    subjDataOverview(i, 4) = true;
                end
            end

            subjDataOverview = array2table(subjDataOverview, 'VariableNames', {'hasRawimagesJson', 'hasRawdata', 'hasUnsortedRawdata', 'hasSourcedata'}, 'RowNames', subjId);
        end

        function rawImages = getRawImages(obj, subjId)
            if isfile(obj.getPrefs(subjId, 'rawimages'))
                rawImages = loadjson(obj.getPrefs(subjId, 'rawimages'));
            else
                rawImages = ea_genrawimagesjson(obj.datasetDir, subjId);
            end

            if isfield(rawImages, 'preop')
                preopAnat = struct2cell(rawImages.preop.anat);
                if ~startsWith(preopAnat{1}, ['sub-', subjId, '_'])
                    wrongSubjId = regexp(preopAnat{1}, '(?<=^sub-)[^\W_]+(?=_)', 'match', 'once');
                    ea_cprintf('CmdWinWarnings', 'Mismatched subjId "%s" found in the rawimages.json! Should be "%s".\n', wrongSubjId, subjId);
                    rawImages = ea_genrawimagesjson(obj.datasetDir, subjId);
                end
            end
        end

        function LeadDBSDirs = getLeadDBSDirs(obj, subjId)
            subjDir = fullfile(obj.datasetDir, 'derivatives', 'leaddbs', ['sub-', subjId]);
            if ~isfolder(subjDir)
                error('Subject ID %s doesn''t exist!', subjId);
            end

            LeadDBSDirs.subjDir = subjDir;
            LeadDBSDirs.atlasDir = fullfile(subjDir, 'atlases');
            LeadDBSDirs.brainshiftDir = fullfile(subjDir, 'brainshift');
            LeadDBSDirs.clinicalDir = fullfile(subjDir, 'clinical');
            LeadDBSDirs.coregDir = fullfile(subjDir, 'coregistration');
            LeadDBSDirs.exportDir = fullfile(subjDir, 'export');
            LeadDBSDirs.logDir = fullfile(subjDir, 'log');
            LeadDBSDirs.normDir = fullfile(subjDir, 'normalization');
            LeadDBSDirs.prefsDir = fullfile(subjDir, 'prefs');
            LeadDBSDirs.preprocDir = fullfile(subjDir, 'preprocessing');
            LeadDBSDirs.reconDir = fullfile(subjDir, 'reconstruction');
            LeadDBSDirs.stimDir = fullfile(subjDir, 'stimulations');
            LeadDBSDirs.warpdriveDir = fullfile(subjDir, 'warpdrive');
            LeadDBSDirs.acpcDir = fullfile(subjDir, 'acpc');
            LeadDBSDirs.freesurferDir = fullfile(fileparts(fileparts(subjDir)),'freesurfer');
        end

        function prefs = getPrefs(obj, subjId, label, format)
            % Get files from prefs folder
            if ~exist('format', 'var') || isempty(format)
                format = '.json';
            end

            if ~startsWith(format, '.')
                format = ['.', format];
            end

            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            prefs = [LeadDBSDirs.prefsDir, filesep, 'sub-', subjId, '_desc-', label, format];
        end

        function [preferMRCT, bothMRCTPresent] = checkModality(obj, subjId, preferMRCT)
            if isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            else
                if ischar(preferMRCT)
                    switch upper(preferMRCT)
                        case {'MRI', 'MR'}
                            preferMRCT = 1;
                        case 'CT'
                            preferMRCT = 2;
                    end
                end
            end

            % Get raw images struct
            rawImages = obj.getRawImages(subjId);

            % Return if no post-op image available
            if ~isfield(rawImages, 'postop')
                preferMRCT = 3;
                bothMRCTPresent = 0;
                return;
            elseif preferMRCT == 3 % update uiprefs modality if postop image is found but preferMRCT is still set to 3 in uiprefs, in this case update uiprefs to defaults
                preferMRCT = obj.settings.preferMRCT;
            end

            % Get post-op image modalities
            modality = fieldnames(rawImages.postop.anat);

            % Check presence of MRI and CT
            MRPresent = ismember('ax_MRI', modality);
            CTPresent = ismember('CT', modality);

            if CTPresent && MRPresent
                bothMRCTPresent = 1;
            else
                bothMRCTPresent = 0;
            end

            % Adapt preferMRCT according to post-op image available
            if MRPresent && (preferMRCT == 1  || preferMRCT == 2 && ~CTPresent)
                preferMRCT = 1;
            elseif CTPresent && (preferMRCT == 2  || preferMRCT == 1 && ~MRPresent)
                preferMRCT = 2;
            else
                error('Post-op images not properly defined!')
            end
        end

        function subj = getSubj(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            [preferMRCT, bothMRCTPresent] = checkModality(obj, subjId, preferMRCT);

            % Set subj dirs
            subj = obj.getLeadDBSDirs(subjId);

            % Set misc fields
            subj.subjId = subjId;
            subj.uiprefs = obj.getPrefs(subjId, 'uiprefs', 'mat');
            subj.methodLog = obj.getLog(subjId, 'methods');
            subj.rawdataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId]);
            subj.sourcedataDir = fullfile(obj.datasetDir, 'sourcedata', ['sub-', subjId]);

            % Set pre-op anat field
            preopAnat = obj.getPreopAnat(subjId);
            preopFields = fieldnames(preopAnat);
            for i=1:length(preopFields)
                subj.preopAnat.(preopFields{i}).raw = [preopAnat.(preopFields{i}), '.gz'];
            end

            % Set pre-op anchor modality
            if ~isempty(preopFields)
                subj.AnchorModality = preopFields{1};
            else
                subj.AnchorModality = obj.settings.prenii_order{1};
            end

            % Set post-op anat field
            if preferMRCT ~= 3
                postopAnat = obj.getPostopAnat(subjId, preferMRCT);
                postopFields = fieldnames(postopAnat);
                for i=1:length(postopFields)
                    subj.postopAnat.(postopFields{i}).raw = [postopAnat.(postopFields{i}), '.gz'];
                end
            end

            % Set post-op modality
            if preferMRCT == 1
                subj.postopModality = 'MRI';
            elseif preferMRCT == 2
                subj.postopModality = 'CT';
            elseif preferMRCT == 3
                subj.postopModality = 'None';
            end

            % Set bothMRCTPresent flag
            subj.bothMRCTPresent = bothMRCTPresent;

            % Set pipeline fields
            if ~isempty(preopFields)
                subj.preproc.anat = obj.getPreprocAnat(subjId, preferMRCT);
                subj.coreg.anat = obj.getCoregAnat(subjId, preferMRCT);
                subj.coreg.transform = obj.getCoregTransform(subjId, preferMRCT);
                subj.coreg.log = obj.getCoregLog(subjId);
                subj.coreg.checkreg = obj.getCoregCheckreg(subjId, preferMRCT);
    
                if preferMRCT ~= 3
                    subj.brainshift.anat = obj.getBrainshiftAnat(subjId, preferMRCT);
                    subj.brainshift.transform = obj.getBrainshiftTransform(subjId);
                    subj.brainshift.log = obj.getBrainshiftLog(subjId);
                    subj.brainshift.checkreg = obj.getBrainshiftCheckreg(subjId, preferMRCT);
                end
    
                subj.norm.anat = obj.getNormAnat(subjId, preferMRCT);
                subj.norm.transform = obj.getNormTransform(subjId);
                subj.norm.log = obj.getNormLog(subjId);
                subj.norm.checkreg = obj.getNormCheckreg(subjId, preferMRCT);
    
                % Set pre-op preprocessed images
                for i=1:length(preopFields)
                    subj.preopAnat.(preopFields{i}).preproc = subj.preproc.anat.preop.(preopFields{i});
                end
    
                % Set pre-op coregistered images
                for i=1:length(preopFields)
                    subj.preopAnat.(preopFields{i}).coreg = subj.coreg.anat.preop.(preopFields{i});
                end
    
                % Set pre-op normalized images
                subj.preopAnat.(preopFields{1}).norm = subj.norm.anat.preop.(preopFields{1});
            end

            if preferMRCT ~= 3 && ~isempty(postopFields)
                if preferMRCT ~= 3
                    % Set post-op preprocessed images
                    for i=1:length(postopFields)
                        subj.postopAnat.(postopFields{i}).preproc = subj.preproc.anat.postop.(postopFields{i});
                    end
    
                    % Set post-op coregistered images
                    for i=1:length(postopFields)
                        subj.postopAnat.(postopFields{i}).coreg = subj.coreg.anat.postop.(postopFields{i});
                    end
        
                    % Set post-op coregistered tone-mapped CT
                    if ismember('CT', postopFields)
                        subj.postopAnat.CT.coregTonemap = subj.coreg.anat.postop.tonemapCT;
                    end
        
                    % Set post-op normalized images
                    for i=1:length(postopFields)
                        subj.postopAnat.(postopFields{i}).norm = subj.norm.anat.postop.(postopFields{i});
                    end
        
                    % Set post-op normalized tone-mapped CT
                    if ismember('CT', postopFields)
                        subj.postopAnat.CT.normTonemap = subj.norm.anat.postop.tonemapCT;
                    end
                end
            end

            % Set reconstruction
            subj.recon = obj.getRecon(subjId, preferMRCT);

            % Set stats
            subj.stats = obj.getStats(subjId);
            
            % Set acpc autodetect
            subj.acpc.acpcAutodetect = fullfile(subj.acpcDir, ['sub-' subj.subjId '_desc-acpcautodetect.mat']);
            subj.acpc.acpcManual     = fullfile(subj.acpcDir, ['sub-' subj.subjId '_desc-acpcmanual.fcsv']);
        end

        function preopAnat = getPreopAnat(obj, subjId)
            % Set dirs
            rawDataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId]);

            % Get raw images struct
            rawImages = obj.getRawImages(subjId);

            % Return in case not found
            if ~isfield(rawImages, 'preop')
                preopAnat = struct;
                return;
            end

            % Get images and modalities
            images = fullfile(rawDataDir, 'ses-preop', 'anat', append(struct2cell(rawImages.preop.anat), obj.settings.niiFileExt));
            modality = fieldnames(rawImages.preop.anat)';

            % Set pre-defined orders
            preniiOrder = obj.settings.prenii_order;
            templateOrder = fieldnames(obj.spacedef.norm_mapping)';
            preopImageOrder = [preniiOrder, setdiff(templateOrder, preniiOrder, 'stable')];
            
            % Set pre-op anat images according to pre-defined orders
            for i=1:length(preopImageOrder)
                % Find the index in the present images
                idx = find(contains(modality, preopImageOrder{i}));
                if ~isempty(idx)
                    [~, ind] = ea_sortalike(lower(regexp(modality(idx), '[^\W_]+(?=_[^\W_]+)', 'match', 'once')), {'iso', 'ax', 'cor', 'sag'});
                    idx = idx(ind);

                    for j=1:length(idx)
                        preopAnat.(modality{idx(j)}) = images{idx(j)};
                    end

                    images(idx) = [];
                    modality(idx) = [];
                end
            end

            % Set other pre-op anat images
            if ~isempty(modality)
                [modality, index] = sort(modality);
                images = images(index);
                for i=1:length(modality)
                    preopAnat.(modality{i}) = images{i};
                end
            end

        end

        function postopAnat = getPostopAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Set dir
            rawDataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId]);

            % Get raw images struct
            rawImages = obj.getRawImages(subjId);

            % Return in case not found
            if ~isfield(rawImages, 'postop')
                postopAnat = struct;
                return;
            end

            % Get images and modalities
            images = fullfile(rawDataDir, 'ses-postop', 'anat', append(struct2cell(rawImages.postop.anat), obj.settings.niiFileExt));
            modality = fieldnames(rawImages.postop.anat);

            % Check presence of CT and MR
            CTPresent = ismember('CT', modality);
            MRPresent = any(contains(modality, 'ax_'));

            if CTPresent && (preferMRCT == 2  || preferMRCT == 1 && ~MRPresent)
                % Check post-op CT
                idx = find(ismember(modality, 'CT'), 1);
                if ~isempty(idx)
                    postopAnat.CT = images{idx};
                end
            elseif MRPresent && (preferMRCT == 1  || preferMRCT == 2 && ~CTPresent)
                % Check post-op axial MRI
                idx = find(contains(modality, 'ax'), 1);
                if ~isempty(idx)
                    postopAnat.(modality{idx}) = images{idx};
                end

                % Check post-op coronal MRI
                idx = find(contains(modality, 'cor'), 1);
                if ~isempty(idx)
                    postopAnat.(modality{idx}) = images{idx};
                end

                % Check post-op sagital MRI
                idx = find(contains(modality, 'sag'), 1);
                if ~isempty(idx)
                    postopAnat.(modality{idx}) = images{idx};
                end
            else
                error('Post-op images not properly defined!')
            end
        end

        function preprocAnat = getPreprocAnat(obj, subjId, preferMRCT)
            % Get pre-op anat images
            preopAnat = obj.getPreopAnat(subjId);

            % Get preprocessing directory
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseDir = fullfile(LeadDBSDirs.preprocDir, 'anat');

            % Get preprocessed pre-op anat images
            baseName = ['sub-', subjId, '_ses-preop_desc-preproc_'];
            fields = fieldnames(preopAnat);
            for i=1:length(fields)
                modality = fields{i};
                parsed = parseBIDSFilePath(preopAnat.(modality));
                preprocAnat.preop.(modality) = strrep(preopAnat.(modality), [parsed.dir, filesep, 'sub-', subjId, '_ses-preop_'], [baseDir, filesep, baseName]);
            end

            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            if preferMRCT == 3
                % No post-op image available
                return;
            else
                postopAnat = obj.getPostopAnat(subjId, preferMRCT);
                % Get preprocessed post-op anat images
                baseName = ['sub-', subjId, '_ses-postop_desc-preproc_'];
                if isfield(postopAnat, 'CT')
                    parsed = parseBIDSFilePath(postopAnat.CT);
                    preprocAnat.postop.CT = fullfile(baseDir, [baseName, 'CT', parsed.ext]);
                else
                    fields = fieldnames(postopAnat);
                    for i=1:length(fields)
                        modality = fields{i};
                        parsed = parseBIDSFilePath(postopAnat.(modality));
                        preprocAnat.postop.(modality) = fullfile(baseDir, [baseName, 'acq-', parsed.acq, '_', parsed.suffix, parsed.ext]);
                    end
                end
            end
        end

        function coregAnat = getCoregAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get preprocessed anat images
            preprocAnat = obj.getPreprocAnat(subjId, preferMRCT);

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set coregistered anat images
            session = fieldnames(preprocAnat);
            for i=1:length(session)
                modality = fieldnames(preprocAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(preprocAnat.(session{i}).(modality{j}), LeadDBSDirs.preprocDir, LeadDBSDirs.coregDir);
                    coregAnat.(session{i}).(modality{j}) = strrep(anat , ['ses-', session{i}, '_'], ['ses-', session{i}, '_space-', obj.anchorSpace, '_']);
                end
            end

            % Set tone-mapped CT
            if isfield(coregAnat, 'postop') && isfield(coregAnat.postop, 'CT')
                coregAnat.postop.tonemapCT = strrep(coregAnat.postop.CT, obj.anchorSpace, [obj.anchorSpace, '_rec-tonemapped']);
            end
        end

        function coregTransform = getCoregTransform(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get LeadDBS dirs and base name
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.coregDir, 'transformations', ['sub-', subjId, '_']);

            % Get coregistered images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Set pre-coregistration transformation
            fields = fieldnames(coregAnat.preop);
            coregTransform.(fields{1}) = [baseName, 'desc-precoreg_', fields{1}, '.mat'];

            % Set pre-op MR transformation
            if numel(fields) >1
                for i = 2:numel(fields)
                    spaceTag = erase(fields{i}, '_');
                    coregTransform.(fields{i}).forwardBaseName = [baseName, 'from-', spaceTag, '_to-', obj.anchorSpace, '_desc-'];
                    coregTransform.(fields{i}).inverseBaseName = [baseName, 'from-', obj.anchorSpace, '_to-', spaceTag, '_desc-'];
                end
            end

            % Set post-op CT transformation
            if isfield(coregAnat, 'postop')
                fields = fieldnames(coregAnat.postop);
                for i = 1:numel(fields)
                    if strcmp(fields{i}, 'CT')
                        coregTransform.CT.forwardBaseName = [baseName, 'from-CT_to-', obj.anchorSpace, '_desc-'];
                        coregTransform.CT.inverseBaseName = [baseName, 'from-', obj.anchorSpace, '_to-CT_desc-'];
                    else % Set post-op MRI transformation
                        spaceTag = erase(fields{i}, '_');
                        coregTransform.(fields{i}).forwardBaseName = [baseName, 'from-', spaceTag, '_to-', obj.anchorSpace, '_desc-'];
                        coregTransform.(fields{i}).inverseBaseName = [baseName, 'from-', obj.anchorSpace, '_to-', spaceTag, '_desc-'];
                    end
                end
            end
        end

        function coregLog = getCoregLog(obj, subjId)
            % Get LeadDBS dirs and base name
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.coregDir, 'log', ['sub-', subjId, '_desc-']);

            % Set coregistration log
            coregLog.method = [baseName, 'coregmethod.json'];
            coregLog.logBaseName = [baseName, 'coreg'];
        end

        function coregCheckreg = getCoregCheckreg(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Remove pre-op anchor anat image
            fields = fieldnames(coregAnat.preop);
            coregAnat.preop = rmfield(coregAnat.preop, fields(1));

            % Remove post-op CT image, will use tone-mapped CT
            if isfield(coregAnat, 'postop') && isfield(coregAnat.postop, 'CT')
                coregAnat.postop = rmfield(coregAnat.postop, 'CT');
            end

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            checkregDir = fullfile(LeadDBSDirs.coregDir, 'checkreg');

            % Set coregistered anat images
            coregCheckreg = struct;
            session = fieldnames(coregAnat);
            for i=1:length(session)
                modality = fieldnames(coregAnat.(session{i}));
                for j=1:length(modality)
                    coregCheckreg.(session{i}).(modality{j}) = setBIDSEntity(coregAnat.(session{i}).(modality{j}), 'dir', checkregDir, 'ext', '.png');
                end
            end
        end

        function brainshiftAnat = getBrainshiftAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set moving post-op image used for brain shift correction
            modality = fieldnames(coregAnat.postop);
            if preferMRCT == 2 && strcmp(obj.settings.scrf.tonemap, 'tp_')
                % Use tone-mapped CT
                brainshiftAnat.moving = coregAnat.postop.(modality{2});
            else
                % Use normal CT or MRI
                brainshiftAnat.moving = coregAnat.postop.(modality{1});
            end
            brainshiftAnat.moving = strrep(brainshiftAnat.moving, LeadDBSDirs.coregDir, LeadDBSDirs.brainshiftDir);

            % Use MRI suffix and ignore real modality
            if preferMRCT == 1 % Post-op MRI detected
                brainshiftAnat.moving = setBIDSEntity(brainshiftAnat.moving, 'acq', '', 'suffix', 'MRI');
            end

            % Set anchor anat image used for brain shift correction
            modality = fieldnames(coregAnat.preop);
            brainshiftAnat.anchor = strrep(coregAnat.preop.(modality{1}), obj.anchorSpace, [obj.anchorSpace, '_rec-brainshift']);
            brainshiftAnat.anchor = strrep(brainshiftAnat.anchor, LeadDBSDirs.coregDir, LeadDBSDirs.brainshiftDir);

            % Set brain shift corrected image
            if preferMRCT == 2 && strcmp(obj.settings.scrf.tonemap, 'tp_')
                % Use tone-mapped CT
                brainshiftAnat.scrf = setBIDSEntity(brainshiftAnat.moving, 'rec', 'tonemappedbrainshift');
            else
                % Use normal CT or MRI
                brainshiftAnat.scrf = strrep(brainshiftAnat.moving, obj.anchorSpace, [obj.anchorSpace, '_rec-brainshift']);
            end

            % Set masks used for brain shift correction
            baseDir = fullfile(LeadDBSDirs.brainshiftDir, 'anat');
            brainshiftAnat.secondstepmask = [baseDir, filesep, 'sub-', subjId, '_space-', obj.anchorSpace, '_desc-secondstepmask', obj.settings.niiFileExt];
            brainshiftAnat.thirdstepmask = [baseDir, filesep, 'sub-', subjId, '_space-', obj.anchorSpace, '_desc-thirdstepmask', obj.settings.niiFileExt];
        end

        function brainshiftTransform = getBrainshiftTransform(obj, subjId)
            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set base dir and base name
            baseDir = fullfile(LeadDBSDirs.brainshiftDir, 'transformations');
            baseName = ['sub-', subjId, '_from-', obj.anchorSpace, '_to-', obj.anchorSpace, 'BSC_desc-'];

            % Set brain shift transformations
            brainshiftTransform.instore = [baseDir, filesep, baseName, 'instore.mat'];
            brainshiftTransform.converted = [baseDir, filesep, baseName, 'converted.mat'];
            brainshiftTransform.scrf = [baseDir, filesep, baseName, 'scrf.mat'];
        end

        function brainshiftLog = getBrainshiftLog(obj, subjId)
            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set base dir
            baseDir = fullfile(LeadDBSDirs.brainshiftDir, 'log');

            % Set brain shift log
            brainshiftLog.method = [baseDir, filesep, 'sub-', subjId, '_desc-brainshiftmethod.json'];
        end

        function brainshiftCheckreg = getBrainshiftCheckreg(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get brain shift anat images
            brainshiftAnat = obj.getBrainshiftAnat(subjId, preferMRCT);

            % Set dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            anatDir = fullfile(LeadDBSDirs.brainshiftDir, 'anat');
            checkregDir = fullfile(LeadDBSDirs.brainshiftDir, 'checkreg');

            % Before brain shift correction
            brainshiftCheckreg.standard = strrep(brainshiftAnat.moving, anatDir, checkregDir);
            brainshiftCheckreg.standard = setBIDSEntity(brainshiftCheckreg.standard, 'ext', '.png');

            % After brain shift correction
            brainshiftCheckreg.scrf = strrep(brainshiftAnat.scrf, anatDir, checkregDir);
            brainshiftCheckreg.scrf = setBIDSEntity(brainshiftCheckreg.scrf,'ext', '.png');
        end

        function normAnat = getNormAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Remove pre-op anat images except for anchor image
            fields = fieldnames(coregAnat.preop);
            coregAnat.preop = rmfield(coregAnat.preop, fields(2:end));

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set normalized anat images
            templateSpace = obj.spacedef.name;
            session = fieldnames(coregAnat);
            for i=1:length(session)
                modality = fieldnames(coregAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(coregAnat.(session{i}).(modality{j}), LeadDBSDirs.coregDir, LeadDBSDirs.normDir);
                    normAnat.(session{i}).(modality{j}) = strrep(anat, obj.anchorSpace, templateSpace);
                end
            end
        end

        function normTransform = getNormTransform(obj, subjId)
            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set base dir and base name
            templateSpace = obj.spacedef.name;
            baseName = fullfile(LeadDBSDirs.normDir, 'transformations', ['sub-', subjId, '_from-']);

            % Set normalization transformations
            normTransform.forwardBaseName = [baseName, obj.anchorSpace, '_to-', templateSpace, '_desc-'];
            normTransform.inverseBaseName = [baseName, templateSpace, '_to-', obj.anchorSpace, '_desc-'];
        end

        function normLog = getNormLog(obj, subjId)
            % Get LeadDBS dirs and base name
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.normDir, 'log', ['sub-', subjId, '_desc-']);

            % Set normalization method and log
            normLog.method = [baseName, 'normmethod.json'];
            normLog.logBaseName = [baseName, 'norm'];
        end

        function normCheckreg = getNormCheckreg(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get normalized anat images
            normAnat = obj.getNormAnat(subjId, preferMRCT);

            % Remove post-op CT image, will use tone-mapped CT
            if isfield(normAnat, 'postop') && isfield(normAnat.postop, 'CT')
                normAnat.postop = rmfield(normAnat.postop, 'CT');
            end

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            checkregDir = fullfile(LeadDBSDirs.normDir, 'checkreg');

            % Set coregistered anat images
            session = fieldnames(normAnat);
            for i=1:length(session)
                modality = fieldnames(normAnat.(session{i}));
                for j=1:length(modality)
                    normCheckreg.(session{i}).(modality{j}) = setBIDSEntity(normAnat.(session{i}).(modality{j}), 'dir', checkregDir, 'ext', '.png');
                end
            end
        end

        function recon = getRecon(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var')
                preferMRCT = obj.settings.preferMRCT;
            end
            preferMRCT = checkModality(obj, subjId, preferMRCT);

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.reconDir, ['sub-', subjId, '_']);

            % Get reconstruction
            recon.recon = [baseName, 'desc-reconstruction.mat'];

            % Mask for CT-based reconstruction
            postopAnat = obj.getPostopAnat(subjId, preferMRCT);
            if isfield(postopAnat, 'CT')
                rawCTSpace = 'rawCT';
                recon.rawCTMask = [baseName, 'space-', rawCTSpace, '_desc-brainmask', obj.settings.niiFileExt];
                recon.anchorNativeMask = [baseName, 'space-', obj.anchorSpace, '_desc-brainmask', obj.settings.niiFileExt];
            end
        end

        function stats = getStats(obj, subjId)
            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Get stats
            stats = fullfile(LeadDBSDirs.subjDir, ['sub-', subjId, '_desc-stats.mat']);
        end

        function log = getLog(obj, subjId, label)
            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.logDir, ['sub-', subjId, '_desc-']);

            % Get log
            log = [baseName, label, '.txt'];
        end
    end

    methods(Static)
        %% Helper functions
        function prefs = leadPrefs(type)
            if ~exist('type', 'var') || isempty(type)
                type = 'm';
            end

            switch type
                case 'json'
                    % Read .ea_prefs.json
                    prefs = loadjson(fullfile(ea_gethome, '.ea_prefs.json'));
                case 'm'
                    % Read .ea_prefs.m and .ea_prefs.mat
                    prefs = ea_prefs;
            end
        end
    end
end

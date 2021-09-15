classdef BIDSFetcher
    % BIDS dataset fetcher

    properties
        settings
        spacedef
        datasetDir
        subjFolderNames
        subjId
    end

    methods
        %% Constructor
        function obj = BIDSFetcher(datasetDir, verbose)
            if ~exist('verbose', 'var') || isempty(verbose)
                verbose = 0;
            end

            % Set up properties
            obj.settings = obj.leadPrefs('m');
            obj.spacedef = ea_getspacedef;
            obj.datasetDir = GetFullPath(datasetDir);
            obj.subjFolderNames = obj.readSubjects;
            obj.subjId = strrep(obj.subjFolderNames, 'sub-', '');

            % TODO: BIDS validation

            % Verbose
            if verbose
                fprintf('\nLoaded BIDS dataset at: %s.\nFound the following subjects:\n', obj.datasetDir);
                fprintf('%s\n', obj.subjId{:});
                fprintf('\n');
            end
        end

        %% Data fetching functions
        function subjFolderNames = readSubjects(obj)
            % Find subject folders: sub-*
            subjDirs = ea_regexpdir([obj.datasetDir, filesep, 'rawdata'], 'sub-.*', 0);
            subjFolderNames = regexp(subjDirs, ['sub-.*(?=\', filesep, '$)'], 'match', 'once');
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
            LeadDBSDirs.stimDir = fullfile(subjDir, 'stimulation');
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

        function subj = getSubj(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Set subj dirs
            subj = obj.getLeadDBSDirs(subjId);

            % Set misc fields
            subj.subjId = subjId;
            subj.uiprefs = obj.getPrefs(subjId, 'uiprefs', 'mat');
            subj.methodLog = obj.getLog(subjId, 'methods');

            if ~isfile(obj.getPrefs(subjId, 'rawimages'))
                subj.rawImageJSONExist = 0;
                return;
            else
                subj.rawImageJSONExist = 1;
            end

            % Set pre-op anat field
            preopAnat = obj.getPreopAnat(subjId);
            preopFields = fieldnames(preopAnat);
            for i=1:length(preopFields)
                subj.preopAnat.(preopFields{i}).raw = preopAnat.(preopFields{i});
            end

            % Set pre-op anchor modality
            subj.AnchorModality = preopFields{1};

            % Set post-op anat field
            [postopAnat, bothMRCTPresent] = obj.getPostopAnat(subjId, preferMRCT);
            postopFields = fieldnames(postopAnat);
            for i=1:length(postopFields)
                subj.postopAnat.(postopFields{i}).raw = postopAnat.(postopFields{i});
            end

            % Set post-op modality
            if ismember('CT', postopFields)
                subj.postopModality = 'CT';
            else
                subj.postopModality = 'MRI';
            end

            % Set bothMRCTPresent flag
            subj.bothMRCTPresent = bothMRCTPresent;

            % Set pipeline fields
            subj.preproc.anat = obj.getPreprocAnat(subjId, preferMRCT);
            subj.coreg.anat = obj.getCoregAnat(subjId, preferMRCT);
            subj.coreg.transform = obj.getCoregTransform(subjId, preferMRCT);
            subj.coreg.log = obj.getCoregLog(subjId);
            subj.coreg.checkreg = obj.getCoregCheckreg(subjId, preferMRCT);
            subj.brainshift.anat = obj.getBrainshiftAnat(subjId, preferMRCT);
            subj.brainshift.transform = obj.getBrainshiftTransform(subjId);
            subj.brainshift.log = obj.getBrainshiftLog(subjId);
            subj.brainshift.checkreg = obj.getBrainshiftCheckreg(subjId, preferMRCT);
            subj.norm.anat = obj.getNormAnat(subjId, preferMRCT);
            subj.norm.transform = obj.getNormTransform(subjId);
            subj.norm.log = obj.getNormLog(subjId);
            subj.norm.checkreg = obj.getNormCheckreg(subjId, preferMRCT);

            % Set pre-op preprocessed images
            for i=1:length(preopFields)
                subj.preopAnat.(preopFields{i}).preproc = subj.preproc.anat.preop.(preopFields{i});
            end

            % Set post-op preprocessed images
            for i=1:length(postopFields)
                subj.postopAnat.(postopFields{i}).preproc = subj.preproc.anat.postop.(postopFields{i});
            end

            % Set pre-op coregistered images
            for i=1:length(preopFields)
                subj.preopAnat.(preopFields{i}).coreg = subj.coreg.anat.preop.(preopFields{i});
            end

            % Set post-op coregistered images
            for i=1:length(postopFields)
                subj.postopAnat.(postopFields{i}).coreg = subj.coreg.anat.postop.(postopFields{i});
            end

            % Set post-op coregistered tone-mapped CT
            if ismember('CT', postopFields)
                subj.postopAnat.CT.coregTonemap = subj.coreg.anat.postop.tonemapCT;
            end

            % Set pre-op normalized images
            subj.preopAnat.(preopFields{1}).norm = subj.norm.anat.preop.(preopFields{1});

            % Set post-op normalized images
            for i=1:length(postopFields)
                subj.postopAnat.(postopFields{i}).norm = subj.norm.anat.postop.(postopFields{i});
            end

            % Set post-op normalized tone-mapped CT
            if ismember('CT', postopFields)
                subj.postopAnat.CT.normTonemap = subj.norm.anat.postop.tonemapCT;
            end

            % Set reconstruction
            subj.recon = obj.getRecon(subjId, preferMRCT);
        end

        function preopAnat = getPreopAnat(obj, subjId)
            % Set dirs
            rawDataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId]);

            % Get raw images struct
            rawImages = loadjson(obj.getPrefs(subjId, 'rawimages'));

            % Get images and modalities
            images = fullfile(rawDataDir, 'ses-preop', 'anat', struct2cell(rawImages.preop.anat));
            modality = fieldnames(rawImages.preop.anat)';

            % Set pre-defined orders
            preniiOrder = obj.settings.prenii_order;
            templateOrder = fieldnames(obj.spacedef.norm_mapping)';
            preopImageOrder = [preniiOrder, setdiff(templateOrder, preniiOrder, 'stable')];

            % Set pre-op anat images according to pre-defined orders
            for i=1:length(preopImageOrder)
                % Find the index in the present images
                idx = find(ismember(modality, preopImageOrder{i}), 1);
                if ~isempty(idx)
                    preopAnat.(preopImageOrder{i}) = images{idx};
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

        function [postopAnat, bothMRCTPresent] = getPostopAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
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
                        case 'MR'
                            preferMRCT = 1;
                        case 'CT'
                            preferMRCT = 2;
                    end
                end
            end

            % Set dirs
            rawDataDir = fullfile(obj.datasetDir, 'rawdata', ['sub-', subjId]);
            subjDir = fullfile(obj.datasetDir, 'derivatives', 'leaddbs', ['sub-', subjId]);

            % Get raw images struct
            rawImages = loadjson(fullfile(subjDir, 'prefs', ['sub-', subjId, '_desc-rawimages.json']));

            % Get images and modalities
            images = fullfile(rawDataDir, 'ses-postop', 'anat', struct2cell(rawImages.postop.anat));
            modality = fieldnames(rawImages.postop.anat);

            % Check presence of CT and MR
            CTPresent = ismember('CT', modality);
            MRPresent = any(contains(modality, 'ax_'));
            if CTPresent && MRPresent
                bothMRCTPresent = 1;
            else
                bothMRCTPresent = 0;
            end

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
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get pre-op and post-op anat images
            preopAnat = obj.getPreopAnat(subjId);
            postopAnat = obj.getPostopAnat(subjId, preferMRCT);

            % Get preprocessing directory
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseDir = fullfile(LeadDBSDirs.preprocDir, 'anat');

            % Get preprocessed pre-op anat images
            baseName = ['sub-', subjId, '_desc-preproc_ses-preop_'];
            fields = fieldnames(preopAnat);
            for i=1:length(fields)
                modality = fields{i};
                parsed = parseBIDSFilePath(preopAnat.(modality));
                preprocAnat.preop.(modality) = fullfile(baseDir, [baseName, parsed.suffix, parsed.ext]);
            end

            % Get preprocessed post-op anat images
            baseName = ['sub-', subjId, '_desc-preproc_ses-postop_'];
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

        function coregAnat = getCoregAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get preprocessed anat images
            preprocAnat = obj.getPreprocAnat(subjId, preferMRCT);

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set coregistered anat images
            anchorSpace = 'anchorNative';
            session = fieldnames(preprocAnat);
            for i=1:length(session)
                modality = fieldnames(preprocAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(preprocAnat.(session{i}).(modality{j}), LeadDBSDirs.preprocDir, LeadDBSDirs.coregDir);
                    coregAnat.(session{i}).(modality{j}) = strrep(anat , [subjId, '_'], [subjId, '_space-', anchorSpace, '_']);
                end
            end

            % Set tone-mapped CT
            if isfield(coregAnat.postop, 'CT')
                coregAnat.postop.tonemapCT = strrep(coregAnat.postop.CT, anchorSpace, [anchorSpace, '_rec-tonemapped']);
            end
        end

        function coregTransform = getCoregTransform(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get LeadDBS dirs and base name
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.coregDir, 'transformations', ['sub-', subjId, '_']);

            % Get coregistered images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Set pre-coregistration transformation
            fields = fieldnames(coregAnat.preop);
            coregTransform.(fields{1}) = [baseName, 'desc-precoreg_', fields{1}, '.mat'];

            % Set post-op CT transformation
            if isfield(coregAnat.postop, 'CT')
                anchorSpace = 'anchorNative';
                coregTransform.CT.forwardBaseName = [baseName, 'from-CT_to-', anchorSpace, '_desc-'];
                coregTransform.CT.inverseBaseName = [baseName, 'from-', anchorSpace, '_to-CT_desc-'];
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
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Remove pre-op anchor anat image
            fields = fieldnames(coregAnat.preop);
            coregAnat.preop = rmfield(coregAnat.preop, fields(1));

            % Remove post-op CT image, will use tone-mapped CT
            if isfield(coregAnat.postop, 'CT')
                coregAnat.postop = rmfield(coregAnat.postop, 'CT');
            end

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            anatDir = fullfile(LeadDBSDirs.coregDir, 'anat');
            checkregDir = fullfile(LeadDBSDirs.coregDir, 'checkreg');

            % Set coregistered anat images
            session = fieldnames(coregAnat);
            for i=1:length(session)
                modality = fieldnames(coregAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(coregAnat.(session{i}).(modality{j}), anatDir, checkregDir);
                    parsed = parseBIDSFilePath(anat);
                    coregCheckreg.(session{i}).(modality{j}) = strrep(anat , parsed.ext, '.png');
                end
            end
        end

        function brainshiftAnat = getBrainshiftAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set anchor anat image used for brain shift correction
            anchorSpace = 'anchorNative';
            modality = fieldnames(coregAnat.preop);
            brainshiftAnat.anchor = strrep(coregAnat.preop.(modality{1}), anchorSpace, [anchorSpace, '_rec-brainshift']);
            brainshiftAnat.anchor = strrep(brainshiftAnat.anchor, LeadDBSDirs.coregDir, LeadDBSDirs.brainshiftDir);

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
                parsed = parseBIDSFilePath(brainshiftAnat.moving);
                brainshiftAnat.moving = strrep(brainshiftAnat.moving, ['_acq-', parsed.acq], '');
                brainshiftAnat.moving = strrep(brainshiftAnat.moving, parsed.suffix, 'MRI');
            end

            % Set masks used for brain shift correction
            baseDir = fullfile(LeadDBSDirs.brainshiftDir, 'anat');
            brainshiftAnat.secondstepmask = [baseDir, filesep, 'sub-', subjId, '_space-', anchorSpace, '_desc-secondstepmask', obj.settings.niiFileExt];
            brainshiftAnat.thirdstepmask = [baseDir, filesep, 'sub-', subjId, '_space-', anchorSpace, '_desc-thirdstepmask', obj.settings.niiFileExt];

            % Set brain shift corrected image
            brainshiftAnat.scrf = strrep(brainshiftAnat.moving, anchorSpace, [anchorSpace, '_rec-brainshift']);
        end

        function brainshiftTransform = getBrainshiftTransform(obj, subjId)
            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set base dir and base name
            anchorSpace = 'anchorNative';
            baseDir = fullfile(LeadDBSDirs.brainshiftDir, 'transformations');
            baseName = ['sub-', subjId, '_from-', anchorSpace, '_to-', anchorSpace, 'BSC_desc-'];

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
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get brain shift anat images
            brainshiftAnat = obj.getBrainshiftAnat(subjId, preferMRCT);

            % Set dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            anatDir = fullfile(LeadDBSDirs.brainshiftDir, 'anat');
            checkregDir = fullfile(LeadDBSDirs.brainshiftDir, 'checkreg');

            % Before brain shift correction
            brainshiftCheckreg.standard = strrep(brainshiftAnat.moving, anatDir, checkregDir);
            parsed = parseBIDSFilePath(brainshiftCheckreg.standard);
            brainshiftCheckreg.standard = strrep(brainshiftCheckreg.standard, parsed.ext, '.png');

            % After brain shift correction
            brainshiftCheckreg.scrf = strrep(brainshiftAnat.scrf, anatDir, checkregDir);
            parsed = parseBIDSFilePath(brainshiftCheckreg.scrf);
            brainshiftCheckreg.scrf = strrep(brainshiftCheckreg.scrf, parsed.ext, '.png');
        end

        function normAnat = getNormAnat(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get coregistered anat images
            coregAnat = obj.getCoregAnat(subjId, preferMRCT);

            % Remove pre-op anat images except for anchor image
            fields = fieldnames(coregAnat.preop);
            coregAnat.preop = rmfield(coregAnat.preop, fields(2:end));

            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set normalized anat images
            anchorSpace = 'anchorNative';
            templateSpace = obj.spacedef.name;
            session = fieldnames(coregAnat);
            for i=1:length(session)
                modality = fieldnames(coregAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(coregAnat.(session{i}).(modality{j}), LeadDBSDirs.coregDir, LeadDBSDirs.normDir);
                    normAnat.(session{i}).(modality{j}) = strrep(anat, anchorSpace, templateSpace);
                end
            end
        end

        function normTransform = getNormTransform(obj, subjId)
            % Get LeadDBS dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);

            % Set base dir and base name
            anchorSpace = 'anchorNative';
            templateSpace = obj.spacedef.name;
            baseName = fullfile(LeadDBSDirs.normDir, 'transformations', ['sub-', subjId, '_from-']);

            % Set normalization transformations
            normTransform.forwardBaseName = [baseName, anchorSpace, '_to-', templateSpace, '_desc-'];
            normTransform.inverseBaseName = [baseName, templateSpace, '_to-', anchorSpace, '_desc-'];
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
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get normalized anat images
            normAnat = obj.getNormAnat(subjId, preferMRCT);

            % Remove post-op CT image, will use tone-mapped CT
            if isfield(normAnat.postop, 'CT')
                normAnat.postop = rmfield(normAnat.postop, 'CT');
            end

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            anatDir = fullfile(LeadDBSDirs.normDir, 'anat');
            checkregDir = fullfile(LeadDBSDirs.normDir, 'checkreg');

            % Set coregistered anat images
            session = fieldnames(normAnat);
            for i=1:length(session)
                modality = fieldnames(normAnat.(session{i}));
                for j=1:length(modality)
                    anat = strrep(normAnat.(session{i}).(modality{j}), anatDir, checkregDir);
                    parsed = parseBIDSFilePath(anat);
                    normCheckreg.(session{i}).(modality{j}) = strrep(anat , parsed.ext, '.png');
                end
            end
        end

        function recon = getRecon(obj, subjId, preferMRCT)
            if ~exist('preferMRCT', 'var') || isempty(preferMRCT)
                uiprefsFile = getPrefs(obj, subjId, 'uiprefs', 'mat');
                % Check if uiprefs exists
                if isfile(uiprefsFile) % Set from uiprefs if it exists
                    uiprefs = load(uiprefsFile);
                    preferMRCT = uiprefs.modality;
                else % Set from LeadDBS prefs
                    preferMRCT = obj.settings.preferMRCT;
                end
            end

            % Get dirs
            LeadDBSDirs = obj.getLeadDBSDirs(subjId);
            baseName = fullfile(LeadDBSDirs.reconDir, ['sub-', subjId, '_']);

            % Get reconstruction
            recon.recon = [baseName, 'desc-reconstruction.mat'];

            % Mask for CT-based reconstruction
            postopAnat = obj.getPostopAnat(subjId, preferMRCT);
            if isfield(postopAnat, 'CT')
                rawCTSpace = 'rawCT';
                anchorSpace = 'anchorNative';
                recon.rawCTMask = [baseName, 'space-', rawCTSpace, '_desc-brainmask', obj.settings.niiFileExt];
                recon.anchorNativeMask = [baseName, 'space-', anchorSpace, '_desc-brainmask', obj.settings.niiFileExt];
            end
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

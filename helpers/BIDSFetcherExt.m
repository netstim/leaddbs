classdef BIDSFetcherExt
    % Extension helpers to add diffusion MRI (dwi) and resting-state fMRI (rs-fMRI)
    % discovery and BIDS-conform filename construction to an existing
    % Lead-DBS BIDS fetcher.  
    %
    %
    % Notes on BIDS naming used here (simplified, v1.9+ compatible):
    %   dMRI: sub-<id>[_ses-<id>][_acq-<label>][_dir-<label>][_run-<int>]_dwi.nii.gz 
    %         + sidecars: .json, .bval, .bvec
    %   rs-fMRI: sub-<id>[_ses-<id>][_acq-<label>][_rec-<label>][_run-<int>]_task-rest_bold.nii.gz
    %            + sidecar .json

    methods(Static)
        function out = findDiffusion(datasetDir, subjId, ses)
            % Scan BIDS tree for preop diffusion data; return paths + parsed entities
            % Also provides canonical BIDS filename for saving/renaming.
            if nargin < 3 || isempty(ses), ses = '';
            end

            base = fullfile(datasetDir, sprintf('sub-%s', subjId));
            if ~isempty(ses)
                base = fullfile(base, sprintf('ses-%s', ses));
            end
            dwiDir = fullfile(base, 'dwi');

            out = struct('present', false, 'files', [], 'entities', [], ...
                         'nii', '', 'json', '', 'bvec', '', 'bval', '', ...
                         'proposed_bids_name', '');

            if ~isfolder(dwiDir)
                return; % nothing found
            end

            % collect all *_dwi.(nii|nii.gz) files
            nii = dir(fullfile(dwiDir, sprintf('sub-%s%s*_dwi.nii*', subjId, BIDSFetcherExt.i_ses(ses))));
            if isempty(nii)
                % be more permissive: any *_dwi.nii*
                nii = dir(fullfile(dwiDir, '*_dwi.nii*'));
            end
            if isempty(nii)
                return;
            end

            % Choose first match (or implement policy for multiple runs)
            [~, idx] = max([nii.datenum]); % most recent
            f = fullfile(nii(idx).folder, nii(idx).name);
            [entities, baseStem] = BIDSFetcherExt.parseEntitiesFromStem(nii(idx).name);

            % derive companions
            json = BIDSFetcherExt.swapExt(f, '.json');
            bvec = BIDSFetcherExt.swapSuffix(f, '_dwi', '.bvec');
            bval = BIDSFetcherExt.swapSuffix(f, '_dwi', '.bval');

            out.present = true;
            out.nii = f;
            out.json = json;
            out.bvec = bvec;
            out.bval = bval;
            out.files = {f, json, bvec, bval};
            out.entities = entities;

            % Construct canonical BIDS name (keeping discovered entities)
            stem = BIDSFetcherExt.composeStem('dwi', entities);
            out.proposed_bids_name = fullfile(dwiDir, [stem '_dwi.nii.gz']);
        end

        function out = findRsfmri(datasetDir, subjId, ses)
            % Scan BIDS tree for resting-state BOLD data
            if nargin < 3 || isempty(ses), ses = '';
            end

            base = fullfile(datasetDir, sprintf('sub-%s', subjId));
            if ~isempty(ses)
                base = fullfile(base, sprintf('ses-%s', ses));
            end
            funcDir = fullfile(base, 'func');

            out = struct('present', false, 'files', [], 'entities', [], ...
                         'nii', '', 'json', '', 'proposed_bids_name', '');

            if ~isfolder(funcDir)
                return;
            end

            % Prefer explicit task-rest; accept common aliases and normalize to rest
            patt = {sprintf('sub-%s%s*_task-rest_bold.nii*', subjId, BIDSFetcherExt.i_ses(ses)), ...
                    sprintf('sub-%s%s*_task-rsfmri_bold.nii*', subjId, BIDSFetcherExt.i_ses(ses)), ...
                    sprintf('sub-%s%s*_task-restingstate_bold.nii*', subjId, BIDSFetcherExt.i_ses(ses)), ...
                    '*_task-rest*_bold.nii*'};
            hits = [];
            for k = 1:numel(patt)
                d = dir(fullfile(funcDir, patt{k}));
                if ~isempty(d), hits = [hits; d]; end %#ok<AGROW>
            end
            if isempty(hits)
                return;
            end

            [~, idx] = max([hits.datenum]);
            f = fullfile(hits(idx).folder, hits(idx).name);
            [entities, ~] = BIDSFetcherExt.parseEntitiesFromStem(hits(idx).name);
            % normalize task label to 'rest' for canonical naming
            entities.task = 'rest';

            json = BIDSFetcherExt.swapExt(f, '.json');

            out.present = true;
            out.nii = f;
            out.json = json;
            out.files = {f, json};
            out.entities = entities;

            stem = BIDSFetcherExt.composeStem('bold', entities);
            out.proposed_bids_name = fullfile(funcDir, [stem '_task-rest_bold.nii.gz']);
        end

        %% ===== Utility helpers =====
        function s = i_ses(ses)
            if isempty(ses), s = ''; else, s = ['_ses-' ses]; end
        end

        function [entities, stemNoSuffix] = parseEntitiesFromStem(filename)
            % Parse common BIDS entities from a filename (stem only)
            % Returns a struct with fields: sub, ses, acq, dir, run, rec, task
            [~, stem, ~] = fileparts(filename);
            stem = regexprep(stem, '\.nii$', ''); % in case of .nii without .gz
            parts = regexp(stem, '([^_]+)', 'match');
            entities = struct('sub','', 'ses','', 'acq','', 'dir','', 'run','', 'rec','', 'task','');
            for i=1:numel(parts)
                tok = parts{i};
                if startsWith(tok, 'sub-'), entities.sub = tok(5:end); end
                if startsWith(tok, 'ses-'), entities.ses = tok(5:end); end
                if startsWith(tok, 'acq-'), entities.acq = tok(5:end); end
                if startsWith(tok, 'dir-'), entities.dir = tok(5:end); end
                if startsWith(tok, 'run-'), entities.run = tok(5:end); end
                if startsWith(tok, 'rec-'), entities.rec = tok(5:end); end
                if startsWith(tok, 'task-'), entities.task = tok(6:end); end
            end
            stemNoSuffix = stem;
        end

        function stem = composeStem(modality, E)
            % Build a BIDS stem using discovered entities
            % modality: 'dwi' or 'bold'
            % Keep order: sub, ses, task, acq, dir/rec, run
            assert(isfield(E,'sub') && ~isempty(E.sub), 'Need sub entity');
            parts = {['sub-' E.sub]};
            if isfield(E,'ses') && ~isempty(E.ses), parts{end+1} = ['ses-' E.ses]; end
            if strcmp(modality,'bold')
                if isfield(E,'task') && ~isempty(E.task), parts{end+1} = ['task-' E.task]; end
            end
            if isfield(E,'acq') && ~isempty(E.acq), parts{end+1} = ['acq-' E.acq]; end
            if strcmp(modality,'dwi')
                if isfield(E,'dir') && ~isempty(E.dir), parts{end+1} = ['dir-' E.dir]; end
            else
                if isfield(E,'rec') && ~isempty(E.rec), parts{end+1} = ['rec-' E.rec]; end
            end
            if isfield(E,'run') && ~isempty(E.run), parts{end+1} = ['run-' E.run]; end
            stem = strjoin(parts, '_');
        end

        function out = swapExt(file, newExt)
            % swap .nii/.nii.gz -> sidecar ext
            if endsWith(file, '.nii.gz')
                out = regexprep(file, '\\.nii\\.gz$', newExt);
            else
                out = regexprep(file, '\\.nii$', newExt);
            end
        end

        function out = swapSuffix(file, suffix, newExt)
            % replace suffix (e.g., _dwi) with a different extension
            [folder, name, ~] = fileparts(file);
            name = regexprep(name, [regexptranslate('escape', suffix) '$'], '');
            out = fullfile(folder, [name suffix newExt]);
        end
    end
end
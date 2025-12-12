% ===============================
% File: ea_list_electrodes.m
% Purpose: Generic discovery of ANY number of electrodes for a patient,
% across models/hemispheres/trajectories. Returns a struct array that the
% GUI (ea_load_pts) can consume without assuming exactly two electrodes.
% ===============================
function electrodes = ea_list_electrodes(subjDir)
% EA_LIST_ELECTRODES  Enumerate electrodes for a given patient directory.
%
%   electrodes = ea_list_electrodes(subjDir)
%
% Output fields per element:
%   .id            : stable identifier (e.g., 'L1', 'R2', 'E3', etc.)
%   .label         : human-readable label
%   .side          : 'L','R', or '' if unknown
%   .model         : electrode model string (e.g., 'BSC 2202', 'Medtronic 3389')
%   .traj_index    : integer trajectory index if available
%   .num_contacts  : detected # contacts
%   .contacts      : [num_contacts x 1] struct with fields:
%                    .index     (1-based)  .name   (e.g., '1A','2','C')
%                    .position  [x y z] in current space (if available)
%   .markers       : struct with fields .head [x y z], .tail [x y z], .orient [x y z]
%   .meta          : arbitrary metadata (json, paths, etc.)
%
% This function is file-system tolerant and works with BIDS-like Lead-DBS
% derivatives as well as legacy recon folders, as long as typical patterns
% are present.

    arguments
        subjDir (1,:) char
    end

    electrodes = struct('id',{},'label',{},'side',{},'model',{}, ...
                        'traj_index',{},'num_contacts',{},'contacts',{}, ...
                        'markers',{},'meta',{});

    if ~isfolder(subjDir); return; end

    % Candidate locations for electrode info
    cand = {
        fullfile(subjDir,'reconstruction')
        fullfile(subjDir,'reconstructions')
        fullfile(subjDir,'stimulations') % sometimes markers are mirrored here
        fullfile(subjDir,'ea_reconstruction')
        fullfile(subjDir,'derivatives','leaddbs',get_basename(subjDir))
    };

    cand = cand(cellfun(@isfolder,cand));

    % Heuristics: JSON/ MAT with markers and model, contact CSV/TXT, etc.
    files = {};
    for i = 1:numel(cand)
        ff = dir(fullfile(cand{i},'**','*.*')); %#ok<DIRND>
        files = [files; fullfile({ff.folder}',{ff.name}')]; %#ok<AGROW>
    end
    files = unique(files);

    % Try to assemble electrode entries by trajectory folder (traj_*),
    % marker pairs (head/tail), or JSON groups that encode a single lead.
    trajIdx = detect_trajectories(files);

    if isempty(trajIdx)
        % Fallback: attempt a single-group parse
        e = parse_electrode_group(files, subjDir, 1);
        if ~isempty(e)
            electrodes = e;
        end
        return
    end

    k = 0;
    for t = 1:numel(trajIdx)
        e = parse_electrode_group(files, subjDir, trajIdx(t));
        if isempty(e), continue; end
        for j = 1:numel(e)
            k = k + 1; electrodes(k) = e(j); %#ok<AGROW>
        end
    end

    % Provide stable IDs and labels
    for i = 1:numel(electrodes)
        if isempty(electrodes(i).side)
            sideTag = '';
        else
            sideTag = electrodes(i).side;
        end
        electrodes(i).id    = compose_id(sideTag, electrodes(i).traj_index, i);
        electrodes(i).label = compose_label(electrodes(i));
    end
end

% ---------- helpers ----------
function nm = get_basename(p)
    [~,nm] = fileparts(p);
end

function out = compose_id(side, trajIdx, fallback)
    if ~isempty(side)
        if ~isempty(trajIdx) && ~isnan(trajIdx)
            out = sprintf('%s%d', upper(side(1)), trajIdx);
            return
        else
            out = sprintf('%s%d', upper(side(1)), fallback);
            return
        end
    end
    if ~isempty(trajIdx) && ~isnan(trajIdx)
        out = sprintf('E%d', trajIdx);
    else
        out = sprintf('E%d', fallback);
    end
end

function out = compose_label(e)
    bits = {};
    if ~isempty(e.side), bits{end+1}=upper(e.side(1)); end %#ok<AGROW>
    if ~isempty(e.traj_index) && ~isnan(e.traj_index)
        bits{end+1}=sprintf('traj-%d',e.traj_index);
    end
    if ~isempty(e.model), bits{end+1}=e.model; end
    if isempty(bits), out = 'Electrode'; else, out = strjoin(bits,' '); end
end

function idx = detect_trajectories(files)
    pat = regexptranslate('wildcard','traj_*');
    d = cellfun(@(f) ~isempty(regexp(f,pat,'once')), files);
    if any(d)
        % Extract numeric suffix
        traj = cellfun(@(f) regexp(f,'traj_(\d+)','tokens','once'), files,'uni',0);
        traj = traj(~cellfun('isempty',traj));
        if ~isempty(traj)
            idx = unique(str2double(cellfun(@(t)t{1},traj,'uni',0)));
            idx = idx(~isnan(idx));
            return
        end
    end
    idx = [];
end

function E = parse_electrode_group(files, subjDir, trajIdx)
    E = struct('id',{},'label',{},'side',{},'model',{},'traj_index',{}, ...
               'num_contacts',{},'contacts',{},'markers',{},'meta',{});

    % Try JSON first
    jsonCandidates = files(endsWith(lower(files),'.json'));
    jsonCandidates = jsonCandidates(contains(lower(jsonCandidates),'electrode') | ...
                                   contains(lower(jsonCandidates),'marker')   | ...
                                   contains(lower(jsonCandidates),'lead'));

    % Filter by trajectory index if present in path
    if ~isempty(trajIdx)
        jsonCandidates = jsonCandidates(contains(jsonCandidates, sprintf('traj_%d', trajIdx)) | ...
                                        ~contains(jsonCandidates,'traj_'));
    end

    e = struct([]);
    for i = 1:numel(jsonCandidates)
        try
            s = jsondecode(fileread(jsonCandidates{i}));
        catch
            continue
        end
        e = [e; json_to_electrodes(s, subjDir, trajIdx, jsonCandidates{i})]; %#ok<AGROW>
    end

    % If still empty, try MAT files
    if isempty(e)
        matCandidates = files(endsWith(lower(files),'.mat'));
        matCandidates = matCandidates(contains(lower(matCandidates),'electrode') | ...
                                     contains(lower(matCandidates),'marker')   | ...
                                     contains(lower(matCandidates),'lead'));
        if ~isempty(trajIdx)
            matCandidates = matCandidates(contains(matCandidates, sprintf('traj_%d', trajIdx)) | ...
                                          ~contains(matCandidates,'traj_'));
        end
        for i = 1:numel(matCandidates)
            try
                S = load(matCandidates{i});
                e = [e; mat_to_electrodes(S, subjDir, trajIdx, matCandidates{i})]; %#ok<AGROW>
            catch
                % ignore
            end
        end
    end

    % Finalize
    for i = 1:numel(e)
        if isempty(e(i).contacts)
            e(i).num_contacts = 0; 
        else
            e(i).num_contacts = numel(e(i).contacts);
        end
    end
    E = e;
end

function e = json_to_electrodes(s, subjDir, trajIdx, srcPath)
    e = struct('id',{},'label',{},'side',{},'model',{},'traj_index',{}, ...
               'num_contacts',{},'contacts',{},'markers',{},'meta',{});
    if isfield(s,'electrodes') && isstruct(s.electrodes)
        names = fieldnames(s.electrodes);
        for k = 1:numel(names)
            e(end+1) = one_from_struct(s.electrodes.(names{k}), names{k}, subjDir, trajIdx, srcPath); %#ok<AGROW>
        end
    elseif isfield(s,'model') || isfield(s,'markers')
        e(1) = one_from_struct(s, '', subjDir, trajIdx, srcPath);
    end
end

function e = mat_to_electrodes(S, subjDir, trajIdx, srcPath)
    e = struct('id',{},'label',{},'side',{},'model',{},'traj_index',{}, ...
               'num_contacts',{},'contacts',{},'markers',{},'meta',{});
    fns = fieldnames(S);
    hit = intersect(fns, {'electrodes','markers','model','side'});
    if ismember('electrodes', hit)
        E = S.electrodes;
        if isstruct(E)
            names = fieldnames(E);
            for k = 1:numel(names)
                e(end+1) = one_from_struct(E.(names{k}), names{k}, subjDir, trajIdx, srcPath); %#ok<AGROW>
            end
        end
    else
        tmp = struct();
        for i = 1:numel(hit)
            tmp.(hit{i}) = S.(hit{i});
        end
        e(1) = one_from_struct(tmp, '', subjDir, trajIdx, srcPath);
    end
end

function e = one_from_struct(s, name, subjDir, trajIdx, srcPath)
    e = struct('id','','label','','side','','model','','traj_index',trajIdx, ...
               'num_contacts',NaN,'contacts',[],'markers',struct(),'meta',struct());

    % side
    if isfield(s,'side') && ischar(s.side)
        if startsWith(lower(s.side),'l'), e.side='L';
        elseif startsWith(lower(s.side),'r'), e.side='R'; end
    elseif contains(lower(srcPath), filesep+'lh'+filesep) || contains(lower(srcPath),'left')
        e.side = 'L';
    elseif contains(lower(srcPath), filesep+'rh'+filesep) || contains(lower(srcPath),'right')
        e.side = 'R';
    end

    % model
    if isfield(s,'model') && ischar(s.model)
        e.model = s.model;
    elseif isfield(s,'electrode_model') && ischar(s.electrode_model)
        e.model = s.electrode_model;
    else
        e.model = '';
    end

    % contacts
    e.contacts = collect_contacts(s);

    % markers
    e.markers = collect_markers(s);

    % meta
    e.meta.source = srcPath;
    e.meta.subjDir = subjDir;
    if ~isempty(name), e.meta.name = name; end
end

function C = collect_contacts(s)
    C = struct('index',{},'name',{},'position',{});
    if isfield(s,'contacts') && isstruct(s.contacts)
        fn = fieldnames(s.contacts);
        for i = 1:numel(fn)
            item = s.contacts.(fn{i});
            idx = NaN; nm = fn{i}; pos = [];
            if isfield(item,'index'), idx = item.index; end
            if isfield(item,'name'),  nm  = item.name;  end
            if isfield(item,'position'), pos = item.position; end
            C(end+1) = struct('index',idx,'name',nm,'position',pos); %#ok<AGROW>
        end
        % sort by index if available
        idxs = [C.index];
        if all(~isnan(idxs))
            [~,ord] = sort(idxs); C = C(ord);
        end
    elseif isfield(s,'num_contacts') && isnumeric(s.num_contacts)
        n = s.num_contacts;
        C(n,1) = struct('index',NaN,'name','', 'position',[]);
        for i = 1:n
            C(i).index = i; C(i).name = sprintf('%d',i);
        end
    end
end

function M = collect_markers(s)
    M = struct('head',[],'tail',[],'orient',[]);
    if isfield(s,'markers') && isstruct(s.markers)
        m = s.markers;
        if isfield(m,'head'),   M.head   = m.head;   end
        if isfield(m,'tail'),   M.tail   = m.tail;   end
        if isfield(m,'orient'), M.orient = m.orient; end
    end
end





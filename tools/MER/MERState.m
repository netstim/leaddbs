classdef MERState < handle
    properties
        Config
        DBSImplants
        Frame
        MERTrajectories
        Markers
        MarkersHistory
        Toggles
        Cache
    end
    properties (SetObservable)
        Trajectory % handle to parent ea_trajectory class
    end
    properties (Constant)
        MarkerTypes = struct(...
            'Generic', 'Generic',...
            'MER', 'MER recording',...
            'LFP', 'LFP recording',...
            'Top', 'Top border',...
            'Bottom', 'Bottom border');
    end
    methods
        function setOptions(obj, options)
            % Take in options struct and keep only necessary parts in
            % obj.Config
            positions = cat(1, options.prefs.mer.tract_info.position);
            labels = {options.prefs.mer.tract_info.label};
            [side_strs, side_ids] = ea_detsidestr('both');
            clear temp
            temp(size(positions, 1) * length(side_strs)) = ...
                struct('label', '', 'offset', [], 'side', '');
            for side_ix = 1:length(side_strs)
                bLeft = strcmpi(side_strs{side_ids(side_ix)}, 'left');
                for pos_ix = 1:size(positions, 1)
                    arr_ix = (side_ix - 1) * size(positions, 1) + pos_ix;
                    temp(arr_ix).label = labels{pos_ix};
                    temp(arr_ix).side = side_strs{side_ix};
                    temp(arr_ix).offset = positions(pos_ix, :) * options.prefs.mer.offset;
                    if bLeft
                        temp(arr_ix).offset(1) = -temp(arr_ix).offset(1);
                    end
                end
            end
            obj.Config.MERTrajectory = temp;
            obj.Config.ImplantLabel = labels{options.prefs.mer.defaulttract};
            obj.Config.MERLength = options.prefs.mer.length;
            obj.Config.MERPnts = options.prefs.mer.n_pnts;
            % For now, some we will just copy directly to struct root
            obj.Config.vis.step_size = options.prefs.mer.step_size;
            obj.Config.vis.tag_visible = options.prefs.mer.tag.visible;
            obj.Config.vis.markersize = options.prefs.mer.markersize;
            obj.Config.uipatdirs = options.uipatdirs;
            obj.Config.native = options.native;
            if ~isfield(options, 'loadnativereco')
                obj.Config.loadnativereco = 0;
            else
                obj.Config.loadnativereco = options.loadnativereco;
            end
            obj.Config.root = options.root;
            obj.Config.patientname = options.patientname;
            obj.Config.sides = options.sides;
            obj.Config.prefs.prenii_searchstring = options.prefs.prenii_searchstring;
            obj.Config.prefs.prenii_order = options.prefs.prenii_order;
            obj.Config.prefs.prenii_unnormalized = options.prefs.prenii_unnormalized;
            obj.Config.elmodel = options.elmodel;
            obj.Config.elspec = options.elspec;
        end
        function clearData(obj)
            obj.Cache = struct('prenii_mat', [], 'native2mni_emp_mat', []);
            obj.DBSImplants = struct(...
                'side', {}, 'depth', {},...
                'implanted_tract_label', {},...
                'coords', {}, 'coords_bottom', {}, 'elec_uv', {});
            obj.Frame = struct('side', {'right', 'left'}, 'yaw_rad', 0,...
                'landmarks', {struct('label', {}, 'coords', {})});
            obj.MERTrajectories = struct('side', {},...
                'label', {}, 'depth', {}, 'translation', {});
            obj.Toggles = struct(...
                'keycontrol', struct('side', {}, 'label', {}, 'value', {}),...
                'togglestates', struct('side', {}, 'label', {}, 'value', {}));
            obj.setMarkersToDefaults();  % Defaults same as clearing.
        end
        function setDataToDefaults(obj)
            obj.setDBSImplantsToDefaults();
            obj.setMERTrajectoriesToDefaults();
            obj.setFrameToDefaults();
            obj.setTogglesToDefaults();
            obj.setMarkersToDefaults();
        end
        function setDBSImplantsToDefaults(obj)
            [side_strs, side_ids] = ea_detsidestr('both');
            for side_ix = 1:length(side_ids)
                sstr = side_strs{side_ids(side_ix)};
                obj.DBSImplants(side_ix).side = sstr;
                obj.DBSImplants(side_ix).depth = 0;
                obj.DBSImplants(side_ix).implanted_tract_label = obj.Config.ImplantLabel;
            end
            obj.loadDBSReconstruction();
        end
        function setFrameToDefaults(obj)
            [side_strs, ~] = ea_detsidestr('both');
            def_frame_landmarks = {'A', 'E', 'Entry'};
            def_frame_coords = [0, 0.5, 0; 0, -0.5, 0; 0, 0, -1];
            landmarks = struct('label', def_frame_landmarks,...
                'coords', mat2cell(def_frame_coords, [1 1 1], 3)');
            % Both sides initialized with same default landmarks.
            for side_ix = 1:length(side_strs)
                obj.updateFrame(side_strs{side_ix}, landmarks);
            end
        end
        function setMarkersToDefaults(obj)
            % No defaults, just clear.
            obj.Markers = struct('side', {}, 'tract_label', {},...
                'depth', {}, 'session', {}, 'type', {}, 'notes', {});
            obj.MarkersHistory = obj.Markers;
        end
        function setMERTrajectoriesToDefaults(obj)
            for traj_ix = 1:length(obj.Config.MERTrajectory)
                obj.MERTrajectories(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.MERTrajectories(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.MERTrajectories(traj_ix).depth = 0;
            end
            % Calc .translation
            obj.calculateMERTranslations();
        end
        function setTogglesToDefaults(obj)
            for traj_ix = 1:length(obj.Config.MERTrajectory)
                obj.Toggles.keycontrol(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.Toggles.keycontrol(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.Toggles.keycontrol(traj_ix).value = 0;

                obj.Toggles.togglestates(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.Toggles.togglestates(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.Toggles.togglestates(traj_ix).value = 1;
            end
        end
        function loadDBSReconstruction(obj)
            opt_native_backup = obj.Config.native;
            obj.Config.native = 1;
            if isempty(obj.Trajectory) % load reconstruction data from patient folder and configure
                [~, ~, dbs_contacts, obj.Config.elmodel] = ea_load_reconstruction(obj.Config);
            else % trajectory object supplied - will relate MER trajectories to information of the object
                obj.Config.elmodel=obj.Trajectory.elmodel;

                switch obj.Trajectory.relateMicro
                    case 'macro' % relate MER fiducials to DBS electrode reconstructed from postoperative data
                        dbs_contacts=obj.Trajectory.elstruct.markers;
                    case 'planning' % relate MER fiducials to planning trajectory
                        % build markers struct from planning fiducial line:
                        dbs_contacts(1).head=obj.Trajectory.target.target;
                        dbs_contacts(1).tail=obj.Trajectory.target.entry;
                        opts.elmodel=obj.Trajectory.elmodel;
                        opts=ea_resolve_elspec(opts);
                        el=load([ea_getearoot,'templates',filesep,'electrode_models',filesep,opts.elspec.matfname,'.mat']);
                        hdist=pdist([el.electrode.head_position;el.electrode.tail_position]);
                        dbs_contacts(1).tail=dbs_contacts(1).head+...
                            ((dbs_contacts(1).tail-dbs_contacts(1).head)/...
                            norm(dbs_contacts(1).tail-dbs_contacts(1).head))*...
                            hdist;
                        [xunitv, yunitv] = ea_calcxy(dbs_contacts(1).head, dbs_contacts(1).tail);
                        dbs_contacts(1).x = dbs_contacts(1).head + xunitv*(opts.elspec.lead_diameter/2);
                        dbs_contacts(1).y = dbs_contacts(1).head + yunitv*(opts.elspec.lead_diameter/2);
                        % end build markers struct from planning fiducial line.
                end
            end

            obj.Config = ea_resolve_elspec(obj.Config);
            dbs_contacts = ea_resolvecoords(dbs_contacts, obj.Config);
            % Get template space
            obj.Config.native = 0;
            [~, ~, dbs_markers_mni] = ea_load_reconstruction(obj.Config);
            [dbs_contacts_mni, ~, ~] = ea_resolvecoords(dbs_markers_mni, obj.Config);
            for side_ix = 1:length(dbs_contacts)
                if ~isempty(dbs_contacts{side_ix})
                    coords = dbs_contacts{side_ix};
                    % Average 3-D distance through electrode contacts.
                    elec_diff = mean(diff(coords));
                    % unit vector through contacts from lowest (?) contact
                    elec_uv = elec_diff / norm(elec_diff);
                    % shift coordinates to bottom of lowest contact.
                    coords_bottom = bsxfun(@minus, coords,...
                        elec_uv * obj.Config.elspec.contact_length / 2);
                    obj.DBSImplants(side_ix).coords = coords;
                    obj.DBSImplants(side_ix).coords_bottom = coords_bottom;
                    obj.DBSImplants(side_ix).elec_uv = elec_uv;
                end
                if ~isempty(dbs_contacts_mni{side_ix})
                    coords = dbs_contacts_mni{side_ix};
                    elec_diff = mean(diff(coords));
                    elec_uv = elec_diff / norm(elec_diff);
                    obj.DBSImplants(side_ix).coords_mni = coords;
                    obj.DBSImplants(side_ix).elec_uv_mni = elec_uv;
                end
            end
            obj.Config.native = opt_native_backup;
        end
        function updateDBSImplantTrack(obj, side, label)
            bSide = strcmpi({obj.DBSImplants.side}, side);
            bTraj = strcmpi({obj.MERTrajectories.side}, side) ...
                & strcmpi({obj.MERTrajectories.label}, label);
            if any(bSide) && any(bTraj)
                obj.DBSImplants(bSide).implanted_tract_label = label;
                obj.calculateMERTranslations();
            else
                warning('DBS implanted_tract_label not updated because matching MER track not found.');
            end
        end
        function updateDBSDepth(obj, side, depth)
            obj.DBSImplants(strcmpi({obj.DBSImplants.side}, side)).depth = depth;
            bTraj = strcmpi({obj.MERTrajectories.side}, side);
            [obj.MERTrajectories(bTraj).depth] = deal(depth);
        end
        function updateFrame(obj, side, landmarks)
            bSide = strcmpi({obj.Frame.side}, side);
            if ~any(bSide)
                obj.Frame(end + 1).side = side;
                bSide(end+1) = true;
            end
            obj.Frame(bSide).landmarks = landmarks;

            % Caclulate frame rotation. Currently only A/E/Entry works.
            % TODO: Support other configurations than requiring A, E, and Entry.
            [~, ~, ib] = intersect({'A', 'E', 'Entry'}, {landmarks.label}, 'stable');
            AEP = cat(1, landmarks(ib).coords);
            m = (AEP(1, :) + AEP(2, :)) / 2;  %midpoint between A and E becomes origin.
            y = AEP(1, :) - m;  % Vector from frame origin to frame +y in scrf space
            z = m - AEP(3, :);  % '' for +z
            x = cross(y, z);    % '' for +x
            x = x / norm(x); y = y / norm(y); z = z / norm(z);  % Make unit vectors.

            % Get the transformation between native/scrf and frame spaces
            def_spc = [1 0 0; 0 1 0; 0 0 1];
            xform = [x; y; z] \ def_spc;  % == inv([x; y; z])
            % xform gives full nexframe transform.
            % However, DBS lead localization constrains pitch and roll,
            % so we only need yaw. Save yaw in merframe
            obj.Frame(bSide).yaw_rad = atan2(xform(2, 1), xform(1, 1));
            obj.calculateMERTranslations();
        end
        function calculateMERTranslations(obj)
            uqsides = unique({obj.Config.MERTrajectory.side}, 'stable');
            for side_ix = 1:length(uqsides)
                side_str = uqsides{side_ix};
                b_side_cfg = strcmpi({obj.Config.MERTrajectory.side}, side_str);
                cfg_traj = obj.Config.MERTrajectory(b_side_cfg);
                cfg_labels = {cfg_traj.label};
                side_transl = cat(1, cfg_traj.offset);  % In native space
                b_side_impl = strcmpi({obj.DBSImplants.side}, side_str);
                impl_label = obj.DBSImplants(b_side_impl).implanted_tract_label;
                % Offset by the MER tract the DBS electrode was implanted in.
                side_transl = bsxfun(@minus, side_transl,...
                    side_transl(strcmpi(cfg_labels, impl_label), :));
                % Apply yaw rotation.
                b_side_frame = strcmpi({obj.Frame.side}, side_str);
                yaw_rad = obj.Frame(b_side_frame).yaw_rad;
                yaw_xform = [cos(yaw_rad) -sin(yaw_rad) 0; sin(yaw_rad) cos(yaw_rad) 0; 0 0 1];
                side_transl = side_transl * yaw_xform;
                % Align to DBS electrode for pitch & roll
                elec_uv = obj.DBSImplants(b_side_impl).elec_uv;
                if ~isequal(elec_uv, [0 0 1])
                    % Find the transform that aligns [0 0 1] with the electrode unit
                    % vector, then apply that transform to each of the MER offsets.
                    U = MERState.align_vectors([0 0 1], elec_uv);
                    side_transl = (U*side_transl')';  % pitch & roll
                end
                % Assign side_transl back into obj.MERTrajectories.translation
                b_side_traj = strcmpi({obj.MERTrajectories.side}, side_str);
                traj_labels = {obj.MERTrajectories.label};
                for traj_ix = 1:length(cfg_traj)
                    b_label = strcmpi(traj_labels, cfg_traj(traj_ix).label);
                    obj.MERTrajectories(b_side_traj & b_label).translation = side_transl(traj_ix, :);
                end
            end
        end
        function updateTrajDepth(obj, side, track, depth)
            bTraj = strcmpi({obj.MERTrajectories.side}, side) ...
                & strcmpi({obj.MERTrajectories.label}, track);
            obj.MERTrajectories(bTraj).depth = depth;
        end
        function coords = getMERTrajectory(obj, traj, spc)
            bSid = strcmpi({obj.DBSImplants.side}, traj.side);
            curr_dist = traj.depth - obj.DBSImplants(bSid).depth;
            im_mm = obj.DBSImplants(bSid).coords_bottom(1, :);
            elec_uv = obj.DBSImplants(bSid).elec_uv;
            startpoint = im_mm + traj.translation + (elec_uv .* curr_dist);
            if strcmpi(spc, 'mni')
                startpoint = obj.native2mni_fast(startpoint, traj.side);
                endpoint = startpoint + obj.DBSImplants(bSid).elec_uv_mni * obj.Config.MERLength;
            else
                endpoint = startpoint + elec_uv * obj.Config.MERLength;
            end
            stepsize = (endpoint - startpoint) / (obj.Config.MERPnts - 1);
            coords = bsxfun(@plus, (0:obj.Config.MERPnts - 1)' * stepsize, startpoint);
        end
        function addMarkerAtDepth(obj, side, label, type, sess_notes, depth)
            bMarkers = strcmpi({obj.Markers.side}, side)...
                & strcmpi({obj.Markers.tract_label}, label)...
                & ([obj.Markers.depth] == depth);
            if any(bMarkers)
                warning('Location along %s - %s at depth %f is already marked as type %s.',...
                    side, label, depth, obj.Markers(bMarkers).type);
            else
                obj.Markers(end + 1).side = side;
                obj.Markers(end).tract_label = label;
                obj.Markers(end).depth = depth;
                obj.Markers(end).type = type;
                obj.Markers(end).session = sess_notes;
            end
        end
        function addMarkersAtTrajs(obj, tstruct, type, sess_notes)
            for traj_ix = 1:length(tstruct)
                ts = tstruct(traj_ix);
                bTraj = strcmpi({obj.MERTrajectories.side}, ts.side) ...
                    & strcmpi({obj.MERTrajectories.label}, ts.label);
                if any(bTraj)
                    traj = obj.MERTrajectories(bTraj);
                    obj.addMarkerAtDepth(ts.side, ts.label, type, sess_notes, traj.depth);
                else
                    fprintf('No MER trajectory found for marker %s - %s.\n', side, label);
                    fprintf('Try running obj.setMERTrajectoriesToDefaults() first.\n');
                end
            end
        end
        function marker_table = exportMarkers(obj)
            col_headers = {'Patient', 'Side', 'Tract', 'Type',...
                'Native_X', 'Native_Y', 'Native_Z',...
                'MNI_X', 'MNI_Y', 'MNI_Z'};
            pts = repmat({obj.Config.patientname}, length(obj.Markers), 1);
            coords_native = nan(length(obj.Markers), 3);
            for m_ix = 1:length(obj.Markers)
                coords_native(m_ix, :) = obj.getMarkerPosition(obj.Markers(m_ix), 'native');
            end
            coords_mni = obj.native2mni_slow(coords_native);
            marker_table = table(...
                pts, {obj.Markers.side}', {obj.Markers.tract_label}', {obj.Markers.type}',...
                coords_native(:, 1), coords_native(:, 2), coords_native(:, 3),...
                coords_mni(:, 1), coords_mni(:, 2), coords_mni(:, 3),...
                'VariableNames', col_headers);
        end
        function coords = getMarkerPosition(obj, marker, spc)
            bTraj = strcmpi({obj.MERTrajectories.side}, marker.side)...
                & strcmpi({obj.MERTrajectories.label}, marker.tract_label);
            if any(bTraj)
                % Calculate coordinate as translation from DBS electrode.
                bSid = strcmpi({obj.DBSImplants.side}, marker.side);
                curr_dist = marker.depth - obj.DBSImplants(bSid).depth;
                traj = obj.MERTrajectories(bTraj);
                im_mm = obj.DBSImplants(bSid).coords_bottom(1, :);
                elec_uv = obj.DBSImplants(bSid).elec_uv;
                coords = im_mm + traj.translation + (elec_uv .* curr_dist);
                if strcmpi(spc, 'mni')
                    coords = obj.native2mni_fast(coords, marker.side);
                end
            else
                warning('No MER trajectory found for marker %s - %s.', marker.side, marker.label);
                fprintf('Try running obj.setMERTrajectoriesToDefaults() first.\n');
                coords = nan(1, 3);
            end
        end
        function undoMarker(obj)
            old = obj.Markers(end);
            obj.MarkersHistory(end+1) = old;
            obj.Markers(end) = [];
        end
        function obj = redoMarker(obj)
            obj.Markers(end + 1) = obj.MarkersHistory(end);
            obj.MarkersHistory(end) = [];
        end
        function translateToggledTrajectories(obj, d)
            %obj = translateToggledTrajectories(obj, d)
            %moves toggled trajectories along their lengths by (+/-) d mm.
            d = sign(d) * obj.Config.vis.step_size(abs(d));
            tog_ids = find([obj.Toggles.keycontrol.value] == 1);
            for tog_ix = 1:length(tog_ids)
                kc = obj.Toggles.keycontrol(tog_ids(tog_ix));
                mer_ids = find(strcmpi({obj.MERTrajectories.side}, kc.side)...
                    & strcmpi({obj.MERTrajectories.label}, kc.label));
                for mer_ix = 1:length(mer_ids)
                    mid = mer_ids(mer_ix);
                    obj.MERTrajectories(mid).depth = obj.MERTrajectories(mid).depth + d;
                end
            end
        end
        function save(obj, varargin)
            fpath = fullfile(obj.Config.root, obj.Config.patientname, 'ea_merstate.mat');
            save_ok = true;
            if exist(fpath, 'file') && ...
                    ((nargin < 2) || varargin{1}(1)~='y')
                overwrite = questdlg({...
                    ['ea_merstate.mat found in ' obj.Config.patientname ' directory.'];...
                    'Would you like to overwrite this file?'}, obj.Config.patientname);
                save_ok = strcmpi(overwrite,'Yes');
            end
            config = obj.Config; %#ok<NASGU>
            dbsimplants = obj.DBSImplants; %#ok<NASGU>
            frame = obj.Frame; %#ok<NASGU>
            markers = obj.Markers; %#ok<NASGU>
            toggles = obj.Toggles; %#ok<NASGU>
            if save_ok
                disp(['Saving: ', fpath]);
                save(fpath, 'config', 'dbsimplants', 'frame', 'markers', 'toggles');
                disp('DONE');
            end
        end
        function load(obj)
            fpath = fullfile(obj.Config.root, obj.Config.patientname, 'ea_merstate.mat');
            if exist(fpath, 'file')
                temp = load(fpath);
                obj.clearData();
                obj.Config = temp.config;
                obj.DBSImplants = temp.dbsimplants;
                obj.Frame = temp.frame;
                obj.Markers = temp.markers;
                obj.Toggles = temp.toggles;
                obj.loadDBSReconstruction();  % In case reco changed since last save.
                obj.setMERTrajectoriesToDefaults();
            else
                disp(['File not found: ' fpath]);
            end
        end
    end
    methods(Access = private)
        function coords_mni = native2mni_fast(obj, coords_native, side)
            % Use native2mni_fast for visualization, native2mni_slow for
            % reporting values accurately.
            if ~isfield(obj.Cache, 'prenii_mat')...
                    || isempty(obj.Cache.prenii_mat)
                options = ea_assignpretra(obj.Config);
                ptdir = fullfile(obj.Config.root, obj.Config.patientname);
                prenii_fname = fullfile(ptdir, options.prefs.prenii_unnormalized);
                nii = ea_load_nii(prenii_fname);
                obj.Cache.prenii_mat = nii(1).mat;
            end
            if ~isfield(obj.Cache, 'native2mni_emp_mat')...
                    || isempty(obj.Cache.native2mni_emp_mat)
                options = ea_assignpretra(obj.Config);
                ptdir = fullfile(obj.Config.root, obj.Config.patientname);
                prenii_fname = fullfile(ptdir, options.prefs.prenii_unnormalized);
                vx_native = cell(1, length(obj.DBSImplants));
                mm_mni = cell(1, length(obj.DBSImplants));
                dbs_coords = arrayfun(@(x)x.coords, obj.DBSImplants, 'UniformOutput', false);
                for sid = 1:length(obj.DBSImplants)
                    % Generate an array of native coordinates spanning the
                    % DBS lead extents
                    dimvec = cell(1, 3);
                    for dim_ix = 1:3
                        this_coords = dbs_coords{sid}(:, dim_ix);
                        this_lims = [min(this_coords), max(this_coords)];
                        this_span = abs(diff(this_lims));
                        this_lims = [this_lims(1)-0.1*this_span this_lims(2)+0.1*this_span];
                        dimvec{dim_ix} = this_lims(1):0.2:this_lims(end);
                        if length(dimvec{dim_ix}) < 10
                            dimvec{dim_ix} = linspace(this_lims(1), this_lims(2), 10);
                        end
                    end
                    [X, Y, Z] = meshgrid(dimvec{:});
                    XYZ_nii_mm = [X(:)'; Y(:)'; Z(:)'; ones(1, numel(X))];
                    % Convert to voxels (voxels are assumed by ea_map_coords)
                    XYZ_nii_vx = obj.Cache.prenii_mat \ XYZ_nii_mm;
                    % Map to template space. Slow, but only once per side.
                    XYZ_mni_mm = ea_map_coords(XYZ_nii_vx(1:3,:), prenii_fname, ...
                        fullfile(ptdir, 'inverseTransform'), '');
                    % Save some values for later.
                    vx_native{sid} = XYZ_nii_vx;
                    mm_mni{sid} = XYZ_mni_mm;
                    % Calculate emperical transform
                    voxnii2mmnorm = (XYZ_nii_vx(1:4, :)' \ [XYZ_mni_mm; ones(1, size(XYZ_mni_mm, 2))]')';
                    % Combine with mm2vx so now it is mm_native_2_mm_mni
                    obj.Cache.native2mni_emp_mat.(obj.DBSImplants(sid).side) = ...
                        voxnii2mmnorm / obj.Cache.prenii_mat; % combine mats
                end
                % Calculate combined transform for when side is unknown.
                XYZ_nii_vx = cat(2, vx_native{:});
                XYZ_mni_mm = cat(2, mm_mni{:});
                voxnii2mmnorm = (XYZ_nii_vx(1:4, :)' \ [XYZ_mni_mm; ones(1, size(XYZ_mni_mm, 2))]')';
                obj.Cache.native2mni_emp_mat.both = ...
                        voxnii2mmnorm / obj.Cache.prenii_mat;
            end
            coords_native = [coords_native, ones(size(coords_native, 1), 1)]';
            if ~exist('side','var')
                side = 'both';
            end
            coords_mni = (obj.Cache.native2mni_emp_mat.(side) * coords_native)';
            coords_mni = coords_mni(:, 1:3);
        end
        function coords_mni = native2mni_slow(obj, coords_native)
            % Use native2mni_fast for visualization, native2mni_slow for
            % reporting values accurately.
            options = ea_assignpretra(obj.Config);
            ptdir = fullfile(obj.Config.root, obj.Config.patientname);
            prenii_fname = fullfile(ptdir, options.prefs.prenii_unnormalized);
            if ~isfield(obj.Cache, 'prenii_mat')...
                    || isempty(obj.Cache.prenii_mat)
                nii = ea_load_nii(prenii_fname);
                obj.Cache.prenii_mat = nii(1).mat;
            end
            cmm = [coords_native, ones(size(coords_native, 1), 1)]';
            cvx = obj.Cache.prenii_mat \ cmm;  % mm2vx
            coords_mni = ea_map_coords(cvx(1:3, :), prenii_fname, ...
                fullfile(ptdir, 'inverseTransform'), '')';
        end
    end
    methods(Static = true)
        function RU = align_vectors(A, B)
            cross_AB = cross(A,B);
            ssc_cross_AB = [...
                0 -cross_AB(3) cross_AB(2);...
                cross_AB(3) 0 -cross_AB(1);...
                -cross_AB(2) cross_AB(1) 0];
            RU = eye(3) + ssc_cross_AB + ssc_cross_AB^2*(1-dot(A,B))/(norm(cross_AB)^2);
        end
    end
end

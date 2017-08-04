classdef MERState
    properties
        Config
        DBSImplants
        Frame
        MERTrajectories
        Markers
        MarkersHistory
        Toggles
    end
    methods
%         function obj = MERState(options)
%             obj = obj.setOptions(options);
%             obj = obj.clearData();
%         end
        function obj = setOptions(obj, options)
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
            obj.Config.loadrecoforviz = options.loadrecoforviz;
            obj.Config.root = options.root;
            obj.Config.patientname = options.patientname;
            obj.Config.sides = options.sides;
            obj.Config.prefs.prenii_searchstring = options.prefs.prenii_searchstring;
            obj.Config.prefs.prenii_order = options.prefs.prenii_order;
            obj.Config.prefs.prenii_unnormalized = options.prefs.prenii_unnormalized;
            obj.Config.elmodel = options.elmodel;
            obj.Config.elspec = options.elspec;            
        end
        function obj = clearData(obj)
            obj.DBSImplants = struct(...
                'side', {}, 'depth', {},...
                'implanted_tract_label', {},...
                'coords', struct('native', [], 'mni', []));
            obj.Frame = struct('side', {}, 'yaw_rad', {},...
                'landmarks', {struct('label', {}, 'coords', struct('native', []))});
            obj.MERTrajectories = struct('side', {},...
                'coords', struct('native', [], 'mni', []),...
                'label', {}, 'depth', {}, 'translation', {});
            obj.Toggles = struct(...
                'keycontrol', struct('side', {}, 'label', {}, 'value', {}),...
                'togglestates', struct('side', {}, 'label', {}, 'value', {}));
            obj = obj.setMarkersToDefaults();  % Defaults same as clearing.
        end
        function obj = setDataToDefaults(obj)
            obj = obj.setDBSImplantsToDefaults();
            obj = obj.setFrameToDefaults();
            obj = obj.setMERTrajectoriesToDefaults();
            obj = obj.setTogglesToDefaults();
            obj = obj.setMarkersToDefaults();
        end
        function obj = setDBSImplantsToDefaults(obj)
            [side_strs, side_ids] = ea_detsidestr('both');
            for side_ix = 1:length(side_ids)
                sstr = side_strs{side_ids(side_ix)};
                obj.DBSImplants(side_ix).side = sstr;
                obj.DBSImplants(side_ix).depth = 0;
                obj.DBSImplants(side_ix).implanted_tract_label = obj.Config.ImplantLabel;
            end
            obj = obj.loadDBSReconstruction();
        end
        function obj = setFrameToDefaults(obj)
            [side_strs, ~] = ea_detsidestr('both');
            def_frame_landmarks = {'A', 'E', 'Entry'};
            def_frame_coords = [0, 0.5, 0; 0, -0.5, 0; 0, 0, -1];
            landmarks = struct('label', def_frame_landmarks);
            for lm_ix = 1:length(def_frame_landmarks)
                landmarks(lm_ix).coords.native = def_frame_coords(lm_ix, :);
            end
            % Both sides initialized with same default coordinates/labels.
            for side_ix = 1:length(side_strs)
                obj = obj.updateFrame(side_strs{side_ix}, landmarks);
            end
        end
        function obj = setMarkersToDefaults(obj)
            % No defaults, just clear.
            obj.Markers = struct('side', {}, 'tract_label', {},...
                'depth', {}, 'session', {}, 'type', {}, 'notes', {},...
                'coords', struct('native', [], 'mni', []));
            obj.MarkersHistory = obj.Markers;
        end
        function obj = setMERTrajectoriesToDefaults(obj)
            for traj_ix = 1:length(obj.Config.MERTrajectory)
                obj.MERTrajectories(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.MERTrajectories(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.MERTrajectories(traj_ix).depth = 0;
            end
            % Calc .translation and .coords
            obj = obj.calculateMERTranslations();
        end
        function obj = setTogglesToDefaults(obj)
            for traj_ix = 1:length(obj.Config.MERTrajectory)
                obj.Toggles.keycontrol(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.Toggles.keycontrol(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.Toggles.keycontrol(traj_ix).value = 0;
                
                obj.Toggles.togglestates(traj_ix).side = obj.Config.MERTrajectory(traj_ix).side;
                obj.Toggles.togglestates(traj_ix).label = obj.Config.MERTrajectory(traj_ix).label;
                obj.Toggles.togglestates(traj_ix).value = 1;
            end
        end
        function obj = updateFrame(obj, side, landmarks)
            bSide = strcmpi({obj.Frame.side}, side);
            if ~any(bSide)
                obj.Frame(end + 1).side = side;
                bSide(end+1) = true;
            end
            obj.Frame(bSide).landmarks = landmarks;
            
            % Caclulate frame rotation. Currently only A/E/Entry works.
            % TODO: Support other configurations than requiring A, E, and Entry.
            AEP_str = {'A', 'E', 'Entry'};
            AEP = nan(3, 3);
            for aep_ix = 1:3
                b_row = strcmpi({landmarks.label}, AEP_str{aep_ix});
                AEP(aep_ix, :) = landmarks(b_row).coords.native;
            end
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
        end
        function obj = loadDBSReconstruction(obj)
            opt_native_backup = obj.Config.native;
            obj.Config.native = 1;
            [~, ~, dbs_contacts, obj.Config.elmodel] = ea_load_reconstruction(obj.Config);
            obj.Config = ea_resolve_elspec(obj.Config);
            dbs_contacts = ea_resolvecoords(dbs_contacts, obj.Config);
            obj = obj.populateImplantCoords(dbs_contacts, 'native');
            if ~opt_native_backup
                obj.Config.native = 0;
                dbs_contacts_mni = ea_load_reconstruction(obj.Config);
                % TODO: Why no ea_resolvecoords for mni space?
                obj = obj.populateImplantCoords(dbs_contacts_mni, 'mni');
            end
            obj.Config.native = opt_native_backup;
        end
        function obj = calculateMERTranslations(obj)
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
                elec_uv = obj.DBSImplants(b_side_impl).elec_uv.native;
                if ~isequal(elec_uv, [0 0 1])
                    % Find the transform that aligns [0 0 1] with the electrode unit
                    % vector, then apply that transform to each of the MER offsets.
                    U = MERState.align_vectors([0 0 1], elec_uv);
                    side_transl = (U*side_transl')';  % pitch & roll
                end
                % Asign side_transl back into obj.MERTrajectories.translation
                b_side_traj = strcmpi({obj.MERTrajectories.side}, side_str);
                traj_labels = {obj.MERTrajectories.label};
                for traj_ix = 1:length(cfg_traj)
                    b_label = strcmpi(traj_labels, cfg_traj(traj_ix).label);
                    obj.MERTrajectories(b_side_traj & b_label).translation.native = side_transl(traj_ix, :);
                end
            end
            
            if ~obj.Config.native
                transl_native = arrayfun(@(x) x.translation.native, obj.MERTrajectories, 'UniformOutput', false);
                transl_native = cat(1, transl_native{:});
                startpoint_native = nan(size(transl_native));
                for side_ix = 1:length(obj.DBSImplants)
                    b_traj = strcmpi({obj.MERTrajectories.side}, obj.DBSImplants(side_ix).side);
                    startpoint_native(b_traj, :) = bsxfun(@plus, transl_native(b_traj, :),...
                        obj.DBSImplants(side_ix).coords_top.native(1, :));
                end
                
                ptdir = fullfile(obj.Config.root, obj.Config.patientname);
                options = ea_assignpretra(obj.Config);
                prenii_fname = fullfile(ptdir, options.prefs.prenii_unnormalized);
                nii = ea_load_nii(prenii_fname);
                
                c = [startpoint_native, ones(size(startpoint_native, 1), 1)]';
                c = nii(1).mat \ c;
                %         try
                %             whichnormmethod = ea_whichnormmethod(ptdir);
                %             if ~ismember(whichnormmethod, ea_getantsnormfuns)
                %                 V = spm_vol(fullfile(ptdir, 'y_ea_inv_normparams.nii'));
                %                 if ~isequal(V.dim, nii.dim)
                %                     ea_redo_inv(ptdir, options);
                %                 end
                %             end
                %         end
                startpoint_mni = ea_map_coords(c(1:3, :), prenii_fname, ...
                    fullfile(ptdir, 'y_ea_inv_normparams.nii'), '')';
                transl_mni = nan(size(startpoint_mni));
                for side_ix = 1:length(obj.DBSImplants)
                    b_traj = strcmpi({obj.MERTrajectories.side}, obj.DBSImplants(side_ix).side);
                    transl_mni(b_traj, :) = bsxfun(@minus, startpoint_mni(b_traj, :),...
                        obj.DBSImplants(side_ix).coords_top.mni(1, :));
                end
                for traj_ix = 1:length(obj.MERTrajectories)
                    obj.MERTrajectories(traj_ix).translation.mni = ...
                        transl_mni(traj_ix, :);
                end
            end
            
            % Everytime the translation is updated, the trajectories need
            % to be udpated too.
            obj = obj.calculateMERTrajectories();
        end
        function obj = calculateMERTrajectories(obj)
            % Collect common values used for each trajectory
            if obj.Config.native
                spc = 'native';
            else
                spc = 'mni';
            end
            % DBS coordinates and unit vector for each hemisphere
            side_strs = {obj.DBSImplants.side};
            im_depths = [obj.DBSImplants.depth];
            
            % Calculate trajectory as translation from DBS electrode.
            for traj_ix = 1:length(obj.MERTrajectories)
                traj = obj.MERTrajectories(traj_ix);
                bSid = strcmpi(side_strs, traj.side);
                curr_dist = traj.depth - im_depths(bSid);
                im_mm = obj.DBSImplants(bSid).coords_top.(spc)(1, :);
                elec_uv = obj.DBSImplants(bSid).elec_uv.(spc);
                startpoint = im_mm + traj.translation.(spc) + (elec_uv .* curr_dist);
                endpoint = startpoint + elec_uv * obj.Config.MERLength;
                stepsize = (endpoint - startpoint) / (obj.Config.MERPnts - 1);
                obj.MERTrajectories(traj_ix).coords.(spc) = ...
                    bsxfun(@plus, (0:obj.Config.MERPnts - 1)' * stepsize, startpoint);
            end
        end
        function obj = addMarker(obj, side, label, type, sess_notes)
            bTraj = strcmpi({obj.MERTrajectories.side}, side) ...
                & strcmpi({obj.MERTrajectories.label}, label);
            if any(bTraj)
                traj = obj.MERTrajectories(bTraj);
                bMarkers = strcmpi({obj.Markers.side}, side)...
                    & strcmpi({obj.Markers.tract_label}, label)...
                    & ([obj.Markers.depth] == traj.depth);
                if any(bMarkers)
                    fprintf('Location along %s - %s at depth %f is already marked as type %s.\n',...
                            side, label, traj.depth, obj.Markers(bMarkers).type);
                else
                    %side, tract_label, depth, session, type, notes,
                    %coords.native/coords.mni
                    obj.Markers(end + 1).side = side;
                    obj.Markers(end).tract_label = label;
                    obj.Markers(end).depth = traj.depth;
                    obj.Markers(end).type = type;
                    obj.Markers(end).session = sess_notes;
                    for fn_cell = fieldnames(traj.coords)
                        fn = fn_cell{:};
                        obj.Markers(end).coords.(fn) = traj.coords.(fn)(1, :);
                    end
                end
            else
                fprintf('No MER trajectory found for marker %s - %s.\n', side, label);
            end
        end
        function obj = undoMarker(obj)
            old = obj.Markers(end);
            obj.MarkersHistory(end+1) = old;
            obj.Markers(end) = [];
        end
        function obj = redoMarker(obj)
            obj.Markers(end + 1) = obj.MarkersHistory(end);
            obj.MarkersHistory(end) = [];
        end
        function obj = translateTrajectories(obj, d)
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
            obj = obj.calculateMERTrajectories();
        end
        function save(obj)
            fpath = fullfile(obj.Config.uipatdirs{1}, 'ea_merstate.mat');
            save_ok = true;
            if exist(fpath, 'file')
                overwrite = ea_questdlg({['ea_merstate.mat found in ' obj.Config.patientname ' directory.'];...
                    'Would you like to overwrite this file?'}, obj.Config.patientname);
                save_ok = strcmpi(overwrite,'Yes');
            end
            config = obj.Config;
            dbsimplants = obj.DBSImplants;
            frame = obj.Frame;
            markers = obj.Markers;
            toggles = obj.Toggles;
            if save_ok
                disp(['Saving: ', fpath]);
                save(fpath, 'config', 'dbsimplants', 'frame', 'markers', 'toggles');
                disp('DONE');
            end
        end
        function obj = load(obj)
            fpath = fullfile(obj.Config.uipatdirs{1}, 'ea_merstate.mat');
            if exist(fpath, 'file')
                temp = load(fpath);
                obj = obj.clearData();
                obj.Config = temp.config;
                obj.DBSImplants = temp.dbsimplants;
                obj.Frame = temp.frame;
                obj.Markers = temp.markers;
                obj.Toggles = temp.toggles;
                obj = obj.loadDBSReconstruction();  % In case lead reconstruction changed since we last saved.
                obj = obj.setMERTrajectoriesToDefaults();
            else
                disp(['File not found: ' fpath]);
            end
        end
    end
    methods(Access = private)
        function obj = populateImplantCoords(obj, dbs_contacts, spc)
            % populateImplantCoords(dbs_contacts, spc)
            % dbs_contacts is a cell array, each with a 4x3 mat
            % spc is a string, either 'native' or 'mni'
            for side_ix = 1:length(dbs_contacts)
                if ~isempty(dbs_contacts{side_ix})
                    coords = dbs_contacts{side_ix};
                    
                    % Average 3-D distance through electrode contacts.
                    elec_diff = mean(diff(coords));
                    % unit vector through contacts from lowest (?) contact
                    elec_uv = elec_diff / norm(elec_diff);
                    % shift coordinates to top of lowest contact.
                    coords_top = bsxfun(@minus, coords,...
                        elec_uv * obj.Config.elspec.contact_length / 2);
                    
                    obj.DBSImplants(side_ix).coords.(spc) = coords;
                    obj.DBSImplants(side_ix).coords_top.(spc) = coords_top;
                    obj.DBSImplants(side_ix).elec_uv.(spc) = elec_uv;
                end
            end
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
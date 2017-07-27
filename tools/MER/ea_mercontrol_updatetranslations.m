function ea_mercontrol_updatetranslations(handles, side_str)
%ea_mercontrol_updatetranslations(handles, side_str)
% The MER translations need to be calculated infrequently...
% only when the frame orientation changes, or when
% the MER tract in which the DBS was inserted changes.
if ~exist('side_str', 'var')
    side_str = 'both';
end
[side_strs, side_ids, ~] = ea_detsidestr(side_str);

options = getappdata(handles.mercontrolfig, 'options');
merframe = getappdata(handles.mercontrolfig, 'merframe');
resultfig = getappdata(handles.mercontrolfig, 'resultfig');
merstruct = getappdata(resultfig, 'merstruct');
elspec = getappdata(resultfig, 'elspec');

% merstruct.tract_info.position assumes electrode vector goes toward -z.
transl_bengun = merstruct.offset * cat(1, merstruct.tract_info.position);

dbs_native = nan(length(side_ids), 3);
for sid = side_ids
    side_str = side_strs{sid};
    im_mm = merstruct.dbs_coords_mm.native{sid};  % Should be in scrf/native
    
    elec_diff = mean(diff(im_mm));  % Average 3-D distance through electrode contacts.
    elec_uv = elec_diff / norm(elec_diff); % unit vector through electrode contacts if tip was at origin.
    im_mm = bsxfun(@minus, im_mm, elec_uv * elspec.contact_length / 2); %shift coordinates to top of contact.
    
    side_transl = transl_bengun;
    if strcmpi(side_str, 'left')
        side_transl(:, 1) = -side_transl(:, 1);
    end
    
    % Offset by the MER tract the DBS electrode was implanted in.
    side_transl = bsxfun(@minus, side_transl, side_transl(merstruct.implant_idx(sid), :));
    
    % Apply yaw rotation.
    yaw_rad = merframe{sid}.yaw_rad;
    yaw_xform = [cos(yaw_rad) -sin(yaw_rad) 0; sin(yaw_rad) cos(yaw_rad) 0; 0 0 1];
    side_transl = side_transl * yaw_xform;
    
    if ~isequal(elec_uv, [0 0 1])
        % Find the transform that aligns [0 0 1] with the electrode unit
        % vector, then apply that transform to each of the MER offsets.
        U = align_vectors([0 0 1], elec_uv);
        side_transl = (U*side_transl')';
    end
    
    % Collect translations and DBS start points.
    merstruct.translations.native{sid} = side_transl;
    dbs_native(sid, :) = im_mm(1, :);
end

% Get transl_mni if needed
if ~options.native  % MNI
    dbs_mni = nan(size(dbs_native));
    for sid = side_ids
        im_mm = merstruct.dbs_coords_mm.mni{sid};
        elec_diff = mean(diff(im_mm));  % Average 3-D distance through electrode contacts.
        elec_uv = elec_diff / norm(elec_diff); % unit vector through electrode contacts if tip was at origin.
        im_mm = bsxfun(@minus, im_mm, elec_uv * elspec.contact_length / 2); %shift coordinates to top of contact.
        dbs_mni(sid, :) = im_mm(1, :);
    end
    
    transl_native = cat(1, merstruct.translations.native{side_ids});
    transl_native = reshape(transl_native, [length(merstruct.tract_info), length(side_ids), 3]);
    startpoint_native = bsxfun(@plus, transl_native, shiftdim(dbs_native, -1));
    
    ptdir = fullfile(options.root, options.patientname);
    options = ea_assignpretra(options);
    prenii_fname = fullfile(ptdir, options.prefs.prenii_unnormalized);
    nii = ea_load_nii(prenii_fname);
    
    [n_pos, n_sides, ~] = size(startpoint_native);
    c = reshape(startpoint_native, n_pos*n_sides, 3);
    c = [c, ones(size(c, 1), 1)]';
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
    c = ea_map_coords(c(1:3, :), prenii_fname, ...
        fullfile(ptdir, 'y_ea_inv_normparams.nii'), '');
    startpoint_mni = reshape(c', n_pos, n_sides, 3);
    clear n_sides n_pos c
    transl_mni = bsxfun(@minus, startpoint_mni, shiftdim(dbs_mni, -1));
    for sid = side_ids
        merstruct.translations.mni{sid} = squeeze(transl_mni(:, sid, :));
    end
end

setappdata(resultfig, 'merstruct', merstruct);

function RU = align_vectors(A, B)
cross_AB = cross(A,B);
ssc_cross_AB = [...
    0 -cross_AB(3) cross_AB(2);...
    cross_AB(3) 0 -cross_AB(1);...
    -cross_AB(2) cross_AB(1) 0];
RU = eye(3) + ssc_cross_AB + ssc_cross_AB^2*(1-dot(A,B))/(norm(cross_AB)^2);
function reco = ea_recalc_angles(reco)
% recalculate angles based on reco markers

if isfile(reco)
    recoFile = reco;
    load(recoFile, 'reco');
end

if ~isfield(reco, 'angles')
    reco.angles = ea_diode_emptyorientation(numel(reco.props));
end

spaces = {'mni', 'native'};
for s = 1:numel(spaces)
    space = spaces{s};
    if ~isfield(reco, space)
        continue;
    end

    m = reco.(space).markers;
    for i = 1:numel(m)
        if ~isempty(m(i).head) && ~isempty(m(i).tail)
            angles = ea_diode_trajectoryToYawPitchPolar(m(i).head, m(i).tail);
            reco.angles.(space)(i).yaw = angles.deg.yaw;
            reco.angles.(space)(i).pitch = angles.deg.pitch;
            reco.angles.(space)(i).polar = angles.deg.polar;
            clear angles
            angles = ea_diode_recomputeRoll(m(i).head, m(i).tail, m(i).y);            
            reco.angles.(space)(i).roll = angles.deg.roll;
            clear angles
            angles = ea_diode_recomputeOrientation(m(i).head, m(i).tail, m(i).y);
            reco.angles.(space)(i).orientation = angles.deg.orientation;
            clear angles
        end
    end
end

if exist('recoFile', 'var')
    save(recoFile, 'reco');
end

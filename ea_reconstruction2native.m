function ea_reconstruction2native(options)

if ~isfield(options.subj, 'norm')
    options = ea_getptopts(fileparts(fileparts(options.subj.recon.recon)));
end

load(options.subj.recon.recon, 'reco');

if ~exist('reco','var') % old format
    reco.mni.coords_mm = coords_mm;
    reco.mni.markers = markers;
    reco.mni.trajectory = trajectory;
    reco.props.elmodel = elmodel;
    reco.props.manually_corrected = manually_corrected;
end

json = loadjson(options.subj.norm.log.method);
if contains(json.method, 'SPM')
    normTemplate = [ea_space, 'TPM.nii'];
else
    spacedef = ea_getspacedef;
    normTemplate = [ea_space, spacedef.templates{1}, '.nii'];
end
nii = ea_load_nii(normTemplate);

if isfile(options.subj.brainshift.transform.scrf)
    usenative = 'scrf';
else
    usenative = 'native';
end

towarp = cell(0);
if ~isfield(options,'sides')
    options = ea_detsides(options);
end

for side = options.sides
    towarp{end+1} = reco.mni.coords_mm{side};
    towarp{end+1} = reco.mni.markers(side).head;
    towarp{end+1} = reco.mni.markers(side).tail;
    towarp{end+1} = reco.mni.trajectory{side};
end

towarp = cell2mat(towarp');
warpedcoord = ea_warpcoord(towarp,nii,options);

cnt = 1;
for side = options.sides
    offset = size(reco.mni.coords_mm{side},1);
    reco.(usenative).coords_mm{side} = warpedcoord(cnt:cnt+offset-1,:);
    cnt = cnt+offset;

    offset = size(reco.mni.markers(side).head,1);
    reco.(usenative).markers(side).head = warpedcoord(cnt:cnt+offset-1,:);
    cnt = cnt+offset;

    offset = size(reco.mni.markers(side).tail,1);
    reco.(usenative).markers(side).tail = warpedcoord(cnt:cnt+offset-1,:);
    cnt = cnt+offset;

    offset = size(reco.mni.trajectory{side},1);
    reco.(usenative).trajectory{side} = warpedcoord(cnt:cnt+offset-1,:);
    cnt = cnt+offset;

    [xunitv, yunitv] = ea_calcxy(reco.(usenative).markers(side).head, reco.(usenative).markers(side).tail);
    reco.(usenative).markers(side).x = reco.(usenative).markers(side).head+xunitv*(options.elspec.lead_diameter/2);
    reco.(usenative).markers(side).y = reco.(usenative).markers(side).head+yunitv*(options.elspec.lead_diameter/2);
end

% apply scrf to native matrix if available
if isfile(options.subj.brainshift.transform.scrf)
    d = load(options.subj.brainshift.transform.scrf);
    mat = inv(d.mat);
    reco.native = ea_applyscrfmat(mat, reco.scrf, options.sides);
else
    if isfield(reco,'scrf')
        reco = rmfield(reco,'scrf'); % delete subcortical transform if user apparently deleted the transform file.
    end
end

save(options.subj.recon.recon,'reco');


function c = ea_warpcoord(c,nii,options)
c = [c,ones(size(c,1),1)]';
c = nii(1).mat\c;

c = ea_map_coords(c(1:3,:), ...
    nii(1).fname, ...
    [options.subj.subjDir,filesep,'forwardTransform'],...
    options.subj.coreg.anat.preop.(options.subj.AnchorModality));
c = c';

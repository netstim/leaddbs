function ea_reconstruction2mni(options)

directory=[options.root,options.patientname,filesep];
load([directory,filesep,'ea_reconstruction.mat'],'reco');

%[~,template]=ea_whichnormmethod(directory);
options=ea_assignpretra(options);
nii=ea_load_nii([directory,options.prefs.prenii_unnormalized]);

if ~isfield(options,'elspec')
    options.elmodel=reco.props.elmodel;
    options=ea_resolve_elspec(options);
end

if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'],'file')
    usenative='scrf';
else
    usenative='native';
end

% apply native to scrf matrix if available
if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'],'file')
    d=load([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat']);
    reco.scrf=ea_applyscrfmat(d.mat,reco.native,options.sides);
elseif exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file') % legacy
    mat=ea_getscrfmat(options);
    save([directory,'scrf',filesep,'scrf_converted.mat'],'mat');
    reco.scrf=ea_applyscrfmat(mat,reco.native,options.sides);
else
    if isfield(reco,'scrf')
        reco=rmfield(reco,'scrf'); % delete subcortical transform if user apparently deleted the transform file.
    end
end

towarp=cell(0);
for side=options.sides
    towarp{end+1}=reco.(usenative).coords_mm{side};
    towarp{end+1}=reco.(usenative).markers(side).head;
    towarp{end+1}=reco.(usenative).markers(side).tail;
    towarp{end+1}=reco.(usenative).trajectory{side};
end
towarp=cell2mat(towarp');
warpedcoord=ea_warpcoord(towarp,nii,options);

cnt=1;
for side=options.sides
    offset=size(reco.(usenative).coords_mm{side},1);
    reco.mni.coords_mm{side}=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.(usenative).markers(side).head,1);
    reco.mni.markers(side).head=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.(usenative).markers(side).tail,1);
    reco.mni.markers(side).tail=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.(usenative).trajectory{side},1);
    reco.mni.trajectory{side}=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    if ~isempty(reco.mni.markers(side).head)
        % Enforce the rotation of x and y markers in MNI space. Use the
        % rotation of y marker calculated from native coordination.
        [~, yvec] = ea_calc_rotation(reco.native.markers(side).y,reco.native.markers(side).head);
        [xunitv, yunitv] = ea_calcxy(reco.mni.markers(side).head, reco.mni.markers(side).tail, yvec);
        reco.mni.markers(side).x = reco.mni.markers(side).head + xunitv*(options.elspec.lead_diameter/2);
        reco.mni.markers(side).y = reco.mni.markers(side).head + yunitv*(options.elspec.lead_diameter/2);
    else
        reco.mni.markers(side).x=[];
        reco.mni.markers(side).y=[];
    end
end

save([directory,filesep,'ea_reconstruction.mat'],'reco');


function c=ea_warpcoord(c,nii,options)
c=[c,ones(size(c,1),1)]';
% to template voxel space:
c=nii(1).mat\c;
try
    whichnormmethod=ea_whichnormmethod([options.root,options.patientname,filesep]);
    if ~ismember(whichnormmethod,ea_getantsnormfuns)
        V=spm_vol([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);
        if ~isequal(V.dim,nii.dim)
            ea_redo_inv([options.root,options.patientname,filesep],options);
        end
    end
end

c=ea_map_coords(c(1:3,:), ...
    [options.root,options.patientname,filesep,options.prefs.prenii_unnormalized], ...
    [options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'], ...
    '');
c=c';

function ea_reconstruction2native(options)

directory=[options.root,options.patientname,filesep];
load([directory,filesep,'ea_reconstruction.mat']);

if ~exist('reco','var') % old format
    reco.mni.coords_mm=coords_mm;
    reco.mni.markers=markers;
    reco.mni.trajectory=trajectory;
    reco.props.elmodel=elmodel;
    reco.props.manually_corrected=manually_corrected;
end

[whichnormmethod,template]=ea_whichnormmethod(directory);
nii=ea_load_nii(template);

if ~ismember(whichnormmethod,ea_getantsnormfuns)
    try
        ea_checkforwardinv(options,'forward')
    end
end

if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'],'file')
    usenative='scrf';
else
    usenative='native';
end

towarp=cell(0);
if ~isfield(options,'sides')
    options=ea_detsides(options);
end

for side=options.sides
    towarp{end+1}=reco.mni.coords_mm{side};
    towarp{end+1}=reco.mni.markers(side).head;
    towarp{end+1}=reco.mni.markers(side).tail;
    towarp{end+1}=reco.mni.trajectory{side};
end
towarp=cell2mat(towarp');
warpedcoord=ea_warpcoord(towarp,nii,options);

cnt=1;
for side=options.sides
    offset=size(reco.mni.coords_mm{side},1);
    reco.(usenative).coords_mm{side}=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.mni.markers(side).head,1);
    reco.(usenative).markers(side).head=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.mni.markers(side).tail,1);
    reco.(usenative).markers(side).tail=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.mni.trajectory{side},1);
    reco.(usenative).trajectory{side}=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    [xunitv, yunitv] = ea_calcxy(reco.(usenative).markers(side).head, reco.(usenative).markers(side).tail);
    reco.(usenative).markers(side).x = reco.(usenative).markers(side).head+xunitv*(options.elspec.lead_diameter/2);
    reco.(usenative).markers(side).y = reco.(usenative).markers(side).head+yunitv*(options.elspec.lead_diameter/2); % corresponding points in reality
end


% apply scrf to native matrix if available
if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat'],'file')
    d=load([options.root,options.patientname,filesep,'scrf',filesep,'scrf_converted.mat']);
    mat=inv(d.mat);
    reco.native=ea_applyscrfmat(mat,reco.scrf,options.sides);
elseif exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file') % legacy
    mat=ea_getscrfmat([options.root,options.patientname,filesep]);
    save([directory,'scrf',filesep,'scrf_converted.mat'],'mat');
    mat=inv(mat);
    reco.native=ea_applyscrfmat(mat,reco.scrf,options.sides);
else
    if isfield(reco,'scrf')
        reco=rmfield(reco,'scrf'); % delete subcortical transform if user apparently deleted the transform file.
    end
end

save([directory,filesep,'ea_reconstruction.mat'],'reco');


function c=ea_warpcoord(c,nii,options)

c=[c,ones(size(c,1),1)]';
% to template voxel space:
c=nii(1).mat\c;

[~,anats]=ea_assignpretra(options);
src=[options.root,options.patientname,filesep,anats{1}]; % assign src image as primary anat image here.
c=ea_map_coords(c(1:3,:), nii(1).fname, [options.root,options.patientname,filesep,'y_ea_normparams.nii'],...
    src);
c=c';

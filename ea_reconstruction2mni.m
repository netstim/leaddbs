function ea_reconstruction2mni(options)

directory=[options.root,options.patientname,filesep];
load([directory,filesep,'ea_reconstruction.mat']);

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
    mat=ea_getscrfmat([options.root,options.patientname,filesep]);
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
    towarp{end+1}=reco.(usenative).markers(side).x;
    towarp{end+1}=reco.(usenative).markers(side).y;
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

    offset=size(reco.(usenative).markers(side).x,1);
    reco.mni.markers(side).x=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.(usenative).markers(side).y,1);
    reco.mni.markers(side).y=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;

    offset=size(reco.(usenative).trajectory{side},1);
    reco.mni.trajectory{side}=warpedcoord(cnt:cnt+offset-1,:); cnt=cnt+offset;
    
    trajvector = diff([reco.mni.markers(side).head;reco.mni.markers(side).tail]);
    normtrajvector = trajvector/norm(trajvector);

    if ~isempty(reco.mni.markers(side).head)
        % calculates x and y using the warped marker.y, projecting it onto
        % the perpendicular plane to normtrajvector and then finding x via
        % the crossproduct.
        y =  reco.mni.markers(side).y - reco.mni.markers(side).head;
        y = y/norm(y);
        y = y - (dot(y,normtrajvector) / (norm(normtrajvector) ^2)) * normtrajvector;
        x = cross(y,normtrajvector);
        reco.mni.markers(side).x = reco.mni.markers(side).head + (x * (options.elspec.lead_diameter/2));
        reco.mni.markers(side).y = reco.mni.markers(side).head + (y * (options.elspec.lead_diameter/2));
        % if strcmp(options.elspec.orientation,'anterior')
        %     % new version which makes y point strictly anterior ~TD
        %     y = [0 normtrajvector(3) normtrajvector(2)];
        %     x = -cross(normtrajvector,y);
        %     y = (y/norm(y)) * (options.elspec.lead_diameter/2);
        %     x = (x/norm(x)) * (options.elspec.lead_diameter/2);
        %     reco.mni.markers(side).x=reco.mni.markers(side).head + x;
        %     reco.mni.markers(side).y=reco.mni.markers(side).head + y;
        % else
        %     orth = fliplr(null(normtrajvector))'*(options.elspec.lead_diameter/2);
        %     orth(1,:) = -orth(1,:);
        %     reco.mni.markers(side).x = reco.mni.markers(side).head+orth(:,1);
        %     reco.mni.markers(side).y = reco.mni.markers(side).head+orth(:,2); % corresponding points in reality
        % end
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

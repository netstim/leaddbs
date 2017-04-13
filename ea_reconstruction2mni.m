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

if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file')
    usenative='scrf';
else
    usenative='native';
end

% apply native to scrf matrix if available
if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file')
load([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'])
mat=ea_antsmat2mat(AffineTransform_float_3_3,fixed);
reco.scrf=ea_applyscrfmat(mat,reco.native);
else
    if isfield(reco,'scrf')
        reco=rmfield(reco,'scrf'); % delete subcortical transform if user apparently deleted the transform file.
    end
end   


for side=1:length(options.sides)
    reco.mni.coords_mm{side}=ea_warpcoord(reco.(usenative).coords_mm{side},nii,options);
    reco.mni.markers(side).head=ea_warpcoord(reco.(usenative).markers(side).head,nii,options);
    reco.mni.markers(side).tail=ea_warpcoord(reco.(usenative).markers(side).tail,nii,options);
    reco.mni.trajectory{side}=ea_warpcoord(reco.(usenative).trajectory{side},nii,options);
    
    normtrajvector{side}=diff([reco.mni.markers(side).head;...
        reco.mni.markers(side).tail])/...
        norm(diff([reco.mni.markers(side).head;...
        reco.mni.markers(side).tail]));    
    orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);
    
    reco.mni.markers(side).x=reco.mni.markers(side).head+orth(:,1)';
    reco.mni.markers(side).y=reco.mni.markers(side).head+orth(:,2)'; % corresponding points in reality    
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

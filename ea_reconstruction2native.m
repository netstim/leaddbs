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


if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file')
    usenative='scrf';
else
    usenative='native';
end







for side=1:length(options.sides)
    
    reco.(usenative).coords_mm{side}=ea_warpcoord(reco.mni.coords_mm{side},nii,options);
    reco.(usenative).markers(side).head=ea_warpcoord(reco.mni.markers(side).head,nii,options);
    reco.(usenative).markers(side).tail=ea_warpcoord(reco.mni.markers(side).tail,nii,options);
    reco.(usenative).trajectory{side}=ea_warpcoord(reco.mni.trajectory{side},nii,options);
    
    normtrajvector{side}=mean(diff(reco.(usenative).trajectory{side}))/norm(mean(diff(reco.(usenative).trajectory{side})));
    orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);
    
    reco.(usenative).markers(side).x=reco.(usenative).markers(side).head+orth(:,1)';
    reco.(usenative).markers(side).y=reco.(usenative).markers(side).head+orth(:,2)'; % corresponding points in reality
    
end

% apply scrf to native matrix if available
if exist([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'],'file')
    load([options.root,options.patientname,filesep,'scrf',filesep,'scrf.mat'])
    mat=inv(ea_antsmat2mat(AffineTransform_float_3_3,fixed));
    reco.native=ea_applyscrfmat(mat,reco.scrf);
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

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
    [whichnormmethod,tempfile]=ea_whichnormmethod(directory);
    nii=ea_load_nii(tempfile);
    if ~ismember(whichnormmethod,ea_getantsnormfuns)
        try
            ea_checkforwardinv(options,'forward')
        end
    end
        for side=options.sides
            
            
            reco.native.coords_mm{side}=ea_warpcoord(reco.mni.coords_mm{side},nii,options);
            reco.native.markers(side).head=ea_warpcoord(reco.mni.markers(side).head,nii,options);
            reco.native.markers(side).tail=ea_warpcoord(reco.mni.markers(side).tail,nii,options);
            reco.native.trajectory{side}=ea_warpcoord(reco.mni.trajectory{side},nii,options);
            
            normtrajvector{side}=mean(diff(reco.native.trajectory{side}))/norm(mean(diff(reco.native.trajectory{side})));
            orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);

            
            reco.native.markers(side).x=reco.native.markers(side).head+orth(:,1)';
            reco.native.markers(side).y=reco.native.markers(side).head+orth(:,2)'; % corresponding points in reality
            
        end
        
        
        save([directory,filesep,'ea_reconstruction.mat'],'reco');
        
    

function c=ea_warpcoord(c,nii,options)
c=[c,ones(size(c,1),1)]';
% to template voxel space:
c=nii(1).mat\c;

c=ea_map_coords(c(1:3,:), nii(1).fname, [options.root,options.patientname,filesep,'y_ea_normparams.nii'],...
     '');
 c=c';
 
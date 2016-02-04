function ea_reconstruction2mni(options)

directory=[options.root,options.patientname,filesep];
load([directory,filesep,'ea_reconstruction.mat']);

    %[~,tempfile]=ea_whichnormmethod(directory);
    nii=ea_load_nii([directory,options.prefs.prenii_unnormalized]);

    
        for side=options.sides
            reco.mni.coords_mm{side}=ea_warpcoord(reco.native.coords_mm{side},nii,options);
            reco.mni.markers(side).head=ea_warpcoord(reco.native.markers(side).head,nii,options);
            reco.mni.markers(side).tail=ea_warpcoord(reco.native.markers(side).tail,nii,options);
            reco.mni.trajectory{side}=ea_warpcoord(reco.native.trajectory{side},nii,options);
            
            normtrajvector{side}=mean(diff(reco.mni.trajectory{side}))/norm(mean(diff(reco.mni.trajectory{side})));
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
V=spm_vol([options.root,options.patientname,filesep,'y_ea_inv_normparams.nii']);
if ~isequal(V.dim,nii.dim)
   ea_redo_inv([options.root,options.patientname,filesep],options); 
end
end

c=ea_map_coords(c(1:3,:), '', [options.root,options.patientname,filesep,'y_ea_inv_normparams.nii'],...
     [options.root,options.patientname,filesep,options.prefs.prenii_unnormalized]);
c=c';
 
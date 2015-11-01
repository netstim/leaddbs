function ea_updatemodel(options)
% markers stay constant, electrode contacts (coords) may not if electrode
% model has changed.
% __________________________________________________________________________________
% Copyright (C) 2015 Charite University Medicine Berlin, Movement Disorders Unit
%
% Andreas Horn

load([options.root,options.patientname,filesep,'ea_reconstruction']);

    if ~exist('markers','var') % backward compatibility to old recon format
        for side=options.sides
            markers(side).head=coords_mm{side}(1,:);
            markers(side).tail=coords_mm{side}(4,:);            
            normtrajvector=(markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
            orth=null(normtrajvector)*(options.elspec.lead_diameter/2);
            markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
            markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
        end
    end
    
coords_mm=ea_resolvecoords(markers,options);
elmodel=options.elmodel;
if ~exist('manually_corrected','var');
    manually_corrected=0;
end
save([options.root,options.patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm','markers','elmodel','manually_corrected');

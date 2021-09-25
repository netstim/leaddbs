function options=ea_detsides(options)

if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
    load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
    if isfield(reco,'mni')
        sides=[];
        if isfield(reco,'native')
            for el=1:length(reco.native.markers)
                if ~isempty(reco.native.markers(el).head)
                    sides(end+1)=el;
                end
            end
        else
            for el=1:length(reco.mni.markers)
                if ~isempty(reco.mni.markers(el).head)
                    sides(end+1)=el;
                end
            end
        end
        options.sides=sides;
    end
end

function options=ea_detsides(options)


if exist([options.root,options.patientname,filesep,'ea_reconstruction.mat'],'file')
    load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
    sides=[];
    for el=1:length(reco.native.coords_mm)
        if ~isempty(reco.native.coords_mm(el))
            sides(end+1)=el;
        end
    end
    options.sides=sides;

end
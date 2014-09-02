function ea_write(options)

if options.elspec.numel>4 % here, electrodes with more than 4 contacts are extrapolated from 4fold data.
    try
    try
        load([options.root,options.patientname,filesep,'ea_reconstruction']);
    catch
        error('No reconstruction information found. Please run reconstruction first.');
    end
    for side=1:length(coords_mm)
    if size(coords_mm{side},1)~=options.elspec.numel
        
        coords_mm{side}=ea_extrapol_coords(coords_mm{side},options);
        
        ea_write_fiducials(coords_mm,fullfile(options.root,options.patientname,'ea_coords.fcsv'),options);
        save([options.root,options.patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm');
    end
    end
    end 
end




% Slice 2D Visualization
if options.d2.write
    
    % Prior Results are loaded here inside the function (this way, function
    % can be called just by giving the patient directory.
    
    cuts=ea_writeplanes(options);
    
end



% Render 3D Visualization
if options.d3.write
    
    % Prior Results are loaded here inside the function (this way, function
    % can be called just by giving the patient directory.
    
    resultfig=ea_render_view(options);
    
    % save scene as matlab figure
    try % if path is not defined, don't save.
    saveas(resultfig,[options.root,options.patientname,filesep,'LEAD_scene.fig']);
    end
    %figure2xhtml([options.root,options.patientname,filesep,'eAuto_scene'],resultfig);
    if options.d3.autoserver
       ea_export_server([],[],options); 
    end
    
end

%% check traject sanity

for side=options.sides
    try
        trajectissane=ea_checktrajectsanity(trajvector{side});
        if ~trajectissane
            disp(['Trajectory of side ',num2str(side),' seems not to have been correctly reconstructed. Check manually.']);
        end
        
    end
end


try
    if isnan(results)
        clear results
    end
    results.coords_mm=coords_mm;
    
    results.realcoords=realcoords;
    
    
    
    for electrode=1:length(coords_mm)
        
        results.distances(electrode)=pdist([coords_mm(electrode,:);realcoords(electrode,:)]);
        
    end
    results.fit=nanmean(results.distances);
    
end


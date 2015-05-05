function ea_autocoord(options)
% This function is the main function of LEAD-DBS. It will generate a 8x3
% vector of coordinates. Rows 1-4 will be the right electrodes, rows 5-8
% the left ones. trajectory{1} will be the right trajectory, trajectory{2} the
% left one.
% For each hemisphere of the brain, this function will call the
% reconstruction routine ea_autocoord_side and lateron call functions for
% manual correction of the results, render and slice views.
% __________________________________________________________________________________
% Copyright (C) 2014 Charite University Medicine Berlin, Movement Disorders Unit
% Andreas Horn




%%
% get accurate electrode specifications and save it in options.
options=ea_resolve_elspec(options);
patientname=options.patientname;
options.prefs=ea_prefs(patientname);




if options.dicomimp % do DICOM-Import.
    
    ea_dicom_import(options);
    
    
end

if ~strcmp(options.patientname,'No Patient Selected') % only 3D-rendering viewer can be opened if no patient is selected.
    



if options.modality==2 % CT support
    options.prefs.tranii=options.prefs.ctnii;
    options.prefs.tranii_unnormalized=options.prefs.rawctnii_unnormalized;
end
results=nan;

if options.coregct.do
    eval([options.coregct.method,'(options)']); % triggers the coregct function and passes the options struct to it.
    try load([options.root,options.patientname,filesep,'ea_coregctmethod_applied']); end
    if exist('coregct_method_applied','var')
        try
            coregct_method_applied{end+1}=options.coregct.method;
        catch
            clear coregct_method_applied
            coregct_method_applied{1}=options.coregct.method;
        end
    else
        coregct_method_applied{1}=options.coregct.method;
    end
    coregct_method_applied=options.coregct.method;
    save([options.root,options.patientname,filesep,'ea_coregctmethod_applied'],'coregct_method_applied');
end

if options.coregctcheck
    
    
    % export "control" niftis with wireframe of normal anatomy..
    
    ea_show_ctcoregistration(options);
    
end


if options.normalize.do
    
    eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.
    try load([options.root,options.patientname,filesep,'ea_normmethod_applied']); end
    if exist('norm_method_applied','var')
        try
            norm_method_applied{end+1}=options.normalize.method;
        catch
            clear norm_method_applied
            norm_method_applied{1}=options.normalize.method;
        end
    else
        norm_method_applied{1}=options.normalize.method;
    end
    save([options.root,options.patientname,filesep,'ea_normmethod_applied'],'norm_method_applied');

end




if options.dolc

    ea_perform_lc(options);
    
end


if options.atl.normalize % normalize patient-specific atlas-set.
    ea_normalize_ptspecific_atl(options)
end
    
% Prepare MR-images

if ~exist([options.root,options.prefs.patientdir,filesep,options.prefs.tranii],'file')
    try
        copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.tranii,'.gz'],[options.root,options.prefs.patientdir,filesep,'c',options.prefs.tranii,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.tranii,'.gz']);
    catch
        system(['gunzip ',options.root,options.prefs.patientdir,filesep,options.prefs.tranii,'.gz']);
        
    end
    try
        movefile([options.root,options.prefs.patientdir,filesep,'c',options.prefs.tranii,'.gz'],[options.root,options.prefs.patientdir,filesep,options.prefs.tranii,'.gz']);
    end
end

if ~exist([options.root,options.prefs.patientdir,filesep,options.prefs.cornii],'file')
    try
        copyfile([options.root,options.prefs.patientdir,filesep,options.prefs.cornii,'.gz'],[options.root,options.prefs.patientdir,filesep,'c',options.prefs.cornii,'.gz']);
    end
    try
        gunzip([options.root,options.prefs.patientdir,filesep,options.prefs.cornii,'.gz']);
    catch
        system(['gunzip ',options.root,options.prefs.patientdir,filesep,options.prefs.cornii,'.gz']);
    end
    try
        movefile([options.root,options.prefs.patientdir,filesep,'c',options.prefs.cornii,'.gz'],[options.root,options.prefs.patientdir,filesep,options.prefs.cornii,'.gz']);
    end
    
end


if options.normalize.check
    
    
    % export "control" niftis with wireframe of normal anatomy..
    
    ea_show_normalization(options);
    
end


if options.atl.genpt % generate patient specific atlas set
    ea_ptspecific_atl(options);
    
end

if options.doreconstruction
    ea_checkfiles(options);
    
    
    for side=options.sides
        %try
        % call main routine reconstructing trajectory for one side.
        [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(patientname,options,side);
        
        
        
        %% refit electrodes starting from first electrode (this is redundant at
        %% this point).
        
        
        
        
        coords_mm{side} = ea_map_coords(coords', [options.root,options.prefs.patientdir,filesep,options.prefs.tranii])';
        
        
        [~,distmm]=ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),[options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
        
        comp = ea_map_coords([0,0,0;trajvector{side}]', [options.root,options.prefs.patientdir,filesep,options.prefs.tranii])'; % (XYZ_mm unaltered)
        
        trajvector{side}=diff(comp);
        
        
        
        normtrajvector{side}=trajvector{side}./norm(trajvector{side});
        
        for electrode=2:4
            
            coords_mm{side}(electrode,:)=coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
            
        end
        
        
        
        
        

    end
    
    
    
%     if length(coords_mm)==4 % only one side was processed.
%         if options.sides==1
%             coords_mm=[nan(4,3);coords_mm];
%         elseif options.sides==2
%             coords_mm=[coords_mm;nan(4,3)];
%         end
%     end
    
    
    
    
    
    try
        %realcoords=load([options.root,patientname,filesep,'L.csv']);
        %realcoords=realcoords(:,1:3);
        
        realcoords=ea_read_fiducials([options.root,patientname,filesep],'L');
    catch
        ea_showdis('No manual survey available.',options.verbose);
        realcoords=[];
    end
    
    % transform trajectory to mm space:
    for side=options.sides
        
        try
            if ~isempty(trajectory{side})
                trajectory{side}=ea_map_coords(trajectory{side}', [options.root,patientname,filesep,options.prefs.tranii])';
            end
            
        end
    end
    
    
    % save reconstruction results
try    ea_write_fiducials(coords_mm,fullfile(options.root,patientname,'ea_coords.fcsv'),options); end
    elmodel=options.elmodel;
    save([options.root,patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm','elmodel');
    
end


if options.manualheightcorrection
    
    % load reconstruction results
    try
        load([options.root,patientname,filesep,'ea_reconstruction']);
    catch
        ea_error('No reconstruction information found. Please run reconstruction first.');
    end
    mcfig=figure('name',[patientname,': Manual Height Correction'],'numbertitle','off');
    ea_manualcorrection(mcfig,coords_mm,trajectory,patientname,options);
    
    
   
else
    ea_write(options)
end

else
    ea_write(options)
end










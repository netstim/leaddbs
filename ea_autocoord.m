function ea_autocoord(options)
% This function is the main function of LEAD-DBS. It will generate a
% vector of coordinates.
% Trajectory{1} will be the right trajectory, trajectory{2} the
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




% check connectome-mapper tags
if isfield(options,'lcm')
    if options.lcm.func.do
        ea_lcm_func(options);
    end
    if options.lcm.struc.do
        ea_lcm_struc(options);
    end
end



if ~strcmp(options.patientname,'No Patient Selected') % only 3D-rendering viewer can be opened if no patient is selected.
    try    ea_compat_patfolder(options); end
    
    [options,presentfiles]=ea_assignpretra(options);
    
    
    if ~isempty(presentfiles)
       if ~exist(ea_niigz([options.root,options.patientname,filesep,'grid.nii']),'file')
           ea_gengrid(options);
       end
    end
    
    try ea_resliceanat(options); end
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
        ea_dumpnormmethod(options,options.normalize.method);
        eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.
    end
    
    
    
    
    if options.dolc % perform lead connectome subroutine..
        ea_perform_lc(options);
    end
    
    
    
    if options.atl.genpt % generate patient specific atlas set
        ea_ptspecific_atl(options);
        
    end
    
    if options.atl.normalize % normalize patient-specific atlas-set.
        ea_normalize_ptspecific_atl(options)
    end
    
    
    
    if options.normalize.check
        
        
        % export "control" niftis with wireframe of normal anatomy..
        
        ea_show_normalization(options);
        
    end
    
    
    
    if options.doreconstruction
        ea_checkfiles(options);
        
        
        for side=1:length(options.sides)
            %try
            % call main routine reconstructing trajectory for one side.
            [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(patientname,options,options.sides(side));
            
            
            
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
            markers(side).head=coords_mm{side}(1,:);
            markers(side).tail=coords_mm{side}(4,:);
            
            
            orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);
            
            markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
            markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
            
            
            
            coords_mm=ea_resolvecoords(markers,options);
            
        end
        
        
        
        %     if length(coords_mm)==4 % only one side was processed.
        %         if options.sides==1
        %             coords_mm=[nan(4,3);coords_mm];
        %         elseif options.sides==2
        %             coords_mm=[coords_mm;nan(4,3)];
        %         end
        %     end
        
        
        
        
        
        
        % transform trajectory to mm space:
        for side=1:length(options.sides)
            
            try
                if ~isempty(trajectory{side})
                    trajectory{side}=ea_map_coords(trajectory{side}', [options.root,patientname,filesep,options.prefs.tranii])';
                end
                
            end
        end
        
        
        % save reconstruction results
        elmodel=options.elmodel;
        ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
        
    end
    
    
    if options.manualheightcorrection
        % load reconstruction results
        %     try
        %         [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
        %     catch
        %         ea_error([patientname,': No reconstruction information found. Please run reconstruction first.']);
        %     end
        %     ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
        
        
        mcfig=figure('name',[patientname,': Manual Height Correction'],'numbertitle','off');
        warning('off');
        try
            ea_maximize(mcfig);
        end
        
        ea_manualreconstruction(mcfig,patientname,options);
        
    else
        ea_write(options)
    end
    

    
else
    ea_write(options)
end










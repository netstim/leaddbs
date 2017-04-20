function  [coords_mm,trajectory,markers,elmodel,manually_corrected,coords_acpc]=ea_load_reconstruction(varargin)

options=varargin{1};

try
    load([options.root,options.patientname,filesep,'ea_reconstruction']);
catch
	try
        coords_mm=ea_read_fiducials([options.root,options.patientname,filesep,'ea_coords.fcsv'],options);
    catch
        warning(['Please localize electrodes of ',options.patientname,' first.']);
	end
end

if exist('reco','var')

	if ~isfield(reco,'native') && isfield(reco,'mni') && options.native
        ea_reconstruction2native(options);
        load([options.root,options.patientname,filesep,'ea_reconstruction']);
	elseif isfield(reco,'native') && ~isfield(reco,'mni') && ~options.native
        ea_reconstruction2mni(options);
        load([options.root,options.patientname,filesep,'ea_reconstruction']);
    end
    
    if options.native && isfield(options,'loadrecoforviz') % if loading reco for visualization, should return scrf.
        if isfield(reco,'scrf')
            coords_mm=reco.scrf.coords_mm;
            trajectory=reco.scrf.trajectory;
            markers=reco.scrf.markers;
        else
            coords_mm=reco.native.coords_mm;
            trajectory=reco.native.trajectory;
            markers=reco.native.markers;
        end
    elseif options.native && ~isfield(options,'loadrecoforviz')
        coords_mm=reco.native.coords_mm;
        trajectory=reco.native.trajectory;
        markers=reco.native.markers;
    else
        coords_mm=reco.mni.coords_mm;
        trajectory=reco.mni.trajectory;
        markers=reco.mni.markers;
    end
    
    try
        coords_acpc=reco.acpc.coords_mm;
    catch
        coords_acpc=nan;
    end
    
    manually_corrected=reco.props.manually_corrected;
    elmodel=reco.props.elmodel;

else % legacy format
    
	if ~exist('markers','var') % backward compatibility to old recon format
        for side=1:length(options.sides)
            markers(side).head=coords_mm{side}(1,:);
            markers(side).tail=coords_mm{side}(4,:);
            normtrajvector=(markers(side).tail-markers(side).head)./norm(markers(side).tail-markers(side).head);
            orth=null(normtrajvector)*(options.elspec.lead_diameter/2);
            markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
            markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
        end
        
        elmodel=options.elmodel;
        % this is still legacy format but has markers variable in it now.
        save([options.root,options.patientname,filesep,'ea_reconstruction'],'trajectory','coords_mm','markers','elmodel');
    end

    try
        load([options.root,options.patientname,filesep,'ea_reconstruction.mat']);
    catch % generate trajectory from coordinates.
        trajectory{1}=ea_fit_line(coords_mm(1:4,:));
        trajectory{2}=ea_fit_line(coords_mm(options.elspec.numel+1:options.elspec.numel+4,:));
    end

    % we have all variables now but need to put them into reco format and
    % they are in MNI.
    options.native=0;
    options.hybridsave=1;
    
    try
        ea_save_reconstruction(coords_mm,trajectory,markers,'Medtronic 3389',0,options);
        options.native=1;
        [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
    end
end
  
if ~exist('manually_corrected','var')
	manually_corrected=0;
end

if ~exist('elmodel','var')
	elmodel=[];
end

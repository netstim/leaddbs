function ea_write(options)

try
    ea_updatemodel(options);
end

% check if sides is specified correctly for visualization
if isfield(options, 'subj')
    options=ea_detsides(options);
end

if options.d2.write || options.d3.write
   if strcmp(options.atlasset,'Segment patient anatomy')
      ea_subcorticalsegmentation(options);

      if ~options.native % normalize 2 MNI space
          ea_normsubcorticalsegm(options);
      end
      options.atlasset=options.patientname;
      options.atl.pt=1;
      options.atl.can=0;

   elseif startsWith(options.atlasset, 'Local atlas: ')
       options.atlasset = erase(options.atlasset, 'Local atlas: ');

   elseif strcmp(options.atlasset,'Use none')
       % do nothing
   end
end

% Slice 2D Visualization
if options.d2.write
    % Prior Results are loaded here inside the function (this way, function
    % can be called just by giving the patient directory.
    ea_writeplanes(options);
end

% Render 3D Visualization
if options.d3.write
    % Prior Results are loaded here inside the function (this way, function
    % can be called just by giving the patient directory.
    resultfig=ea_elvis(options);
    set(0,'CurrentFigure',resultfig);

    % zoom out contacts when patient is selected and electroddes are reconstructed
    if ~strcmp(options.patientname, 'No Patient Selected')
        try
            % calculate the center of the contacts
            coords = ea_load_reconstruction(options);
            center = mean(cell2mat(cellfun(@mean, coords, 'Uniformoutput',0)'), 1);
            % zoom on the center of the contacts
            ea_zoomcenter(resultfig.CurrentAxes, center, 3);
        catch
            zoom(3);
        end
    end

    % save scene as matlab figure
    % try % if path is not defined, don't save.
    %     savefig(resultfig,[options.root,options.patientname,filesep,'LEAD_scene.fig'],'compact');
    % end
    % figure2xhtml([options.root,options.patientname,filesep,'eAuto_scene'],resultfig);

    if options.d3.autoserver
        ea_export_server([],[],options);
        close(resultfig);
    end
end

%% check traject sanity

for iside=1:length(options.sides)
    side=options.sides(iside);
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
        results.distances(electrode)=ea_pdist([coords_mm(electrode,:);realcoords(electrode,:)]);
    end
    results.fit=ea_nanmean(results.distances);
end

% chirp on completed task.
ea_chirp(options);


function y = ea_nanmean(varargin)
if nargin==2
    x=varargin{1};
    dim=varargin{2};
elseif nargin==1
    x=varargin{1};
    dim=1;
end

N = sum(~isnan(x), dim);
y = nansum(x, dim) ./ N;

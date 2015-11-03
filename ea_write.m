function ea_write(options)
try
ea_updatemodel(options);
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
    
    resultfig=ea_elvis(options);
    
    % save scene as matlab figure
    try % if path is not defined, don't save.
    saveas(resultfig,[options.root,options.patientname,filesep,'LEAD_scene.fig']);
    end
    %figure2xhtml([options.root,options.patientname,filesep,'eAuto_scene'],resultfig);
    if options.d3.autoserver
       ea_export_server([],[],options);
       close(resultfig);
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
    results.fit=ea_nanmean(results.distances);
    
end


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

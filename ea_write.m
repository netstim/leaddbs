function ea_write(options)


if options.scrf
       ea_subcorticalrefine(options); 
end


try
    ea_updatemodel(options);
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
      

   elseif length(options.atlasset)>=13 && strcmp(options.atlasset(1:13),'Local atlas: ')
       options.atlasset=options.atlasset(14:end);
   
   elseif strcmp(options.atlasset,'Use none')
       % do nothing
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
    
    resultfig=ea_elvis(options);
    
    % save scene as matlab figure
%     try % if path is not defined, don't save.
%         savefig(resultfig,[options.root,options.patientname,filesep,'LEAD_scene.fig'],'compact');
%     end
    %figure2xhtml([options.root,options.patientname,filesep,'eAuto_scene'],resultfig);
    if options.d3.autoserver
       ea_export_server([],[],options);
       close(resultfig);
    end
end

%% check traject sanity

for side=1:length(options.sides)
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

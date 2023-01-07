function ea_screenshots(uipatdir,target)
% Wrapper of ea_screenshot used for exporting ZIP

options=ea_getptopts(uipatdir);
options.fiberthresh=1;
options.d3.verbose='on';
options.native=0;
options=ea_detsides(options);

[~,~,~,elmodel]=ea_load_reconstruction(options);
options.leadprod='dbs';
options.elmodel=elmodel;
options=ea_resolve_elspec(options);
options=ea_detsides(options);
options.d3.elrendering=1;
options.d3.hlactivecontacts=0;
options.d3.writeatlases=1;
options.d3.showisovolume=0;
viewsets=load([ea_getearoot,'helpers',filesep,'export',filesep,'ea_exportviews']);
if target=="VIM"
    viewsets.VIM.atlas='Essential Tremor Hypointensity (Neudorfer 2022)';
    viewsets.VIM.views(1).structures={'VIM','DRTT','Hypointensity'};
    viewsets.VIM.views(2).structures={'VIM','DRTT','Hypointensity'};
    viewsets.VIM.views(3).structures={'VIM','DRTT','Hypointensity'};
    viewsets.VIM.views(4).structures={'VIM','DRTT','Hypointensity'};
elseif target=="GPi"
    viewsets.GPi.views(1).structures={'GPi'};
    viewsets.GPi.views(2).structures={'GPi'};
    viewsets.GPi.views(3).structures={'GPi'};
    viewsets.GPi.views(4).structures={'GPi'};
elseif target=="STN"
    viewsets.STN.views(1).structures={'STN'};
    viewsets.STN.views(2).structures={'STN'};
    viewsets.STN.views(3).structures={'STN'};
    viewsets.STN.views(4).structures={'STN'};
end
options.atlasset=viewsets.(target).atlas;
options.sidecolor=1;
options.writeoutstats=0;
options.patientname = options.subj.subjId;
resultfig=ea_elvis(options);

if ~exist(fullfile(options.subj.subjDir, 'export', 'views'),'dir')
    mkdir(fullfile(options.subj.subjDir, 'export', 'views'));
end

views=viewsets.(target).views;
for hh=1:length(views)
    views(hh).v.camva=0.83; %changing zoom
end 


for view=1:length(views)
    set(0,'CurrentFigure',resultfig);
    if target=="VIM"
        ea_keepatlaslabels(views(view).structures{1,2:3});
        ea_meshatlaslabels(views(view).structures{1,1});
    else 
        ea_keepatlaslabels(views(view).structures{:});
    end
    set(0,'CurrentFigure',resultfig);
    ea_setplanes(views(view).planes.x,views(view).planes.y,views(view).planes.z);
    set(0,'CurrentFigure',resultfig);
    ea_view(views(view).v,resultfig);
    set(0,'CurrentFigure',resultfig);
    png_path=fullfile(options.subj.subjDir, 'export', 'views', ['view_',sprintf('%03.0f',view),'.png']);
    ea_screenshot(png_path,'ld', resultfig);
    J = imresize(imread(png_path),0.4,"triangle");
    imwrite(J,png_path);
end
close(resultfig);

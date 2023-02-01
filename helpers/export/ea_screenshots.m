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
    viewsets.VIM.views(1).structures={'VIM'}; %,'DRTT','Hypointensity'};
    viewsets.VIM.views(2).structures={'VIM'}; %^,'DRTT','Hypointensity'};
    viewsets.VIM.views(3).structures={'VIM'}; %,'DRTT','Hypointensity'};
    viewsets.VIM.views(4).structures={'VIM'}; %,'DRTT','Hypointensity'};
elseif target=="GPi"
    viewsets.GPi.views(1).structures={'GPi'};
    viewsets.GPi.views(2).structures={'GPi'};
    viewsets.GPi.views(3).structures={'GPi'};
    viewsets.GPi.views(4).structures={'GPi'};
    viewsets.STN.views(5).structures={'STN'};
    viewsets.STN.views(6).structures={'STN'};
elseif target=="STN"
    viewsets.STN.views(1).structures={'STN'};
    viewsets.STN.views(2).structures={'STN'};
    viewsets.STN.views(3).structures={'STN'};
    viewsets.STN.views(4).structures={'STN'};
    viewsets.STN.views(5).structures={'STN'};
    viewsets.STN.views(6).structures={'STN'};
end
options.atlasset=viewsets.(target).atlas;
options.sidecolor=1;
options.writeoutstats=0;
resultfig=ea_elvis(options);

if ~exist([options.root,options.patientname,filesep,'export',filesep,'views'],'dir')
    mkdir([options.root,options.patientname,filesep,'export',filesep,'views']);
end

views=viewsets.(target).views;

% Adjusting Dorsal View
views(4).v.az=180;
views(4).v.el=74.08;
views(4).v.camva=1.0197;
views(4).v.camup=[0 0 -1];
views(4).v.camtarget=[-0.2429 -8.1055 -8.0924];
views(4).v.campos=[-0.1791 967.6553 1.8607e+03];
views(4).planes.y=-100;

% Adding Left Oblique View
views(5).v.az=-120;
views(5).v.el=36;
views(5).v.camva=1.0197;
views(5).v.camup=[0 0 1];%[0.2337 -0.3866 0.8921];
views(5).v.camtarget=[-0.2429 -8.1055 -8.0924];
views(5).v.camproj = 'orthographic';
views(5).v.campos=[-1.2278e+03 894.6170 847.5294];
views(5).planes=views(4).planes;
views(5).planes.x=200;

%Adding Right Oblique View
views(6).v.az=60;
views(6).v.el=36;
views(6).v.camva=1.0197;
views(6).v.camup=[0 0 1];
views(6).v.camtarget=[-0.2429 -8.1055 -8.0924];
views(6).v.campos=[1.2278e+03 894.6170 847.5294];
views(6).v.camproj = 'orthographic';
views(6).planes=views(4).planes;
views(6).planes.x=-100;

for view=1:length(views)
    set(0,'CurrentFigure',resultfig);
    ea_keepatlaslabels(views(view).structures{:});
    %ea_keepatlaslabels(target);
    set(0,'CurrentFigure',resultfig);
    ea_setplanes(views(view).planes.x,views(view).planes.y,views(view).planes.z);
    set(0,'CurrentFigure',resultfig);
    ea_view(views(view).v,resultfig);
    set(0,'CurrentFigure',resultfig);
    png_path=fullfile(options.subj.subjDir, 'export', 'views', ['view_',sprintf('%03.0f',view),'.png']);
    ea_screenshot(png_path,'ld', resultfig);
end
close(resultfig);

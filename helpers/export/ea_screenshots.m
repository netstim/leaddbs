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
options.atlasset=viewsets.(target).atlas;
options.sidecolor=1;
options.writeoutstats=0;
resultfig=ea_elvis(options);

if ~exist([options.root,options.patientname,filesep,'export',filesep,'views'],'dir')
    mkdir([options.root,options.patientname,filesep,'export',filesep,'views']);
end

views=viewsets.(target).views;
for view=1:length(views)
    set(0,'CurrentFigure',resultfig);
    ea_keepatlaslabels( views(view).structures{:});
    set(0,'CurrentFigure',resultfig);
    ea_setplanes(views(view).planes.x,views(view).planes.y,views(view).planes.z);
    set(0,'CurrentFigure',resultfig);
    ea_view(views(view).v,resultfig);
    set(0,'CurrentFigure',resultfig);
    ea_screenshot([options.root,options.patientname,filesep,'export',filesep,'views',filesep,'view_',sprintf('%03.0f',view),'.png'],'ld', resultfig);
end
close(resultfig);

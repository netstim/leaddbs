function ea_openatlascontrol(hobj,ev,atlases,resultfig,options)
atlassurfs=getappdata(resultfig,'atlassurfs');
colorbuttons=getappdata(resultfig,'colorbuttons');
aswin=ea_atlasselect(colorbuttons,atlassurfs,atlases,options,resultfig);
set(aswin,'visible',options.d3.verbose);

setappdata(resultfig,'aswin',aswin);

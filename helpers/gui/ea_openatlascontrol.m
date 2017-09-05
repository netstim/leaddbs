function ea_openatlascontrol(hobj, ev, atlases, resultfig, options)
atlassurfs = getappdata(resultfig, 'atlassurfs');
colorbuttons = getappdata(resultfig, 'colorbuttons');
labelbutton = getappdata(resultfig,'labelbutton');
atlaslabels = getappdata(resultfig,'atlaslabels');
aswin = ea_atlasselect(colorbuttons, atlassurfs, atlases, options, resultfig, labelbutton, atlaslabels);
set(aswin, 'visible', options.d3.verbose);

setappdata(resultfig, 'aswin', aswin);

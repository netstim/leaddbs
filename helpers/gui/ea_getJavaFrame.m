function jFrame = ea_getJavaFrame(hFig)
% Get JavaFrame from hFig

warnState1 = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warnState2 = warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');
jFrame = get(hFig,'JavaFrame');
warning(warnState2);
warning(warnState1);

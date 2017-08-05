function [jContentPane, jFrame] = getJavaFrame(hFig)
% Get javaContentPane and javaFrame from hFig

wrn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(hFig,'JavaFrame');
jContentPane = jFrame.fHG2Client().getContentPane();
warning(wrn);

function [jContentPane, jFrame] = getJavaContentPane(hFig)
% Get javaContentPane and javaFrame from hFig

wrn = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(hFig,'JavaFrame');
jContentPane = jFrame.fHG2Client(1).getContentPane();
warning(wrn);

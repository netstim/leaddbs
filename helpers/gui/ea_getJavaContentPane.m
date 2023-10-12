function jContentPane = ea_getJavaContentPane(hFig)
% Get JavaContentPane from hFig

jFrame = ea_getJavaFrame(hFig);
jContentPane = jFrame.fHG2Client(1).getContentPane();

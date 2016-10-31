% getMonitorPositions.m
%
%      usage: monitorPositions = getMonitorPositions
%         by: julien besle (provided by The Mathworks)
%       date: 12/01/2011
%        $Id$
%    purpose: get correct MonitorPositions property or root object
%             using a Java workaround to replace get(0,'monitorPositions') 
%             which bug on Macs from version 7.10b on
%             Note that it returns (almost) the correct positions
%             consistent with  the 'position' property of figures)
%             There is no need to call correctMonitorPosition
%             althought there is a slight problem with the vertical position coordinates that needs to be fixed

function monitorPositions = getMonitorPositions

ge = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment();
gs = ge.getScreenDevices();

monitorPositions=zeros(length(gs),4);
for iter=1:length(gs)
    gb = gs(iter).getDefaultConfiguration().getBounds();
    monitorPositions(iter,:) = [gb.getX(), gb.getY(), gb.getWidth(), gb.getHeight()];
end


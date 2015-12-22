function h = mesocreateStopButton
hmain = gcf;
pos = get(hmain,'Position');

h = figure(539375677);
%set(h,'WindowStyle','modal');

%set(handle,'title','FiberGT Progress Control');
set(h,'NumberTitle','off');
set(h,'resize','off');
set(h,'MenuBar','none')
set(h,'ToolBar','none')
set(h,'Visible','off')
set(h,'Tag','stopbuttonpresent');
set(h,'Name','Progress Control');
set(h,'Position',[pos(1)+pos(3)/2 pos(2)+pos(4)/2 200 50]);
bh1 = uicontrol(h,'Position',[10 10 180 30],...
                'String','Stop Tracking',...
                'Callback',@closefig);
drawnow;
function closefig(hObject,eventdata)
close(539375677);



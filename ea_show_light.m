function [resultfig,lightbulb]=ea_show_light(resultfig)
        set(0,'CurrentFigure',resultfig); 

ceiling_lamp=light('Position',[0 0 100]); % not modifiable, infinite light.
right_lamp=light('Position',[100 0 0]); % not modifiable, infinite light.
left_lamp=light('Position',[-100 0 0]); % not modifiable, infinite light.

lightobj=light('Position',[30 30 30],'Style','local');
lightbulb=plot3(30, 30, 30,'o','MarkerSize',20,'MarkerFaceColor','y','MarkerEdgeColor','k');
setappdata(resultfig,'lightbulb',lightbulb);
setappdata(resultfig,'lightobj',lightobj);

setappdata(resultfig,'right_lamp',right_lamp);
setappdata(resultfig,'left_lamp',left_lamp);
setappdata(resultfig,'ceiling_lamp',ceiling_lamp);

ea_draggable(lightbulb,@movelight);


function movelight(lightbulb)
xx=get(lightbulb,'xdata');
yy=get(lightbulb,'ydata');
zz=get(lightbulb,'zdata');
lightobj=getappdata(gcf,'lightobj');
set(lightobj,'Position',[xx yy zz]);

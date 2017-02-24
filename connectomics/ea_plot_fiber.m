function [h,fv]=ea_plot_fiber(thisfib,numpol,show,options)
% small function selecting to either draw polygonized tube or thin line.

switch options.prefs.d3.fiberstyle
    case 'tube'
        [h,fv]=ea_plot3t(thisfib(1,:),thisfib(2,:),thisfib(3,:),options.prefs.d3.fiberdiameter,thisfib(4:end,:),numpol,show);
    case 'line'
        h=surface([thisfib(1,:);thisfib(1,:)],...
            [thisfib(2,:);thisfib(2,:)],...
            [thisfib(3,:);thisfib(3,:)],...
            [thisfib(4,:);thisfib(4,:)],'facecol','no','edgecol','interp','linew',options.prefs.d3.fiberdiameter);
        fv=nan;
    otherwise
        ea_error('Please set ea_prefs.d3.fiberstyle to either tube or line');
end
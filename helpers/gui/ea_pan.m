function ea_pan(h,~,cmd)

% pan(cmd);
% if strcmp(cmd,'on')
%     ea_distogzoomin;
%     ea_distogzoomout;
%     ea_distogrotate;
%     ea_distogslide;
% else
%     ea_view; % reset bars to center.
% end

cameratoolbar(gcf,'SetMode','dollyhv');

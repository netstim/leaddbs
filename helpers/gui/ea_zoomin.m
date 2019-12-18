function ea_zoomin(h,~,cmd)

% if strcmp(cmd,'on')
%     ea_distogpan;
%     ea_distogrotate;
%     ea_distogzoomout;
%     ea_distogslide;
% end
% h=zoom;
% h.Enable=cmd;
% h.Motion='both';
% h.Direction='in';

cameratoolbar(gcf,'SetMode','zoom');

function ea_zoomout(h,~,cmd)

if strcmp(cmd,'on')
    ea_distogpan;
    ea_distogrotate;
    ea_distogzoomin;
    ea_distogslide;
end
h=zoom;
h.Enable=cmd;
h.Motion='both';
h.Direction='out';
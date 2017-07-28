function ea_pan(h,~,cmd)

pan(cmd);
if strcmp(cmd,'on')
    ea_distogzoomin;
    ea_distogzoomout;
    ea_distogrotate;
    ea_distogslide;
end
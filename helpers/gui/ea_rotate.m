function ea_rotate(h,~,cmd)

h = rotate3d;
h.RotateStyle = 'orbit';
h.Enable = cmd;

if strcmp(cmd,'on')
    ea_distogzoomin;
    ea_distogzoomout;
    ea_distogpan;
    ea_distogslide;
end
function [scit,lcit]=ea_getspacedefcit

sd=ea_getspacedef;
if isfield(sd,'citation')
    scit=['(',sd.citation{1},')'];
    lcit=sd.citation{2};
else
    scit=''; lcit='';
end
function ea_lighteffect(resultfig,factor,reset)

if ~exist('factor','var')
    factor=1;
end
if ~exist('reset','var')
    reset=0;
end
ea_alter_reflection(resultfig,factor,reset);

if rand(1)>0.9
   ea_flicker_light(resultfig,factor,reset); 
end



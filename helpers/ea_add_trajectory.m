function ea_add_trajectory(~,~,pobj)

if exist('pobj','var')
    ea_trajectory(pobj);
else
    ea_trajectory;
end

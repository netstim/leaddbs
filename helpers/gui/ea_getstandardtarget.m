function target=ea_getstandardtarget(side)

switch side
    case 1
        target.entry=[27.58,1.59,49.12];
        target.target=[12.58,-13.41,-5.87];
    case 2
        target.entry=[-27.58,1.59,49.13];
        target.target=[-12.58,-13.41,-5.87];
end
target.offset=0;
target.hemisphere=1; % right
end
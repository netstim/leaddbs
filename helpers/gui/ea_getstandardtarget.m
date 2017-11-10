function target=ea_getstandardtarget(site)

switch site
    case 1
        target.entry=[36,39,-50];
        target.target=[12.02,-1.53,1.91];
    case 2
        target.entry=[-36,39,-50];
        target.target=[-12.02,-1.53,1.91];
end
target.offset=0;
target.hemisphere=1; % right
end
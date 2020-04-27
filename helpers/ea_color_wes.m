function pal=ea_color_wes(set)
% source: http://opencolor.tools/palettes/wesanderson/
% similar or same palettes found in R here: https://github.com/karthik/wesanderson
% different sets here https://prafter.com/color/
sets={'lifeaquatic','budapest','moonrisekingdom','royaltenenbaums','fantasticmrfox','darjeeling','hotelchevalier','rushmore'};
if ~exist('set','var')
    s=randperm(length(sets));
    set=sets{s(1)};
end
switch set
    case 'all'
        pal=[];
        for s=1:length(sets)
            pal=[pal;ea_color_wes(sets{s})];
        end
    case 'lifeaquatic'
        pal(1,:)=ea_hex2rgb('#1DACE8');
        pal(2,:)=ea_hex2rgb('#1C366B');
        pal(3,:)=ea_hex2rgb('#F24D29');
        pal(4,:)=ea_hex2rgb('#E5C4A1');
        pal(5,:)=ea_hex2rgb('#C4CFD0');
    case 'royaltenenbaums'
        pal(1,:)=ea_hex2rgb('#9A872D');
        pal(2,:)=ea_hex2rgb('#F5CDB6');
        pal(3,:)=ea_hex2rgb('#F7B0AA');
        pal(4,:)=ea_hex2rgb('#FDDDA4');
        pal(5,:)=ea_hex2rgb('#76A08A');
    case 'budapest'
        pal(1,:)=ea_hex2rgb('#D8A49B');
        pal(2,:)=ea_hex2rgb('#C7CEF6');
        pal(3,:)=ea_hex2rgb('#7496D2');
    case 'moonrisekingdom'
        pal(1,:)=ea_hex2rgb('#B62A3D');
        pal(2,:)=ea_hex2rgb('#EDCB64');
        pal(3,:)=ea_hex2rgb('#B5966D');
        pal(4,:)=ea_hex2rgb('#DAECED');
        pal(5,:)=ea_hex2rgb('#CECD7B');
    case 'fantasticmrfox'
        pal(1,:)=ea_hex2rgb('#F8DF4F');
        pal(2,:)=ea_hex2rgb('#A35E60');
        pal(3,:)=ea_hex2rgb('#541F12');
        pal(4,:)=ea_hex2rgb('#CC8B3C');
        pal(5,:)=ea_hex2rgb('#E8D2B9');
    case 'darjeeling'
        pal(1,:)=ea_hex2rgb('#AEA8A8');
        pal(2,:)=ea_hex2rgb('#CB9E23');
        pal(3,:)=ea_hex2rgb('#957A6D');
        pal(4,:)=ea_hex2rgb('#AC6E49');
    case 'hotelchevalier'
        pal(1,:)=ea_hex2rgb('#456355');
        pal(2,:)=ea_hex2rgb('#FCD16B');
        pal(3,:)=ea_hex2rgb('#D3DDDC');
        pal(4,:)=ea_hex2rgb('#C6B19D');
    case 'rushmore'
        pal(1,:)=ea_hex2rgb('#DBB165');
        pal(2,:)=ea_hex2rgb('#DEB18B');
        pal(3,:)=ea_hex2rgb('#2E604A');
        pal(4,:)=ea_hex2rgb('#27223C');
        pal(5,:)=ea_hex2rgb('#D1362F');
end

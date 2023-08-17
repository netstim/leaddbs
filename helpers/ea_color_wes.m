function colors=ea_color_wes(set, colnum)
% source: http://opencolor.tools/palettes/wesanderson/
% similar or same palettes found in R here: https://github.com/karthik/wesanderson
% different sets here https://prafter.com/color/
sets={'lifeaquatic','budapest','moonrisekingdom','royaltenenbaums','fantasticmrfox','darjeeling','hotelchevalier','rushmore'};
if ~exist('set','var') || isempty(set)
    s=randperm(length(sets));
    set=sets{s(1)};
end

switch set
    case 'all'
        colors=[];
        for s=1:length(sets)
            colors=[colors;ea_color_wes(sets{s})];
        end
    case 'lifeaquatic'
        colors(1,:)=ea_hex2rgb('#1DACE8');
        colors(2,:)=ea_hex2rgb('#1C366B');
        colors(3,:)=ea_hex2rgb('#F24D29');
        colors(4,:)=ea_hex2rgb('#E5C4A1');
        colors(5,:)=ea_hex2rgb('#C4CFD0');
    case 'royaltenenbaums'
        colors(1,:)=ea_hex2rgb('#9A872D');
        colors(2,:)=ea_hex2rgb('#F5CDB6');
        colors(3,:)=ea_hex2rgb('#F7B0AA');
        colors(4,:)=ea_hex2rgb('#FDDDA4');
        colors(5,:)=ea_hex2rgb('#76A08A');
    case 'budapest'
        colors(1,:)=ea_hex2rgb('#D8A49B');
        colors(2,:)=ea_hex2rgb('#C7CEF6');
        colors(3,:)=ea_hex2rgb('#7496D2');
    case 'moonrisekingdom'
        colors(1,:)=ea_hex2rgb('#B62A3D');
        colors(2,:)=ea_hex2rgb('#EDCB64');
        colors(3,:)=ea_hex2rgb('#B5966D');
        colors(4,:)=ea_hex2rgb('#DAECED');
        colors(5,:)=ea_hex2rgb('#CECD7B');
    case 'fantasticmrfox'
        colors(1,:)=ea_hex2rgb('#F8DF4F');
        colors(2,:)=ea_hex2rgb('#A35E60');
        colors(3,:)=ea_hex2rgb('#541F12');
        colors(4,:)=ea_hex2rgb('#CC8B3C');
        colors(5,:)=ea_hex2rgb('#E8D2B9');
    case 'darjeeling'
        colors(1,:)=ea_hex2rgb('#AEA8A8');
        colors(2,:)=ea_hex2rgb('#CB9E23');
        colors(3,:)=ea_hex2rgb('#957A6D');
        colors(4,:)=ea_hex2rgb('#AC6E49');
    case 'hotelchevalier'
        colors(1,:)=ea_hex2rgb('#456355');
        colors(2,:)=ea_hex2rgb('#FCD16B');
        colors(3,:)=ea_hex2rgb('#D3DDDC');
        colors(4,:)=ea_hex2rgb('#C6B19D');
    case 'rushmore'
        colors(1,:)=ea_hex2rgb('#DBB165');
        colors(2,:)=ea_hex2rgb('#DEB18B');
        colors(3,:)=ea_hex2rgb('#2E604A');
        colors(4,:)=ea_hex2rgb('#27223C');
        colors(5,:)=ea_hex2rgb('#D1362F');
end

% Set number of colors
if exist('colnum', 'var')
    if colnum > size(colors,1)
        % Repeat colors in case pre-defined ones are not enough
        colors = repmat(colors, ceil(colnum/size(colors,1)),1);
    end
    colors = colors(1:colnum,:);
end


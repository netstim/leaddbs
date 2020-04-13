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
        pal(1,:)=hex2rgb('#1DACE8');
        pal(2,:)=hex2rgb('#1C366B');
        pal(3,:)=hex2rgb('#F24D29');
        pal(4,:)=hex2rgb('#E5C4A1');
        pal(5,:)=hex2rgb('#C4CFD0');
    case 'royaltenenbaums'
        pal(1,:)=hex2rgb('#9A872D');
        pal(2,:)=hex2rgb('#F5CDB6');
        pal(3,:)=hex2rgb('#F7B0AA');
        pal(4,:)=hex2rgb('#FDDDA4');
        pal(5,:)=hex2rgb('#76A08A');
    case 'budapest'
        pal(1,:)=hex2rgb('#D8A49B');
        pal(2,:)=hex2rgb('#C7CEF6');
        pal(3,:)=hex2rgb('#7496D2');
    case 'moonrisekingdom'
        pal(1,:)=hex2rgb('#B62A3D');
        pal(2,:)=hex2rgb('#EDCB64');
        pal(3,:)=hex2rgb('#B5966D');
        pal(4,:)=hex2rgb('#DAECED');
        pal(5,:)=hex2rgb('#CECD7B');
    case 'fantasticmrfox'
        pal(1,:)=hex2rgb('#F8DF4F');
        pal(2,:)=hex2rgb('#A35E60');
        pal(3,:)=hex2rgb('#541F12');
        pal(4,:)=hex2rgb('#CC8B3C');
        pal(5,:)=hex2rgb('#E8D2B9');
    case 'darjeeling'
        pal(1,:)=hex2rgb('#AEA8A8');
        pal(2,:)=hex2rgb('#CB9E23');
        pal(3,:)=hex2rgb('#957A6D');
        pal(4,:)=hex2rgb('#AC6E49');
    case 'hotelchevalier'
        pal(1,:)=hex2rgb('#456355');
        pal(2,:)=hex2rgb('#FCD16B');
        pal(3,:)=hex2rgb('#D3DDDC');
        pal(4,:)=hex2rgb('#C6B19D');
    case 'rushmore'
        pal(1,:)=hex2rgb('#DBB165');
        pal(2,:)=hex2rgb('#DEB18B');
        pal(3,:)=hex2rgb('#2E604A');
        pal(4,:)=hex2rgb('#27223C');
        pal(5,:)=hex2rgb('#D1362F');
end





function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1.
%
%
% * * * * * * * * * * * * * * * * * * * *
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default.
%
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255.
%
%
% * * * * * * * * * * * * * * * * * * * *
% EXAMPLES:
%
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
%
%
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional
%    = 0.2000    0.3020    0.4000
%
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
%
%
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
%
%
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
%
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'};
% rgbvals = hex2rgb(HexValsAsACharacterArray)
%
% * * * * * * * * * * * * * * * * * * * *
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input.
%
% * * * * * * * * * * * * * * * * * * * *
% See also rgb2hex, dec2hex, hex2num, and ColorSpec.
%
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.')
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary:
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex';
    end
    
    % If input is cell, convert to matrix:
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1;
end
%% Convert from hex to rgb:
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
        
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end


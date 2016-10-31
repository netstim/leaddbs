% green.m
%
%        $Id: green.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = green(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR green
%

function colorMap = green(numberColors)

minValue = .2;
hsvGreen = rgb2hsv(color2RGB('green'));
colorMap = repmat(hsvGreen,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);
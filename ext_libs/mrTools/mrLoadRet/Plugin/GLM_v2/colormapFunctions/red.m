% red.m
%
%        $Id: red.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = red(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR red
%

function colorMap = red(numberColors)

minValue = .2;
hsvRed = rgb2hsv(color2RGB('red'));
colorMap = repmat(hsvRed,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);
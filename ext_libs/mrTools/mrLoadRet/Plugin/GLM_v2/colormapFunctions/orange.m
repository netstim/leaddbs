% orange.m
%
%        $Id: orange.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = orange(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR orange
%

function colorMap = orange(numberColors)

minValue = .2;
hsvOrange = rgb2hsv(color2RGB('orange'));
colorMap = repmat(hsvOrange,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);
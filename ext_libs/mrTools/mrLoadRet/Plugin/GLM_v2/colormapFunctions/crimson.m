% crimson.m
%
%        $Id: crimson.m 1950 2010-12-18 10:12:48Z julien $ 
%      usage: colorMap = crimson(numberColors)
%         by: julien besle
%       date: 13/02/11
%    purpose: returns a colorMap from dark to bright MLR crimson
%

function colorMap = crimson(numberColors)

minValue = .2;
hsvRed = rgb2hsv([1 0 .3]);
colorMap = repmat(hsvRed,numberColors,1);
colorMap(:,3) = (1:numberColors)/numberColors*(1-minValue)+minValue;
colorMap = hsv2rgb(colorMap);